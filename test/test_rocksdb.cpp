// Copyright (c) 2013, Facebook, Inc.  All rights reserved.
// This source code is licensed under the BSD-style license found in the
// LICENSE file in the root directory of this source tree. An additional grant
// of patent rights can be found in the PATENTS file in the same directory.
#include <cstdio>
#include <string>
#include <iostream>

#include "rocksdb/db.h"
#include "rocksdb/env.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

using namespace rocksdb;

std::string kDBPath = "/home/hannese/user/git/gyper/build/test.db";

int main() {
  DB* db;
  Options options;

  options.IncreaseParallelism();
  options.OptimizeLevelStyleCompaction();
  std::cout << options.env << std::endl;
  // Env* default_env = Env::Default();q
  options.env = rocksdb::NewMemEnv(options.env);
  // Optimize RocksDB. This is the easiest way to get RocksDB to perform well
  
  // create the DB if it's not already present
  options.create_if_missing = true;

  Status s = options.env->CreateDir("/dir");
  if (!s.ok())
    std::cerr << s.ToString() << std::endl;

  std::vector<std::string> children;
  s = options.env->GetChildren("/dir", &children);
  if (!s.ok())
    std::cerr << s.ToString() << std::endl;

  std::cout << "children.size() = " << children.size() << std::endl;
  // open DB
  s = DB::Open(options, "/dir/db", &db);
  if (!s.ok())
  {
    std::cerr << s.ToString() << std::endl;
    return 1;
  }
  
  // Put key-value
  s = db->Put(WriteOptions(), "key1", "value");
  assert(s.ok());
  std::string value;
  // get value
  s = db->Get(ReadOptions(), "key1", &value);
  assert(s.ok());
  assert(value == "value");

  // atomically apply a set of updates
  {
    WriteBatch batch;
    batch.Delete("key1");
    batch.Put("key2", value);
    s = db->Write(WriteOptions(), &batch);
  }

  s = db->Get(ReadOptions(), "key1", &value);
  assert(s.IsNotFound());

  db->Get(ReadOptions(), "key2", &value);
  assert(value == "value");

  std::cout << value << std::endl;

  delete db;

  return 0;
}
