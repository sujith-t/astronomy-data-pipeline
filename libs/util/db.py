#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May  13 18:30:21 2026

@author: sujith-t
"""

import mysql.connector

class MySQLUtil:

    connection = None

    def __init__(self, op_system):
        self.query = list()
        self.params = list()

        if MySQLUtil.connection is None:
            MySQLUtil.connection = mysql.connector.connect(host=op_system.getenv("DB_HOST"), user=op_system.getenv("DB_USER"), password=op_system.getenv("DB_PASSWORD"), database=op_system.getenv("DB_NAME"))

    def execute(self, sql:str, params=None, commit=False):
        self.query.append(sql)
        self.params.append(params)

        if commit:
            cursor = MySQLUtil.connection.cursor()
            for q, p in zip(self.query, self.params):
                if p is None:
                    cursor.execute(q)
                else:
                    cursor.execute(q, p)

            MySQLUtil.connection.commit()
            self.query.clear()
            self.params.clear()
            cursor.close()

    def fetch_all(self, sql:str, params=None, col_names:bool=False):
        cursor = MySQLUtil.connection.cursor(dictionary=col_names)
        if params is None:
            cursor.execute(sql)
        else:
            cursor.execute(sql, params)

        rows = cursor.fetchall()
        cursor.close()
        return rows

    def fetch_one(self, sql:str, params=None, col_names:bool=False):
        cursor = MySQLUtil.connection.cursor(dictionary=col_names)
        if params is None:
            cursor.execute(sql)
        else:
            cursor.execute(sql, params)

        row = cursor.fetchone()
        cursor.close()
        return row

    def close(self):
        if MySQLUtil.connection is not None:
            MySQLUtil.connection.close()
            MySQLUtil.connection = None
