## scipy 编译 ios arm64

Fortran代码问题
采用F2C 转成C文件
http://www.netlib.org/f2c/
转换后的C代码运行有些问题。

大部分C模块都可以运行，只有F2C模块有问题，后期慢慢修整。
假如您有Fortran在IOS编译的方法提供，我将非常感激您的贡献。

当前scipy版本 1.2.2
* main.m 模块导出函数
* ScipyImporter.py 模块导入包装方法
* scipy目录 脚本文件版本 1.2.2
* scipylib目录  编译IOS库，ARM64格式
* scipy_build   编译源码，只提供了F2C文件，其它C模块只提供o文件，需要可下载源码编译

APPStore上线测试应用
* CN  https://apps.apple.com/cn/app/id1471351733
* US  https://apps.apple.com/us/app/id1471351733

使用集成方法参考
* https://github.com/ColdGrub1384/Pyto

部分F2C模块不能运行，否则会崩溃
但是scipy依赖库，貌似不影响其它库运行
比如：scikit-image、scikit-learn、statsmodels 在测试应用中，我做了一些测试，可能正常运行。

我把方法公开出来，希望大家可以共同完美scipy在IOS上正确运行。
为手机上实现AI入门基础提供帮助，推动人类人工智能的实现一大步。
