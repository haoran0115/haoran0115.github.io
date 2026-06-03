# haoran0115.github.io

一个极简的个人学术主页骨架，适合直接部署到 `haoran0115.github.io`。

## 文件说明

- `index.md`: 主页内容，平时主要改这个文件
- `_layouts/default.html`: 极简页面布局，内置少量样式和 MathJax
- `_config.yml`: GitHub Pages / Jekyll 的站点配置
- `Gemfile`: 本地预览所需依赖
- `Gemfile.lock`: 锁定本地依赖版本，尽量和 GitHub Pages 保持一致

## 本地生成与预览

这台机器上我已经安装了 Homebrew 的 `ruby@3.3`。如果你在另一台 macOS 机器上从头配置，可以先运行：

```bash
brew install ruby@3.3
echo 'export PATH="/opt/homebrew/opt/ruby@3.3/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

确认版本：

```bash
ruby -v
bundle -v
```

建议看到的是 Ruby 3.3.x。

第一次使用时，在仓库根目录运行：

```bash
bundle config set --local path vendor/bundle
bundle install
```

启动本地预览：

```bash
bundle exec jekyll serve --livereload
```

然后打开：

- `http://127.0.0.1:4000/`

如果只想先生成静态文件，不启动本地服务器，可以运行：

```bash
bundle exec jekyll build
```

生成结果会在 `_site/` 目录。

如果命令找不到新版 Ruby，也可以临时这样运行：

```bash
PATH="/opt/homebrew/opt/ruby@3.3/bin:$PATH" bundle exec jekyll serve --livereload
```

## 日常修改

1. 编辑 `index.md`
2. 把示例里的姓名、单位、邮箱、论文列表替换成你的信息
3. 保存后刷新本地预览页面

## 部署到 GitHub Pages

这个仓库名已经是 `haoran0115.github.io`，所以它会被部署成用户主页站点。

1. 提交并推送代码：

```bash
git add .
git commit -m "Update homepage"
git push origin main
```

2. 打开 GitHub 仓库页面，进入：
   `Settings -> Pages`
3. 在 `Build and deployment` 下选择：
   `Source: Deploy from a branch`
4. Branch 选择 `main`，Folder 选择 `/(root)`，点击 `Save`
5. 等待 GitHub Pages 完成发布

发布后网址通常是：

- `https://haoran0115.github.io/`

如果页面没有立刻更新，通常等待几分钟再刷新即可。GitHub 官方说明里提到，Pages 发布可能需要几分钟。

## 公式示例

行内公式：

```md
$\nabla_\theta \mathcal{L}$
```

块公式：

```md
$$
\min_{\theta} \sum_{i=1}^{n} \ell(f_{\theta}(x_i), y_i)
$$
```
