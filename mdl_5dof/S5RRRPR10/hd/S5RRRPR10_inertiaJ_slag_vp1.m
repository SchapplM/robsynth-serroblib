% Calculate joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:40
% EndTime: 2019-12-31 21:26:48
% DurationCPUTime: 2.75s
% Computational Cost: add. (14498->493), mult. (25571->707), div. (0->0), fcn. (31991->12), ass. (0->238)
t219 = cos(pkin(5));
t224 = sin(qJ(1));
t227 = cos(qJ(2));
t264 = t224 * t227;
t223 = sin(qJ(2));
t228 = cos(qJ(1));
t265 = t223 * t228;
t201 = t219 * t265 + t264;
t250 = qJ(3) + pkin(10);
t216 = sin(t250);
t218 = sin(pkin(5));
t240 = cos(t250);
t234 = t218 * t240;
t171 = t201 * t216 + t228 * t234;
t262 = t227 * t228;
t266 = t223 * t224;
t203 = -t219 * t266 + t262;
t173 = t203 * t216 - t224 * t234;
t272 = t218 * t223;
t188 = t216 * t272 - t219 * t240;
t271 = t218 * t224;
t174 = t203 * t240 + t216 * t271;
t202 = t219 * t264 + t265;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t137 = -t174 * t221 + t202 * t225;
t138 = t174 * t225 + t202 * t221;
t268 = t218 * t228;
t172 = t201 * t240 - t216 * t268;
t200 = -t219 * t262 + t266;
t135 = -t172 * t221 + t200 * t225;
t136 = t172 * t225 + t200 * t221;
t74 = Icges(6,5) * t136 + Icges(6,6) * t135 + Icges(6,3) * t171;
t76 = Icges(6,4) * t136 + Icges(6,2) * t135 + Icges(6,6) * t171;
t78 = Icges(6,1) * t136 + Icges(6,4) * t135 + Icges(6,5) * t171;
t26 = t137 * t76 + t138 * t78 + t173 * t74;
t75 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t173;
t77 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t173;
t79 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t173;
t27 = t137 * t77 + t138 * t79 + t173 * t75;
t189 = t216 * t219 + t223 * t234;
t269 = t218 * t227;
t169 = -t189 * t221 - t225 * t269;
t170 = t189 * t225 - t221 * t269;
t97 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t188;
t98 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t188;
t99 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t188;
t36 = t137 * t98 + t138 * t99 + t173 * t97;
t2 = t171 * t26 + t173 * t27 + t188 * t36;
t281 = t2 / 0.2e1;
t280 = t171 / 0.2e1;
t279 = t173 / 0.2e1;
t278 = t188 / 0.2e1;
t277 = pkin(4) * t172;
t226 = cos(qJ(3));
t215 = pkin(3) * t226 + pkin(2);
t276 = -pkin(2) + t215;
t235 = -rSges(6,1) * t136 - rSges(6,2) * t135;
t80 = rSges(6,3) * t171 - t235;
t275 = pkin(9) * t171 + t277 + t80;
t81 = t138 * rSges(6,1) + t137 * rSges(6,2) + t173 * rSges(6,3);
t274 = t174 * pkin(4) + pkin(9) * t173 + t81;
t220 = -qJ(4) - pkin(8);
t273 = t200 * t220;
t270 = t218 * t226;
t222 = sin(qJ(3));
t267 = t219 * t222;
t263 = t224 * t228;
t100 = rSges(6,1) * t170 + rSges(6,2) * t169 + rSges(6,3) * t188;
t261 = pkin(4) * t189 + pkin(9) * t188 + t100;
t196 = t200 * pkin(8);
t245 = t222 * t268;
t207 = pkin(3) * t245;
t113 = t201 * t276 - t196 - t207 - t273;
t149 = pkin(3) * t267 + ((pkin(8) + t220) * t227 + t276 * t223) * t218;
t260 = t113 * t269 + t200 * t149;
t164 = t203 * pkin(2) + pkin(8) * t202;
t246 = t222 * t271;
t242 = pkin(3) * t246 - t202 * t220 + t203 * t215;
t114 = -t164 + t242;
t162 = t219 * t164;
t259 = t219 * t114 + t162;
t112 = t174 * rSges(5,1) - t173 * rSges(5,2) + t202 * rSges(5,3);
t258 = -t112 - t114;
t163 = pkin(2) * t201 + t196;
t257 = -t113 - t163;
t177 = -t201 * t222 - t226 * t268;
t178 = t201 * t226 - t245;
t122 = rSges(4,1) * t178 + rSges(4,2) * t177 + rSges(4,3) * t200;
t256 = -t122 - t163;
t142 = Icges(5,4) * t189 - Icges(5,2) * t188 - Icges(5,6) * t269;
t143 = Icges(5,1) * t189 - Icges(5,4) * t188 - Icges(5,5) * t269;
t255 = -t188 * t142 + t189 * t143;
t198 = t219 * t226 - t222 * t272;
t199 = t223 * t270 + t267;
t147 = Icges(4,4) * t199 + Icges(4,2) * t198 - Icges(4,6) * t269;
t148 = Icges(4,1) * t199 + Icges(4,4) * t198 - Icges(4,5) * t269;
t254 = t198 * t147 + t199 * t148;
t144 = rSges(5,1) * t189 - rSges(5,2) * t188 - rSges(5,3) * t269;
t253 = -t144 - t149;
t252 = t163 * t271 + t164 * t268;
t251 = t228 * pkin(1) + pkin(7) * t271;
t40 = t169 * t98 + t170 * t99 + t188 * t97;
t32 = t169 * t76 + t170 * t78 + t188 * t74;
t35 = t135 * t98 + t136 * t99 + t171 * t97;
t249 = t32 / 0.2e1 + t35 / 0.2e1;
t33 = t169 * t77 + t170 * t79 + t188 * t75;
t248 = t36 / 0.2e1 + t33 / 0.2e1;
t247 = -t114 - t274;
t244 = -t149 - t261;
t179 = -t203 * t222 + t224 * t270;
t180 = t203 * t226 + t246;
t123 = t180 * rSges(4,1) + t179 * rSges(4,2) + t202 * rSges(4,3);
t185 = Icges(3,3) * t219 + (Icges(3,5) * t223 + Icges(3,6) * t227) * t218;
t186 = Icges(3,6) * t219 + (Icges(3,4) * t223 + Icges(3,2) * t227) * t218;
t187 = Icges(3,5) * t219 + (Icges(3,1) * t223 + Icges(3,4) * t227) * t218;
t243 = t219 * t185 + t186 * t269 + t187 * t272;
t158 = t203 * rSges(3,1) - t202 * rSges(3,2) + rSges(3,3) * t271;
t241 = -t224 * pkin(1) + pkin(7) * t268;
t150 = rSges(4,1) * t199 + rSges(4,2) * t198 - rSges(4,3) * t269;
t204 = (pkin(2) * t223 - pkin(8) * t227) * t218;
t239 = t218 * (-t150 - t204);
t238 = t113 * t271 + t114 * t268 + t252;
t237 = t218 * (-t204 + t253);
t236 = -rSges(5,1) * t172 + rSges(5,2) * t171;
t233 = t242 + t251;
t232 = t218 * (-t204 + t244);
t231 = -t201 * t215 + t207 + t241;
t157 = rSges(3,1) * t201 - rSges(3,2) * t200 - rSges(3,3) * t268;
t105 = Icges(5,5) * t172 - Icges(5,6) * t171 + Icges(5,3) * t200;
t107 = Icges(5,4) * t172 - Icges(5,2) * t171 + Icges(5,6) * t200;
t109 = Icges(5,1) * t172 - Icges(5,4) * t171 + Icges(5,5) * t200;
t50 = -t105 * t269 - t107 * t188 + t109 * t189;
t116 = Icges(4,5) * t178 + Icges(4,6) * t177 + Icges(4,3) * t200;
t118 = Icges(4,4) * t178 + Icges(4,2) * t177 + Icges(4,6) * t200;
t120 = Icges(4,1) * t178 + Icges(4,4) * t177 + Icges(4,5) * t200;
t60 = -t116 * t269 + t118 * t198 + t120 * t199;
t141 = Icges(5,5) * t189 - Icges(5,6) * t188 - Icges(5,3) * t269;
t63 = t141 * t200 - t142 * t171 + t143 * t172;
t146 = Icges(4,5) * t199 + Icges(4,6) * t198 - Icges(4,3) * t269;
t67 = t146 * t200 + t147 * t177 + t148 * t178;
t230 = t67 / 0.2e1 + t60 / 0.2e1 + t50 / 0.2e1 + t63 / 0.2e1 + t249;
t106 = Icges(5,5) * t174 - Icges(5,6) * t173 + Icges(5,3) * t202;
t108 = Icges(5,4) * t174 - Icges(5,2) * t173 + Icges(5,6) * t202;
t110 = Icges(5,1) * t174 - Icges(5,4) * t173 + Icges(5,5) * t202;
t51 = -t106 * t269 - t108 * t188 + t110 * t189;
t117 = Icges(4,5) * t180 + Icges(4,6) * t179 + Icges(4,3) * t202;
t119 = Icges(4,4) * t180 + Icges(4,2) * t179 + Icges(4,6) * t202;
t121 = Icges(4,1) * t180 + Icges(4,4) * t179 + Icges(4,5) * t202;
t61 = -t117 * t269 + t119 * t198 + t121 * t199;
t64 = t141 * t202 - t142 * t173 + t143 * t174;
t68 = t146 * t202 + t147 * t179 + t148 * t180;
t229 = t61 / 0.2e1 + t51 / 0.2e1 + t68 / 0.2e1 + t64 / 0.2e1 + t248;
t209 = rSges(2,1) * t228 - rSges(2,2) * t224;
t208 = -rSges(2,1) * t224 - rSges(2,2) * t228;
t190 = t219 * rSges(3,3) + (rSges(3,1) * t223 + rSges(3,2) * t227) * t218;
t156 = Icges(3,1) * t203 - Icges(3,4) * t202 + Icges(3,5) * t271;
t155 = Icges(3,1) * t201 - Icges(3,4) * t200 - Icges(3,5) * t268;
t154 = Icges(3,4) * t203 - Icges(3,2) * t202 + Icges(3,6) * t271;
t153 = Icges(3,4) * t201 - Icges(3,2) * t200 - Icges(3,6) * t268;
t152 = Icges(3,5) * t203 - Icges(3,6) * t202 + Icges(3,3) * t271;
t151 = Icges(3,5) * t201 - Icges(3,6) * t200 - Icges(3,3) * t268;
t140 = t158 + t251;
t139 = -t157 + t241;
t127 = -t157 * t219 - t190 * t268;
t126 = t158 * t219 - t190 * t271;
t115 = t243 * t219;
t111 = rSges(5,3) * t200 - t236;
t96 = (t157 * t224 + t158 * t228) * t218;
t95 = t185 * t271 - t186 * t202 + t187 * t203;
t94 = -t185 * t268 - t186 * t200 + t187 * t201;
t93 = t202 * t113;
t91 = t164 + t123 + t251;
t90 = t241 + t256;
t87 = -t123 * t269 - t150 * t202;
t86 = t122 * t269 + t150 * t200;
t85 = t219 * t152 + (t154 * t227 + t156 * t223) * t218;
t84 = t219 * t151 + (t153 * t227 + t155 * t223) * t218;
t83 = t233 + t112;
t82 = (-rSges(5,3) + t220) * t200 + t231 + t236;
t73 = -t146 * t269 + t254;
t72 = t73 * t219;
t71 = t122 * t202 - t123 * t200;
t70 = t219 * t256 + t228 * t239;
t69 = t123 * t219 + t224 * t239 + t162;
t66 = -t141 * t269 + t255;
t65 = t66 * t219;
t62 = (t122 * t224 + t123 * t228) * t218 + t252;
t59 = t233 + t274;
t58 = -t277 + t273 + (-rSges(6,3) - pkin(9)) * t171 + t231 + t235;
t57 = -t100 * t173 + t188 * t81;
t56 = t100 * t171 - t188 * t80;
t55 = t117 * t202 + t119 * t179 + t121 * t180;
t54 = t116 * t202 + t118 * t179 + t120 * t180;
t53 = t117 * t200 + t119 * t177 + t121 * t178;
t52 = t116 * t200 + t118 * t177 + t120 * t178;
t49 = t202 * t253 + t258 * t269;
t48 = t111 * t269 + t144 * t200 + t260;
t47 = t106 * t202 - t108 * t173 + t110 * t174;
t46 = t105 * t202 - t107 * t173 + t109 * t174;
t45 = t106 * t200 - t108 * t171 + t110 * t172;
t44 = t105 * t200 - t107 * t171 + t109 * t172;
t43 = (-t111 + t257) * t219 + t228 * t237;
t42 = t112 * t219 + t224 * t237 + t259;
t41 = -t171 * t81 + t173 * t80;
t39 = t40 * t219;
t38 = t40 * t188;
t37 = t111 * t202 + t200 * t258 + t93;
t34 = (t111 * t224 + t112 * t228) * t218 + t238;
t31 = t202 * t244 + t247 * t269;
t30 = t200 * t261 + t269 * t275 + t260;
t29 = (t257 - t275) * t219 + t228 * t232;
t28 = t219 * t274 + t224 * t232 + t259;
t25 = t135 * t77 + t136 * t79 + t171 * t75;
t24 = t135 * t76 + t136 * t78 + t171 * t74;
t23 = t200 * t247 + t202 * t275 + t93;
t22 = (t224 * t275 + t228 * t274) * t218 + t238;
t21 = t72 + (t61 * t224 - t60 * t228) * t218;
t20 = t60 * t200 + t61 * t202 - t269 * t73;
t19 = t68 * t219 + (t224 * t55 - t228 * t54) * t218;
t18 = t67 * t219 + (t224 * t53 - t228 * t52) * t218;
t17 = t65 + (t51 * t224 - t50 * t228) * t218;
t16 = t200 * t54 + t202 * t55 - t269 * t68;
t15 = t200 * t52 + t202 * t53 - t269 * t67;
t14 = t50 * t200 + t51 * t202 - t269 * t66;
t13 = t64 * t219 + (t224 * t47 - t228 * t46) * t218;
t12 = t63 * t219 + (t224 * t45 - t228 * t44) * t218;
t11 = t200 * t46 + t202 * t47 - t269 * t64;
t10 = t200 * t44 + t202 * t45 - t269 * t63;
t9 = t39 + (t33 * t224 - t32 * t228) * t218;
t8 = t32 * t200 + t33 * t202 - t269 * t40;
t7 = t32 * t171 + t33 * t173 + t38;
t6 = t36 * t219 + (t224 * t27 - t228 * t26) * t218;
t5 = t35 * t219 + (t224 * t25 - t228 * t24) * t218;
t4 = t200 * t26 + t202 * t27 - t269 * t36;
t3 = t200 * t24 + t202 * t25 - t269 * t35;
t1 = t171 * t24 + t173 * t25 + t188 * t35;
t88 = [Icges(2,3) + (-t141 - t146) * t269 + m(6) * (t58 ^ 2 + t59 ^ 2) + m(5) * (t82 ^ 2 + t83 ^ 2) + m(4) * (t90 ^ 2 + t91 ^ 2) + m(3) * (t139 ^ 2 + t140 ^ 2) + m(2) * (t208 ^ 2 + t209 ^ 2) + t243 + t40 + t254 + t255; t39 + t72 + t65 + t115 + m(6) * (t28 * t59 + t29 * t58) + m(5) * (t42 * t83 + t43 * t82) + m(4) * (t69 * t91 + t70 * t90) + m(3) * (t126 * t140 + t127 * t139) + ((-t84 / 0.2e1 - t94 / 0.2e1 - t230) * t228 + (t95 / 0.2e1 + t85 / 0.2e1 + t229) * t224) * t218; (t9 + t21 + t17 + t115) * t219 + m(6) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t34 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(4) * (t62 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2 + t96 ^ 2) + (-t228 * t5 + t224 * t6 - t228 * t18 - t228 * t12 + t224 * t19 + t224 * t13 + (-t228 * ((-t154 * t200 + t156 * t201) * t224 - (-t153 * t200 + t155 * t201) * t228) + t224 * ((-t154 * t202 + t156 * t203) * t224 - (-t153 * t202 + t155 * t203) * t228) + (-t228 * (t151 * t228 ^ 2 - t152 * t263) + t224 * (t152 * t224 ^ 2 - t151 * t263)) * t218) * t218 + ((-t84 - t94) * t228 + (t85 + t95) * t224) * t219) * t218; (-t40 - t66 - t73) * t269 + m(6) * (t30 * t58 + t31 * t59) + m(5) * (t48 * t82 + t49 * t83) + m(4) * (t86 * t90 + t87 * t91) + t229 * t202 + t230 * t200; (t8 / 0.2e1 + t20 / 0.2e1 + t14 / 0.2e1) * t219 + (t6 / 0.2e1 + t19 / 0.2e1 + t13 / 0.2e1) * t202 + (t5 / 0.2e1 + t18 / 0.2e1 + t12 / 0.2e1) * t200 + m(6) * (t22 * t23 + t28 * t31 + t29 * t30) + m(5) * (t34 * t37 + t42 * t49 + t43 * t48) + m(4) * (t62 * t71 + t69 * t87 + t70 * t86) + ((-t3 / 0.2e1 - t10 / 0.2e1 - t15 / 0.2e1) * t228 + (-t9 / 0.2e1 - t17 / 0.2e1 - t21 / 0.2e1) * t227 + (t4 / 0.2e1 + t16 / 0.2e1 + t11 / 0.2e1) * t224) * t218; (-t14 - t20 - t8) * t269 + (t4 + t11 + t16) * t202 + (t3 + t10 + t15) * t200 + m(6) * (t23 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t37 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t71 ^ 2 + t86 ^ 2 + t87 ^ 2); m(6) * (t200 * t59 + t202 * t58) + m(5) * (t200 * t83 + t202 * t82); m(6) * (t200 * t28 + t202 * t29 - t22 * t269) + m(5) * (t200 * t42 + t202 * t43 - t269 * t34); m(6) * (t200 * t31 + t202 * t30 - t23 * t269) + m(5) * (t200 * t49 + t202 * t48 - t269 * t37); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t218 ^ 2 * t227 ^ 2 + t200 ^ 2 + t202 ^ 2); m(6) * (t56 * t58 + t57 * t59) + t38 + t248 * t173 + t249 * t171; m(6) * (t22 * t41 + t28 * t57 + t29 * t56) + t5 * t280 + t9 * t278 + t6 * t279 + t219 * t7 / 0.2e1 + (t224 * t281 - t228 * t1 / 0.2e1) * t218; t200 * t1 / 0.2e1 + t4 * t279 + m(6) * (t23 * t41 + t30 * t56 + t31 * t57) + t3 * t280 + t8 * t278 + t202 * t281 - t7 * t269 / 0.2e1; m(6) * (t200 * t57 + t202 * t56 - t269 * t41); m(6) * (t41 ^ 2 + t56 ^ 2 + t57 ^ 2) + t173 * t2 + t171 * t1 + t188 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t88(1), t88(2), t88(4), t88(7), t88(11); t88(2), t88(3), t88(5), t88(8), t88(12); t88(4), t88(5), t88(6), t88(9), t88(13); t88(7), t88(8), t88(9), t88(10), t88(14); t88(11), t88(12), t88(13), t88(14), t88(15);];
Mq = res;
