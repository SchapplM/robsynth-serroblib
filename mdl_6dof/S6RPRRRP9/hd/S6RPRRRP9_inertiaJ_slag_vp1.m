% Calculate joint inertia matrix for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:10
% EndTime: 2019-03-09 06:26:20
% DurationCPUTime: 4.27s
% Computational Cost: add. (8090->432), mult. (11445->600), div. (0->0), fcn. (12280->8), ass. (0->214)
t218 = qJ(4) + qJ(5);
t206 = sin(t218);
t207 = cos(t218);
t224 = cos(qJ(1));
t220 = sin(qJ(3));
t221 = sin(qJ(1));
t275 = t220 * t221;
t167 = -t206 * t275 + t207 * t224;
t168 = t206 * t224 + t207 * t275;
t223 = cos(qJ(3));
t272 = t221 * t223;
t95 = Icges(7,5) * t168 + Icges(7,6) * t167 - Icges(7,3) * t272;
t97 = Icges(6,5) * t168 + Icges(6,6) * t167 - Icges(6,3) * t272;
t323 = t95 + t97;
t274 = t220 * t224;
t169 = t206 * t274 + t207 * t221;
t170 = t206 * t221 - t207 * t274;
t270 = t223 * t224;
t96 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t270;
t98 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t270;
t322 = t96 + t98;
t101 = Icges(6,4) * t168 + Icges(6,2) * t167 - Icges(6,6) * t272;
t99 = Icges(7,4) * t168 + Icges(7,2) * t167 - Icges(7,6) * t272;
t321 = t101 + t99;
t100 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t270;
t102 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t270;
t320 = t100 + t102;
t103 = Icges(7,1) * t168 + Icges(7,4) * t167 - Icges(7,5) * t272;
t105 = Icges(6,1) * t168 + Icges(6,4) * t167 - Icges(6,5) * t272;
t319 = t103 + t105;
t104 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t270;
t106 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t270;
t318 = t104 + t106;
t225 = -pkin(9) - pkin(8);
t214 = -qJ(6) + t225;
t311 = -rSges(7,3) + t214;
t317 = t321 * t167 + t319 * t168 - t272 * t323;
t316 = t167 * t320 + t168 * t318 - t272 * t322;
t315 = t321 * t169 + t319 * t170 + t270 * t323;
t314 = t169 * t320 + t170 * t318 + t270 * t322;
t144 = Icges(7,3) * t220 + (Icges(7,5) * t207 - Icges(7,6) * t206) * t223;
t146 = Icges(7,6) * t220 + (Icges(7,4) * t207 - Icges(7,2) * t206) * t223;
t148 = Icges(7,5) * t220 + (Icges(7,1) * t207 - Icges(7,4) * t206) * t223;
t65 = -t144 * t272 + t146 * t167 + t148 * t168;
t145 = Icges(6,3) * t220 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t223;
t147 = Icges(6,6) * t220 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t223;
t149 = Icges(6,5) * t220 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t223;
t66 = -t145 * t272 + t147 * t167 + t149 * t168;
t313 = t65 + t66;
t67 = t144 * t270 + t146 * t169 + t148 * t170;
t68 = t145 * t270 + t147 * t169 + t149 * t170;
t312 = t68 + t67;
t310 = Icges(4,5) * t223;
t309 = Icges(4,6) * t220;
t308 = (t148 + t149) * t207 * t223 + (t144 + t145) * t220;
t307 = (-t146 - t147) * t206;
t222 = cos(qJ(4));
t205 = t222 * pkin(4) + pkin(3);
t184 = pkin(5) * t207 + t205;
t219 = sin(qJ(4));
t286 = pkin(4) * t219;
t187 = pkin(5) * t206 + t286;
t306 = t168 * rSges(7,1) + t167 * rSges(7,2) + t184 * t275 + t224 * t187 + t311 * t272;
t305 = t310 / 0.2e1 - t309 / 0.2e1;
t304 = (-t317 * t221 + t316 * t224) * t223 + t313 * t220;
t303 = (-t315 * t221 + t314 * t224) * t223 + t312 * t220;
t302 = t316 * t221 + t317 * t224;
t301 = t314 * t221 + t315 * t224;
t48 = t220 * t95 + (t103 * t207 - t206 * t99) * t223;
t50 = t220 * t97 + (-t101 * t206 + t105 * t207) * t223;
t300 = t48 + t50;
t49 = t220 * t96 + (-t100 * t206 + t104 * t207) * t223;
t51 = t220 * t98 + (-t102 * t206 + t106 * t207) * t223;
t299 = t49 + t51;
t298 = (t223 * t307 + t308) * t220;
t269 = t223 * t225;
t276 = t219 * t224;
t248 = pkin(4) * t276 + t205 * t275 + t221 * t269;
t297 = -t248 + t306;
t263 = (rSges(7,1) * t207 - rSges(7,2) * t206 + t184 - t205) * t223 + (t225 - t311) * t220;
t296 = (rSges(4,1) * t220 + rSges(4,2) * t223) * t224;
t215 = t221 ^ 2;
t217 = t224 ^ 2;
t295 = -pkin(1) - pkin(7);
t292 = t220 / 0.2e1;
t291 = t221 / 0.2e1;
t289 = t224 / 0.2e1;
t194 = rSges(4,1) * t223 - rSges(4,2) * t220;
t288 = m(4) * t194;
t287 = m(7) * t223;
t238 = -rSges(7,1) * t170 - rSges(7,2) * t169;
t256 = t205 * t274 + t224 * t269;
t281 = t184 * t220;
t284 = rSges(7,3) * t270 - t238 + (-t214 * t223 - t281) * t224 + (t187 - t286) * t221 + t256;
t108 = t168 * rSges(6,1) + t167 * rSges(6,2) - rSges(6,3) * t272;
t151 = rSges(6,3) * t220 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t223;
t81 = t220 * t108 + t151 * t272;
t160 = Icges(5,6) * t220 + (Icges(5,4) * t222 - Icges(5,2) * t219) * t223;
t278 = t219 * t160;
t277 = t219 * t221;
t273 = t221 * t222;
t271 = t222 * t224;
t202 = pkin(3) * t275;
t181 = -pkin(8) * t272 + t202;
t124 = -t181 + t248;
t268 = -t108 - t124;
t110 = rSges(6,1) * t170 + rSges(6,2) * t169 + rSges(6,3) * t270;
t204 = pkin(3) * t274;
t241 = pkin(8) * t270 - t204;
t125 = pkin(4) * t277 - t241 - t256;
t267 = -t110 - t125;
t143 = (-pkin(3) + t205) * t223 + (-pkin(8) - t225) * t220;
t266 = t220 * t124 + t143 * t272;
t174 = t224 * t241;
t265 = t224 * t125 + t174;
t264 = t263 * t270;
t196 = t223 * pkin(3) + t220 * pkin(8);
t183 = t221 * t196;
t260 = t221 * t143 + t183;
t157 = Icges(5,3) * t220 + (Icges(5,5) * t222 - Icges(5,6) * t219) * t223;
t163 = Icges(5,5) * t220 + (Icges(5,1) * t222 - Icges(5,4) * t219) * t223;
t259 = t223 * t222 * t163 + t220 * t157;
t258 = -t143 - t196;
t177 = -t219 * t275 + t271;
t178 = t220 * t273 + t276;
t257 = t178 * rSges(5,1) + t177 * rSges(5,2);
t255 = t224 * pkin(1) + t221 * qJ(2);
t254 = t215 + t217;
t116 = Icges(5,5) * t178 + Icges(5,6) * t177 - Icges(5,3) * t272;
t118 = Icges(5,4) * t178 + Icges(5,2) * t177 - Icges(5,6) * t272;
t120 = Icges(5,1) * t178 + Icges(5,4) * t177 - Icges(5,5) * t272;
t58 = t116 * t220 + (-t118 * t219 + t120 * t222) * t223;
t74 = -t157 * t272 + t160 * t177 + t163 * t178;
t253 = -t58 / 0.2e1 - t74 / 0.2e1;
t179 = t219 * t274 + t273;
t180 = -t220 * t271 + t277;
t117 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t270;
t119 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t270;
t121 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t270;
t59 = t117 * t220 + (-t119 * t219 + t121 * t222) * t223;
t75 = t157 * t270 + t160 * t179 + t163 * t180;
t252 = t59 / 0.2e1 + t75 / 0.2e1;
t251 = -t124 - t297;
t250 = -t125 - t284;
t247 = rSges(4,1) * t275 + rSges(4,2) * t272 + t224 * rSges(4,3);
t246 = t224 * pkin(7) + t255;
t245 = (-rSges(5,3) - pkin(8)) * t223;
t244 = -t272 / 0.2e1;
t243 = t270 / 0.2e1;
t242 = t303 * t270 + ((-t221 * t300 + t224 * t299) * t223 + t298) * t220;
t34 = t220 * t297 + t263 * t272;
t239 = -rSges(5,1) * t180 - rSges(5,2) * t179;
t29 = t34 + t266;
t134 = t143 * t270;
t30 = t220 * t250 + t134 + t264;
t237 = t221 * t30 - t224 * t29;
t35 = -t220 * t284 + t264;
t236 = t221 * t35 - t224 * t34;
t70 = t221 * t263 + t260;
t71 = (t258 - t263) * t224;
t235 = t221 * t70 - t224 * t71;
t209 = t224 * qJ(2);
t72 = t209 + (t311 * t223 + t281) * t224 + (-t187 + t295) * t221 + t238;
t73 = t246 + t306;
t234 = t221 * t72 - t224 * t73;
t231 = Icges(4,5) * t220 + Icges(4,6) * t223;
t228 = -t272 * t304 + t242;
t227 = (t221 * t299 + t224 * t300) * t292 + t303 * t291 + t304 * t289 + t302 * t244 + t301 * t243;
t226 = (t300 + t313) * t244 + (t299 + t312) * t243 + t298;
t195 = rSges(2,1) * t224 - rSges(2,2) * t221;
t193 = -rSges(2,1) * t221 - rSges(2,2) * t224;
t189 = -t309 + t310;
t172 = -rSges(3,2) * t224 + rSges(3,3) * t221 + t255;
t171 = rSges(3,3) * t224 + t209 + (rSges(3,2) - pkin(1)) * t221;
t166 = rSges(5,3) * t220 + (rSges(5,1) * t222 - rSges(5,2) * t219) * t223;
t159 = Icges(4,3) * t221 - t224 * t231;
t158 = Icges(4,3) * t224 + t221 * t231;
t138 = t151 * t270;
t130 = t246 + t247;
t129 = t209 + t296 + (-rSges(4,3) + t295) * t221;
t128 = (-t166 - t196) * t224;
t127 = t166 * t221 + t183;
t123 = rSges(5,3) * t270 - t239;
t122 = -rSges(5,3) * t272 + t257;
t111 = -t221 * t247 + (t221 * rSges(4,3) - t296) * t224;
t92 = (-t151 + t258) * t224;
t91 = t151 * t221 + t260;
t90 = t221 * t245 + t202 + t246 + t257;
t89 = t221 * t295 + t224 * t245 + t204 + t209 + t239;
t86 = -t123 * t220 + t166 * t270;
t85 = t122 * t220 + t166 * t272;
t83 = (-t223 * t278 + t259) * t220;
t82 = -t110 * t220 + t138;
t80 = t108 + t246 + t248;
t79 = t209 + (-t286 + t295) * t221 - t110 + t256;
t76 = (-t122 * t224 - t123 * t221) * t223;
t69 = (-t108 * t224 - t110 * t221) * t223;
t60 = t123 * t224 + t174 + (-t122 - t181) * t221;
t57 = t220 * t267 + t134 + t138;
t56 = t266 + t81;
t55 = t117 * t270 + t119 * t179 + t121 * t180;
t54 = t116 * t270 + t118 * t179 + t120 * t180;
t53 = -t117 * t272 + t119 * t177 + t121 * t178;
t52 = -t116 * t272 + t118 * t177 + t120 * t178;
t33 = (t221 * t267 + t224 * t268) * t223;
t32 = t110 * t224 + (-t181 + t268) * t221 + t265;
t31 = (-t221 * t284 - t224 * t297) * t223;
t28 = t221 * t55 + t224 * t54;
t27 = t221 * t53 + t224 * t52;
t24 = (t221 * t250 + t224 * t251) * t223;
t23 = t284 * t224 + (-t181 + t251) * t221 + t265;
t14 = t220 * t75 + (-t221 * t54 + t224 * t55) * t223;
t13 = t220 * t74 + (-t221 * t52 + t224 * t53) * t223;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t223 - t278 + t307) * t223 + m(7) * (t72 ^ 2 + t73 ^ 2) + m(6) * (t79 ^ 2 + t80 ^ 2) + m(5) * (t89 ^ 2 + t90 ^ 2) + m(4) * (t129 ^ 2 + t130 ^ 2) + m(3) * (t171 ^ 2 + t172 ^ 2) + m(2) * (t193 ^ 2 + t195 ^ 2) + t259 + (-0.2e1 * Icges(4,4) * t223 + Icges(4,2) * t220) * t220 + t308; m(7) * t234 + m(6) * (t221 * t79 - t224 * t80) + m(5) * (t221 * t89 - t224 * t90) + m(4) * (t129 * t221 - t130 * t224) + m(3) * (t171 * t221 - t172 * t224); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t254; m(7) * (t70 * t72 + t71 * t73) + m(6) * (t79 * t91 + t80 * t92) + m(5) * (t127 * t89 + t128 * t90) + (t50 / 0.2e1 + t48 / 0.2e1 + t66 / 0.2e1 + t65 / 0.2e1 - t130 * t288 + t189 * t289 - t253 + t305 * t224) * t224 + (t51 / 0.2e1 + t49 / 0.2e1 + t68 / 0.2e1 + t67 / 0.2e1 + t129 * t288 + t189 * t291 + t252 + t305 * t221) * t221; m(5) * (t127 * t221 - t128 * t224) + m(6) * (t221 * t91 - t224 * t92) + m(7) * t235 + t254 * t288; m(7) * (t23 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t32 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t127 ^ 2 + t128 ^ 2 + t60 ^ 2) + m(4) * (t194 ^ 2 * t254 + t111 ^ 2) + (t217 * t158 + t27 + t302) * t224 + (t215 * t159 + t28 + (t221 * t158 + t224 * t159) * t224 + t301) * t221; (t221 * t253 + t224 * t252) * t223 + m(7) * (t29 * t73 + t30 * t72) + m(6) * (t56 * t80 + t57 * t79) + m(5) * (t85 * t90 + t86 * t89) + t83 + t226; m(5) * (t221 * t86 - t224 * t85) + m(6) * (t221 * t57 - t224 * t56) + m(7) * t237; m(7) * (t23 * t24 + t29 * t71 + t30 * t70) + m(6) * (t33 * t32 + t56 * t92 + t57 * t91) + m(5) * (t127 * t86 + t128 * t85 + t60 * t76) + (t28 * t289 - t221 * t27 / 0.2e1) * t223 + (t59 * t221 + t58 * t224) * t292 + t13 * t289 + t14 * t291 + t227; t220 * t83 + m(7) * (t24 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t33 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t76 ^ 2 + t85 ^ 2 + t86 ^ 2) + ((t220 * t59 + t14) * t224 + (-t220 * t58 - t13 - t304) * t221) * t223 + t242; m(7) * (t34 * t73 + t35 * t72) + m(6) * (t79 * t82 + t80 * t81) + t226; m(6) * (t221 * t82 - t224 * t81) + m(7) * t236; m(7) * (t23 * t31 + t34 * t71 + t35 * t70) + m(6) * (t69 * t32 + t81 * t92 + t82 * t91) + t227; m(7) * (t24 * t31 + t29 * t34 + t30 * t35) + m(6) * (t69 * t33 + t56 * t81 + t57 * t82) + t228; m(7) * (t31 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t69 ^ 2 + t81 ^ 2 + t82 ^ 2) + t228; -t234 * t287; -t254 * t287; m(7) * (t220 * t23 - t223 * t235); m(7) * (t220 * t24 - t223 * t237); m(7) * (t220 * t31 - t223 * t236); m(7) * (t223 ^ 2 * t254 + t220 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
