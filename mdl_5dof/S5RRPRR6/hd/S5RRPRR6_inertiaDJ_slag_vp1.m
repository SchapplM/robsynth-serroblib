% Calculate time derivative of joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:41
% DurationCPUTime: 5.22s
% Computational Cost: add. (16061->524), mult. (15018->740), div. (0->0), fcn. (14002->10), ass. (0->248)
t238 = cos(qJ(4));
t221 = pkin(4) * t238 + pkin(3);
t233 = qJ(1) + qJ(2);
t226 = cos(t233);
t235 = cos(pkin(9));
t301 = t226 * t235;
t224 = sin(t233);
t236 = sin(qJ(4));
t304 = t224 * t236;
t323 = pkin(4) * t304 + t221 * t301;
t231 = qJD(1) + qJD(2);
t234 = sin(pkin(9));
t240 = -pkin(8) - pkin(7);
t296 = t234 * t240;
t271 = t231 * t296;
t281 = qJD(4) * t238;
t322 = -pkin(4) * t281 - t271;
t321 = -t224 * rSges(3,1) - t226 * rSges(3,2);
t320 = 2 * m(3);
t319 = 2 * m(4);
t318 = 2 * m(5);
t317 = 2 * m(6);
t229 = t234 ^ 2;
t316 = -pkin(3) + t221;
t315 = pkin(1) * qJD(1);
t232 = qJ(4) + qJ(5);
t225 = cos(t232);
t230 = qJD(4) + qJD(5);
t299 = t230 * t234;
t223 = sin(t232);
t310 = Icges(6,4) * t223;
t146 = (-Icges(6,2) * t225 - t310) * t299;
t157 = -Icges(6,5) * t235 + (Icges(6,1) * t225 - t310) * t234;
t309 = Icges(6,4) * t225;
t156 = -Icges(6,6) * t235 + (-Icges(6,2) * t223 + t309) * t234;
t272 = t230 * t225 * t156;
t244 = -t272 + (-t157 * t230 - t146) * t223;
t147 = (-Icges(6,1) * t223 - t309) * t299;
t144 = t234 * t225 * t147;
t145 = (-Icges(6,5) * t223 - Icges(6,6) * t225) * t299;
t266 = -t235 * t145 + t144;
t53 = t234 * t244 + t266;
t314 = t53 * t235;
t158 = -rSges(6,3) * t235 + (rSges(6,1) * t225 - rSges(6,2) * t223) * t234;
t298 = t231 * t234;
t275 = t224 * t298;
t263 = t230 * t235 - t231;
t252 = t263 * t226;
t297 = t231 * t235;
t264 = -t230 + t297;
t255 = t223 * t264;
t119 = -t224 * t255 + t225 * t252;
t254 = t264 * t225;
t120 = t223 * t252 + t224 * t254;
t256 = -rSges(6,1) * t120 - rSges(6,2) * t119;
t70 = rSges(6,3) * t275 - t256;
t313 = t158 * t275 + t235 * t70;
t312 = Icges(5,4) * t236;
t311 = Icges(5,4) * t238;
t148 = (-rSges(6,1) * t223 - rSges(6,2) * t225) * t299;
t308 = t148 * t224;
t307 = t224 * t231;
t306 = t224 * t234;
t305 = t224 * t235;
t303 = t226 * t231;
t302 = t226 * t234;
t300 = t226 * t236;
t295 = t235 * t236;
t294 = t235 * t238;
t283 = qJD(4) * t234;
t176 = (-Icges(5,2) * t238 - t312) * t283;
t293 = t236 * t176;
t159 = -t223 * t305 - t225 * t226;
t160 = -t223 * t226 + t225 * t305;
t117 = t160 * rSges(6,1) + t159 * rSges(6,2) + rSges(6,3) * t306;
t161 = t223 * t301 - t224 * t225;
t162 = -t223 * t224 - t225 * t301;
t118 = t162 * rSges(6,1) + t161 * rSges(6,2) - rSges(6,3) * t302;
t69 = t117 * t302 + t118 * t306;
t258 = t221 * t305 - t224 * t296;
t280 = pkin(4) * t300;
t286 = -pkin(3) * t305 - pkin(7) * t306;
t129 = t258 - t280 + t286;
t292 = -t117 - t129;
t151 = (pkin(7) + t240) * t235 + t316 * t234;
t291 = -t151 - t158;
t290 = t323 * t231;
t273 = t226 * t297;
t274 = t226 * t298;
t289 = -pkin(3) * t273 - pkin(7) * t274;
t288 = pkin(2) * t303 + qJ(3) * t307;
t287 = qJ(3) * t303 + qJD(3) * t224;
t285 = pkin(3) * t301 + pkin(7) * t302;
t284 = t226 * pkin(2) + t224 * qJ(3);
t282 = qJD(4) * t236;
t253 = t263 * t224;
t121 = -t225 * t253 - t226 * t255;
t122 = -t223 * t253 + t226 * t254;
t71 = t122 * rSges(6,1) + t121 * rSges(6,2) + rSges(6,3) * t274;
t279 = t118 * t274 + t71 * t302 + t70 * t306;
t237 = sin(qJ(1));
t278 = t237 * t315;
t277 = pkin(4) * t282;
t172 = t224 * t294 - t300;
t173 = -t224 * t238 + t226 * t295;
t137 = -qJD(4) * t172 - t173 * t231;
t171 = -t224 * t295 - t226 * t238;
t247 = t171 * qJD(4);
t250 = t226 * t294 + t304;
t138 = t231 * t250 + t247;
t84 = t138 * rSges(5,1) + t137 * rSges(5,2) + rSges(5,3) * t274;
t131 = t172 * rSges(5,1) + t171 * rSges(5,2) + rSges(5,3) * t306;
t270 = t306 / 0.2e1;
t269 = -t302 / 0.2e1;
t268 = t298 / 0.2e1;
t186 = t226 * rSges(3,1) - rSges(3,2) * t224;
t267 = t291 * t234;
t177 = (-Icges(5,1) * t236 - t311) * t283;
t154 = t234 * t238 * t177;
t175 = (-Icges(5,5) * t236 - Icges(5,6) * t238) * t283;
t265 = -t235 * t175 + t154;
t262 = t229 * t277;
t261 = t235 * t277;
t260 = t224 * t268;
t259 = t226 * t268;
t169 = rSges(3,1) * t303 - rSges(3,2) * t307;
t135 = qJD(4) * t250 + t171 * t231;
t136 = qJD(4) * t173 + t172 * t231;
t257 = -rSges(5,1) * t136 - rSges(5,2) * t135;
t251 = t226 * t296 - t323;
t132 = -rSges(5,1) * t250 + rSges(5,2) * t173 - rSges(5,3) * t302;
t168 = t321 * t231;
t114 = Icges(6,4) * t162 + Icges(6,2) * t161 - Icges(6,6) * t302;
t116 = Icges(6,1) * t162 + Icges(6,4) * t161 - Icges(6,5) * t302;
t63 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t275;
t65 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t275;
t67 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t275;
t10 = -t235 * t63 + ((-t114 * t230 + t67) * t225 + (-t116 * t230 - t65) * t223) * t234;
t155 = -Icges(6,3) * t235 + (Icges(6,5) * t225 - Icges(6,6) * t223) * t234;
t20 = t119 * t156 + t120 * t157 + t146 * t161 + t147 * t162 + (-t145 * t226 + t155 * t307) * t234;
t21 = t121 * t156 + t122 * t157 + t146 * t159 + t147 * t160 + (t145 * t224 + t155 * t303) * t234;
t111 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t306;
t113 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t306;
t115 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t306;
t43 = -t111 * t235 + (-t113 * t223 + t115 * t225) * t234;
t112 = Icges(6,5) * t162 + Icges(6,6) * t161 - Icges(6,3) * t302;
t44 = -t112 * t235 + (-t114 * t223 + t116 * t225) * t234;
t72 = t155 * t306 + t156 * t159 + t157 * t160;
t73 = -t155 * t302 + t156 * t161 + t157 * t162;
t64 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t274;
t66 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t274;
t68 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t274;
t9 = -t235 * t64 + ((-t113 * t230 + t68) * t225 + (-t115 * t230 - t66) * t223) * t234;
t249 = (t21 + t9) * t270 + (t10 + t20) * t269 + (t44 + t73) * t260 + (t43 + t72) * t259;
t142 = rSges(4,1) * t301 - rSges(4,2) * t302 + t224 * rSges(4,3) + t284;
t27 = t111 * t306 + t113 * t159 + t115 * t160;
t28 = t112 * t306 + t114 * t159 + t116 * t160;
t29 = -t111 * t302 + t113 * t161 + t115 * t162;
t30 = -t112 * t302 + t114 * t161 + t116 * t162;
t248 = -(-t20 * t235 + ((t119 * t113 + t120 * t115 + t161 * t66 + t162 * t68) * t224 + t29 * t303 - (t119 * t114 + t120 * t116 + t161 * t65 + t162 * t67) * t226 + t30 * t307 + ((t111 * t307 - t226 * t64) * t224 - (t112 * t307 - t226 * t63) * t226) * t234) * t234) * t302 - t235 * (-t314 + ((t231 * t43 - t10) * t226 + (t231 * t44 + t9) * t224) * t234) + (-t21 * t235 + ((t121 * t113 + t122 * t115 + t159 * t66 + t160 * t68) * t224 + t27 * t303 - (t121 * t114 + t122 * t116 + t159 * t65 + t160 * t67) * t226 + t28 * t307 + ((t111 * t303 + t224 * t64) * t224 - (t112 * t303 + t224 * t63) * t226) * t234) * t234) * t306 + (-t235 * t73 + (t224 * t29 - t226 * t30) * t234) * t275 + (-t235 * t72 + (t224 * t27 - t226 * t28) * t234) * t274;
t246 = t322 * t224 + t226 * t261 - t231 * t280;
t219 = t224 * pkin(2);
t99 = -t226 * qJ(3) + t131 + t219 - t286;
t245 = t249 - t314;
t141 = -rSges(4,2) * t306 + rSges(4,1) * t305 + t219 + (-rSges(4,3) - qJ(3)) * t226;
t57 = -qJD(3) * t226 + t288 - t289 + t84;
t100 = -t132 + t284 + t285;
t166 = -Icges(5,6) * t235 + (-Icges(5,2) * t236 + t311) * t234;
t167 = -Icges(5,5) * t235 + (Icges(5,1) * t238 - t312) * t234;
t243 = -t293 + (-t166 * t238 - t167 * t236) * qJD(4);
t107 = rSges(4,2) * t275 + rSges(4,3) * t303 + (-rSges(4,1) * t235 - pkin(2)) * t307 + t287;
t108 = rSges(4,1) * t273 + rSges(4,3) * t307 + (-rSges(4,2) * t298 - qJD(3)) * t226 + t288;
t90 = -t118 - t251 + t284;
t89 = t219 + (-pkin(4) * t236 - qJ(3)) * t226 + t258 + t117;
t56 = (-pkin(3) * t235 - pkin(2) + (-rSges(5,3) - pkin(7)) * t234) * t307 + t257 + t287;
t125 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t306;
t127 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t306;
t78 = Icges(5,5) * t138 + Icges(5,6) * t137 + Icges(5,3) * t274;
t80 = Icges(5,4) * t138 + Icges(5,2) * t137 + Icges(5,6) * t274;
t82 = Icges(5,1) * t138 + Icges(5,4) * t137 + Icges(5,5) * t274;
t13 = -t235 * t78 + (-t236 * t80 + t238 * t82 + (-t125 * t238 - t127 * t236) * qJD(4)) * t234;
t126 = -Icges(5,4) * t250 + Icges(5,2) * t173 - Icges(5,6) * t302;
t128 = -Icges(5,1) * t250 + Icges(5,4) * t173 - Icges(5,5) * t302;
t77 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t275;
t79 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t275;
t81 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t275;
t14 = -t235 * t77 + (-t236 * t79 + t238 * t81 + (-t126 * t238 - t128 * t236) * qJD(4)) * t234;
t165 = -Icges(5,3) * t235 + (Icges(5,5) * t238 - Icges(5,6) * t236) * t234;
t24 = t135 * t166 + t136 * t167 + t173 * t176 - t250 * t177 + (t165 * t307 - t175 * t226) * t234;
t25 = t137 * t166 + t138 * t167 + t171 * t176 + t172 * t177 + (t165 * t303 + t175 * t224) * t234;
t123 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t306;
t49 = -t123 * t235 + (-t125 * t236 + t127 * t238) * t234;
t124 = -Icges(5,5) * t250 + Icges(5,6) * t173 - Icges(5,3) * t302;
t50 = -t124 * t235 + (-t126 * t236 + t128 * t238) * t234;
t74 = t234 * t243 + t265;
t85 = t165 * t306 + t166 * t171 + t167 * t172;
t86 = -t165 * t302 + t166 * t173 - t167 * t250;
t242 = (-t74 - t53) * t235 + t249 + (t25 + t13) * t270 + (t24 + t14) * t269 + (t50 + t86) * t260 + (t49 + t85) * t259;
t34 = -t224 * t261 + (-qJD(3) + t322) * t226 + t71 + t288 + t290;
t33 = (-rSges(6,3) * t234 - t221 * t235 - pkin(2)) * t307 - t246 + t256 + t287;
t241 = t144 + t154 + (-t175 - t145) * t235 + (t243 + t244) * t234;
t239 = cos(qJ(1));
t228 = t239 * pkin(1);
t227 = t237 * pkin(1);
t222 = t239 * t315;
t180 = t186 + t228;
t179 = t227 - t321;
t178 = (-rSges(5,1) * t236 - rSges(5,2) * t238) * t283;
t170 = -rSges(5,3) * t235 + (rSges(5,1) * t238 - rSges(5,2) * t236) * t234;
t150 = t169 + t222;
t149 = t168 - t278;
t140 = t142 + t228;
t139 = t227 + t141;
t130 = t251 + t285;
t106 = t235 * t118;
t105 = t222 + t108;
t104 = t107 - t278;
t98 = t132 * t235 - t170 * t302;
t97 = -t131 * t235 - t170 * t306;
t96 = t100 + t228;
t95 = t227 + t99;
t94 = pkin(4) * t247 - t226 * t271 + t289 + t290;
t93 = (-pkin(7) * t234 + t316 * t235) * t307 + t246;
t92 = -t158 * t302 + t106;
t91 = -t117 * t235 - t158 * t306;
t88 = t228 + t90;
t87 = t227 + t89;
t83 = rSges(5,3) * t275 - t257;
t55 = t222 + t57;
t54 = t56 - t278;
t52 = -t235 * t84 + (-t170 * t303 - t178 * t224) * t234;
t51 = t235 * t83 + (t170 * t307 - t178 * t226) * t234;
t48 = t130 * t235 + t226 * t267 + t106;
t47 = t224 * t267 + t235 * t292;
t42 = -t124 * t302 + t126 * t173 - t128 * t250;
t41 = -t123 * t302 + t125 * t173 - t127 * t250;
t40 = t124 * t306 + t126 * t171 + t128 * t172;
t39 = t123 * t306 + t125 * t171 + t127 * t172;
t38 = -t235 * t71 + (-t158 * t303 - t308) * t234;
t37 = -t148 * t302 + t313;
t32 = t222 + t34;
t31 = t33 - t278;
t26 = (t129 * t226 + t130 * t224) * t234 + t69;
t17 = t224 * t262 + (-t71 - t94) * t235 + (t291 * t303 - t308) * t234;
t16 = t151 * t275 + t235 * t93 + (-t148 * t234 + t262) * t226 + t313;
t15 = -t117 * t275 + t279;
t4 = ((t130 * t231 + t94) * t226 + (t231 * t292 + t93) * t224) * t234 + t279;
t1 = [(t149 * t180 + t150 * t179) * t320 + (t104 * t140 + t105 * t139) * t319 + (t54 * t96 + t55 * t95) * t318 - t223 * t157 * t299 + (t31 * t88 + t32 * t87) * t317 + t265 + t266 + (-t223 * t146 - t166 * t281 - t167 * t282 - t272 - t293) * t234; m(3) * (t149 * t186 - t150 * t321 + t168 * t180 + t169 * t179) + m(4) * (t104 * t142 + t105 * t141 + t107 * t140 + t108 * t139) + m(5) * (t100 * t54 + t55 * t99 + t56 * t96 + t57 * t95) + m(6) * (t31 * t90 + t32 * t89 + t33 * t88 + t34 * t87) + t241; (t33 * t90 + t34 * t89) * t317 + (t100 * t56 + t57 * t99) * t318 + (t107 * t142 + t108 * t141) * t319 + (t168 * t186 - t169 * t321) * t320 + t241; m(4) * ((-t139 * t231 - t104) * t226 + (t140 * t231 - t105) * t224) + m(5) * ((-t231 * t95 - t54) * t226 + (t231 * t96 - t55) * t224) + m(6) * ((-t231 * t87 - t31) * t226 + (t231 * t88 - t32) * t224); m(6) * ((-t231 * t89 - t33) * t226 + (t231 * t90 - t34) * t224) + m(5) * ((-t231 * t99 - t56) * t226 + (t100 * t231 - t57) * t224) + m(4) * ((-t141 * t231 - t107) * t226 + (t142 * t231 - t108) * t224); 0; t242 + m(5) * (t51 * t96 + t52 * t95 + t54 * t98 + t55 * t97) + m(6) * (t16 * t88 + t17 * t87 + t31 * t48 + t32 * t47); t242 + m(6) * (t16 * t90 + t17 * t89 + t33 * t48 + t34 * t47) + m(5) * (t100 * t51 + t52 * t99 + t56 * t98 + t57 * t97); m(5) * ((-t231 * t97 - t51) * t226 + (t231 * t98 - t52) * t224) + m(6) * ((-t231 * t47 - t16) * t226 + (t231 * t48 - t17) * t224); (t98 * t51 + t97 * t52 + (t131 * t226 + t132 * t224) * ((t132 * t231 + t84) * t226 + (-t131 * t231 + t83) * t224) * t229) * t318 - t235 * (-t74 * t235 + ((t231 * t49 - t14) * t226 + (t231 * t50 + t13) * t224) * t234) + (-t235 * t85 + (t224 * t39 - t226 * t40) * t234) * t274 + (-t25 * t235 + ((t137 * t125 + t138 * t127 + t171 * t80 + t172 * t82) * t224 + t39 * t303 - (t137 * t126 + t138 * t128 + t171 * t79 + t172 * t81) * t226 + t40 * t307 + ((t123 * t303 + t224 * t78) * t224 - (t124 * t303 + t224 * t77) * t226) * t234) * t234) * t306 + (-t235 * t86 + (t224 * t41 - t226 * t42) * t234) * t275 - (-t24 * t235 + ((t135 * t125 + t136 * t127 + t173 * t80 - t250 * t82) * t224 + t41 * t303 - (t135 * t126 + t136 * t128 + t173 * t79 - t250 * t81) * t226 + t42 * t307 + ((t123 * t307 - t226 * t78) * t224 - (t124 * t307 - t226 * t77) * t226) * t234) * t234) * t302 + (t16 * t48 + t17 * t47 + t26 * t4) * t317 + t248; m(6) * (t31 * t92 + t32 * t91 + t37 * t88 + t38 * t87) + t245; m(6) * (t33 * t92 + t34 * t91 + t37 * t90 + t38 * t89) + t245; m(6) * ((-t231 * t91 - t37) * t226 + (t231 * t92 - t38) * t224); m(6) * (t15 * t26 + t16 * t92 + t17 * t91 + t37 * t48 + t38 * t47 + t4 * t69) + t248; (t15 * t69 + t37 * t92 + t38 * t91) * t317 + t248;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
