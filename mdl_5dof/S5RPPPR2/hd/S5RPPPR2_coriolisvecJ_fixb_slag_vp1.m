% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:20
% EndTime: 2020-01-03 11:22:39
% DurationCPUTime: 7.98s
% Computational Cost: add. (8798->432), mult. (23999->570), div. (0->0), fcn. (27434->10), ass. (0->213)
t214 = sin(pkin(9));
t217 = cos(pkin(9));
t218 = cos(pkin(7));
t216 = sin(pkin(7));
t317 = cos(pkin(8));
t267 = t216 * t317;
t173 = t214 * t267 + t218 * t217;
t148 = qJD(5) * t173 + qJD(1);
t220 = sin(qJ(1));
t266 = t220 * t317;
t215 = sin(pkin(8));
t222 = cos(qJ(1));
t305 = t222 * t215;
t176 = t218 * t266 - t305;
t309 = t216 * t220;
t141 = t176 * t217 + t214 * t309;
t265 = t222 * t317;
t306 = t220 * t215;
t175 = t218 * t306 + t265;
t219 = sin(qJ(5));
t221 = cos(qJ(5));
t108 = -t141 * t219 + t175 * t221;
t109 = t141 * t221 + t175 * t219;
t140 = t176 * t214 - t217 * t309;
t55 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t140;
t316 = Icges(6,4) * t109;
t58 = Icges(6,2) * t108 + Icges(6,6) * t140 + t316;
t104 = Icges(6,4) * t108;
t61 = Icges(6,1) * t109 + Icges(6,5) * t140 + t104;
t14 = t108 * t58 + t109 * t61 + t140 * t55;
t178 = t218 * t265 + t306;
t308 = t216 * t222;
t144 = -t178 * t214 + t217 * t308;
t143 = t178 * t217 + t214 * t308;
t177 = t218 * t305 - t266;
t111 = t143 * t221 + t177 * t219;
t112 = t143 * t219 - t177 * t221;
t56 = Icges(6,5) * t111 - Icges(6,6) * t112 - Icges(6,3) * t144;
t315 = Icges(6,4) * t111;
t59 = -Icges(6,2) * t112 - Icges(6,6) * t144 + t315;
t105 = Icges(6,4) * t112;
t62 = Icges(6,1) * t111 - Icges(6,5) * t144 - t105;
t15 = -t108 * t59 - t109 * t62 - t140 * t56;
t252 = t14 * t140 + t144 * t15;
t174 = -t218 * t214 + t217 * t267;
t310 = t215 * t216;
t138 = -t174 * t219 + t221 * t310;
t139 = t174 * t221 + t219 * t310;
t88 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t173;
t314 = Icges(6,4) * t139;
t89 = Icges(6,2) * t138 + Icges(6,6) * t173 + t314;
t133 = Icges(6,4) * t138;
t90 = Icges(6,1) * t139 + Icges(6,5) * t173 + t133;
t28 = t108 * t89 + t109 * t90 + t140 * t88;
t5 = qJD(5) * t252 + t28 * t148;
t290 = qJD(5) * t144;
t269 = t290 / 0.2e1;
t313 = qJ(2) * t220;
t188 = pkin(1) * t222 + t313;
t294 = qJD(1) * t222;
t295 = qJD(1) * t220;
t298 = pkin(1) * t294 + qJ(2) * t295;
t362 = -qJD(1) * t188 + t298;
t16 = -t111 * t61 + t112 * t58 + t144 * t55;
t347 = -t111 * t62 + t112 * t59 + t144 * t56;
t251 = t140 * t16 - t347 * t144;
t21 = -t138 * t59 - t139 * t62 - t173 * t56;
t326 = pkin(2) * t218;
t255 = qJ(3) * t216 + t326;
t180 = t255 * t222;
t292 = qJD(3) * t216;
t197 = t220 * t292;
t293 = qJD(2) * t222;
t263 = t197 - t293;
t276 = t218 * t294;
t277 = t216 * t294;
t281 = pkin(2) * t276 + qJ(3) * t277 + t197;
t354 = -qJD(1) * t180 - t263 + t281 + t362;
t65 = rSges(6,1) * t111 - rSges(6,2) * t112 - rSges(6,3) * t144;
t91 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t173;
t353 = t148 * t65 + t91 * t290;
t163 = qJD(4) * t175;
t159 = -qJD(1) * t266 + t215 * t276;
t160 = t178 * qJD(1);
t239 = t160 * pkin(3) + qJ(4) * t159 + t163;
t331 = pkin(3) * t178 + t177 * qJ(4);
t351 = -qJD(1) * t331 - t163 + t239 + t354;
t129 = t160 * t214 - t217 * t277;
t130 = t160 * t217 + t214 * t277;
t350 = t130 * pkin(4) + pkin(6) * t129 - t293;
t346 = t111 * t90 - t112 * t89 - t144 * t88;
t345 = t140 * t91;
t336 = t160 * rSges(4,1) - t159 * rSges(4,2) + rSges(4,3) * t277 - t293;
t335 = t130 * rSges(5,1) - t129 * rSges(5,2) + t159 * rSges(5,3) - t293;
t198 = t222 * t292;
t307 = t218 * t220;
t299 = pkin(2) * t307 + qJ(3) * t309;
t211 = qJD(2) * t220;
t213 = t220 * pkin(1);
t302 = qJD(1) * (-qJ(2) * t222 + t213) - t211;
t262 = qJD(1) * t299 - t198 + t302;
t264 = t176 * pkin(3) + qJ(4) * t175;
t289 = t177 * qJD(4);
t334 = -qJD(1) * t264 - t262 + t289;
t202 = rSges(3,1) * t307;
t333 = qJD(1) * (rSges(3,2) * t309 + rSges(3,3) * t222 - t202);
t320 = rSges(3,2) * t216;
t332 = (-qJD(1) * t320 - qJD(2)) * t222 + rSges(3,1) * t276 + rSges(3,3) * t295;
t260 = pkin(4) * t143 - pkin(6) * t144;
t236 = rSges(5,1) * t143 + rSges(5,2) * t144 + rSges(5,3) * t177;
t157 = t175 * qJD(1);
t158 = t176 * qJD(1);
t103 = pkin(3) * t158 + t157 * qJ(4) - t289;
t244 = pkin(1) + t255;
t297 = qJ(2) * t294 + t211;
t280 = t198 + t297;
t330 = -t244 * t295 - t103 + t280;
t101 = t141 * pkin(4) + pkin(6) * t140;
t64 = t109 * rSges(6,1) + t108 * rSges(6,2) + t140 * rSges(6,3);
t328 = -qJD(1) * t101 - t148 * t64 + t334;
t92 = t141 * rSges(5,1) - t140 * rSges(5,2) + t175 * rSges(5,3);
t327 = -qJD(1) * t92 + t334;
t323 = qJD(5) / 0.2e1;
t322 = Icges(6,1) * t138 - t314 - t89;
t321 = -Icges(6,2) * t139 + t133 + t90;
t318 = rSges(3,3) + qJ(2);
t172 = pkin(1) * t295 - t297;
t304 = -t255 * t295 - t172 + t198;
t303 = t180 + t188;
t206 = qJD(1) * t211;
t275 = qJD(1) * t292;
t301 = t222 * t275 + t206;
t257 = rSges(3,1) * t218 - t320;
t156 = rSges(3,3) * t220 + t222 * t257;
t118 = -t293 + (t156 + t188) * qJD(1);
t296 = qJD(1) * t118;
t291 = qJD(5) * t140;
t73 = -qJD(5) * t109 - t130 * t219 + t159 * t221;
t74 = qJD(5) * t108 + t130 * t221 + t159 * t219;
t38 = t74 * rSges(6,1) + t73 * rSges(6,2) + t129 * rSges(6,3);
t288 = -t103 + t304;
t161 = qJD(1) * (-t293 + t298);
t286 = qJD(1) * t281 + t220 * t275 + t161;
t285 = qJD(4) * t159 + t301;
t283 = t176 * rSges(4,1) - t175 * rSges(4,2) + rSges(4,3) * t309;
t282 = t331 + t303;
t279 = t213 + t299;
t278 = t216 * t295;
t127 = t158 * t214 - t217 * t278;
t274 = t127 * t323;
t273 = t129 * t323;
t271 = t291 / 0.2e1;
t270 = -t290 / 0.2e1;
t128 = t158 * t217 + t214 * t278;
t261 = -pkin(4) * t128 - pkin(6) * t127;
t258 = t163 + t263;
t256 = -rSges(4,1) * t158 + rSges(4,2) * t157;
t131 = t138 * qJD(5);
t132 = t139 * qJD(5);
t96 = rSges(6,1) * t131 - rSges(6,2) * t132;
t254 = t127 * t91 + t144 * t96;
t253 = -t129 * t91 - t140 * t96;
t250 = -t140 * t65 - t144 * t64;
t249 = t140 * (Icges(6,5) * t108 - Icges(6,6) * t109) + t144 * (Icges(6,5) * t112 + Icges(6,6) * t111);
t248 = qJD(1) * t239 + qJD(4) * t157 + t286;
t245 = pkin(1) + t257;
t241 = rSges(4,1) * t178 - rSges(4,2) * t177 + rSges(4,3) * t308;
t240 = t331 + t313;
t238 = t264 + t279;
t71 = qJD(5) * t111 - t128 * t219 + t157 * t221;
t72 = qJD(5) * t112 + t128 * t221 + t157 * t219;
t37 = rSges(6,1) * t72 + rSges(6,2) * t71 + rSges(6,3) * t127;
t237 = -rSges(5,1) * t128 + rSges(5,2) * t127 - rSges(5,3) * t157;
t233 = (Icges(6,1) * t108 - t316 - t58) * t140 + (Icges(6,1) * t112 + t315 + t59) * t144;
t232 = (-Icges(6,2) * t109 + t104 + t61) * t140 + (Icges(6,2) * t111 + t105 - t62) * t144;
t31 = Icges(6,5) * t72 + Icges(6,6) * t71 + Icges(6,3) * t127;
t32 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t129;
t33 = Icges(6,4) * t72 + Icges(6,2) * t71 + Icges(6,6) * t127;
t34 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t129;
t35 = Icges(6,1) * t72 + Icges(6,4) * t71 + Icges(6,5) * t127;
t36 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t129;
t229 = (-t111 * t36 + t112 * t34 + t127 * t55 + t144 * t32 + t58 * t71 + t61 * t72) * t140 - t127 * t347 + t129 * t16 + t144 * (-t111 * t35 + t112 * t33 - t127 * t56 + t144 * t31 - t59 * t71 - t62 * t72);
t228 = t127 * t15 + t129 * t14 + t140 * (t108 * t34 + t109 * t36 + t129 * t55 + t140 * t32 + t58 * t73 + t61 * t74) + t144 * (t108 * t33 + t109 * t35 - t129 * t56 + t140 * t31 - t59 * t73 - t62 * t74);
t20 = t138 * t58 + t139 * t61 + t173 * t55;
t7 = t131 * t61 - t132 * t58 + t138 * t34 + t139 * t36 + t173 * t32;
t8 = -t131 * t62 + t132 * t59 + t138 * t33 + t139 * t35 + t173 * t31;
t227 = t127 * t21 + t129 * t20 + t140 * t7 + t144 * t8;
t226 = -t127 * t64 - t129 * t65 + t140 * t37 - t144 * t38;
t115 = t332 * qJD(1) + t161;
t114 = t206 + (-t172 + t333) * qJD(1);
t97 = Icges(6,5) * t138 - Icges(6,6) * t139;
t95 = Icges(6,1) * t131 - Icges(6,4) * t132;
t94 = Icges(6,4) * t131 - Icges(6,2) * t132;
t93 = Icges(6,5) * t131 - Icges(6,6) * t132;
t87 = (t241 + t303) * qJD(1) + t263;
t82 = rSges(6,1) * t112 + rSges(6,2) * t111;
t81 = rSges(6,1) * t108 - rSges(6,2) * t109;
t70 = t336 * qJD(1) + t286;
t69 = (-rSges(4,3) * t278 + t256 + t304) * qJD(1) + t301;
t51 = (t236 + t282) * qJD(1) + t258;
t40 = t335 * qJD(1) + t248;
t39 = (t237 + t288) * qJD(1) + t285;
t30 = -qJD(3) * t218 + qJD(4) * t310 + qJD(5) * t250;
t26 = (t260 + t282) * qJD(1) + t258 + t353;
t25 = -t291 * t91 - t328;
t19 = t131 * t90 - t132 * t89 + t138 * t94 + t139 * t95 + t173 * t93;
t18 = t19 * t148;
t13 = -t148 * t37 + t254 * qJD(5) + (t261 + t288) * qJD(1) + t285;
t12 = qJD(1) * t350 + t253 * qJD(5) + t148 * t38 + t248;
t11 = t108 * t94 + t109 * t95 + t129 * t88 + t140 * t93 + t73 * t89 + t74 * t90;
t10 = -t111 * t95 + t112 * t94 + t127 * t88 + t144 * t93 + t71 * t89 + t72 * t90;
t9 = t226 * qJD(5);
t1 = [t5 * t270 + t18 + (t21 - t346) * t274 + (t20 + t28) * t273 + (t7 + t11) * t271 + (t13 * (t240 + t260 + t65) + t12 * (t101 + t238 + t64) + (-t12 * qJ(2) + t13 * t244) * t222 + (-t345 * qJD(5) + t261 - t328 + t330 - t37) * t26 + (-qJD(1) * t260 + t350 + t351 - t353 + t38) * t25) * m(6) + (t39 * (t236 + t240) + t40 * (t238 + t92) + (-t40 * qJ(2) + t244 * t39) * t222 + (t237 + t330 - t327) * t51 - (-qJD(1) * t236 + t335 + t351) * t327) * m(5) + (t69 * (t241 + t313) + t70 * (t279 + t283) + (-t70 * qJ(2) + t244 * t69) * t222 + (t256 + t280 + (-t326 - pkin(1) + (-rSges(4,3) - qJ(3)) * t216) * t295) * t87 + (-qJD(1) * t241 + t336 + t354 + t87) * (qJD(1) * t283 + t262)) * m(4) + (t118 * t297 + t115 * (t202 + t213) + (t114 * t318 - t115 * t320 - t245 * t296) * t220 + (rSges(3,3) * t296 + t114 * t245 - t115 * t318) * t222 + (-qJD(1) * t156 + t118 + t293 + t332 + t362) * (t302 - t333)) * m(3) + (t5 + t8 + t10) * t269; m(3) * (-t114 * t222 - t115 * t220) + m(4) * (-t220 * t70 - t222 * t69) + m(5) * (-t40 * t220 - t39 * t222) + m(6) * (-t12 * t220 - t13 * t222); -t9 * t218 * m(6) + 0.2e1 * (m(4) * (t220 * t69 - t222 * t70) / 0.2e1 + m(5) * (t220 * t39 - t222 * t40) / 0.2e1 + m(6) * (-t12 * t222 + t13 * t220) / 0.2e1) * t216; m(5) * (-t157 * t327 + t159 * t51 + t175 * t39 - t177 * t40) + m(6) * (-t12 * t177 + t13 * t175 + t157 * t25 + t159 * t26 + t9 * t310) + 0.2e1 * (-m(5) * (-t175 * t327 + t177 * t51) / 0.2e1 - m(6) * (t175 * t25 + t177 * t26) / 0.2e1) * qJD(1); t173 * (qJD(5) * t227 + t18) / 0.2e1 + t148 * (t173 * t19 + t227) / 0.2e1 + t129 * t5 / 0.2e1 + t140 * (qJD(5) * t228 + t11 * t148) / 0.2e1 + (t173 * t28 + t252) * t273 + (t11 * t173 + t228) * t271 + t127 * (qJD(5) * t251 - t346 * t148) / 0.2e1 + t144 * (qJD(5) * t229 + t10 * t148) / 0.2e1 + (-t173 * t346 + t251) * t274 + (t10 * t173 + t229) * t269 - t148 * ((t321 * t138 + t322 * t139 + t173 * t97) * t148 + (t138 * t232 + t139 * t233 + t173 * t249) * qJD(5)) / 0.2e1 - ((t108 * t321 + t109 * t322 + t140 * t97) * t148 + (t108 * t232 + t109 * t233 + t140 * t249) * qJD(5)) * t291 / 0.2e1 + ((-t111 * t322 + t112 * t321 + t144 * t97) * t148 + (-t111 * t233 + t112 * t232 + t144 * t249) * qJD(5)) * t270 + (t9 * t250 + t30 * t226 + t13 * (t144 * t91 + t173 * t65) + t26 * (-t173 * t37 + t254) + t12 * (t173 * t64 - t345) + t25 * (t173 * t38 + t253) - (t25 * t81 - t26 * t82) * t148 - (t30 * (t140 * t82 - t144 * t81) + (-t140 * t25 + t144 * t26) * (rSges(6,1) * t138 - rSges(6,2) * t139)) * qJD(5)) * m(6);];
tauc = t1(:);
