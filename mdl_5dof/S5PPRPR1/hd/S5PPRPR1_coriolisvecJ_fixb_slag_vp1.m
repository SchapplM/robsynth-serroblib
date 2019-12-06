% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:01:10
% DurationCPUTime: 12.32s
% Computational Cost: add. (11409->496), mult. (12438->783), div. (0->0), fcn. (11761->8), ass. (0->266)
t209 = sin(pkin(7));
t204 = t209 ^ 2;
t211 = cos(pkin(7));
t205 = t211 ^ 2;
t363 = t204 + t205;
t282 = qJD(3) * t363;
t207 = pkin(8) + qJ(3);
t201 = sin(t207);
t203 = cos(t207);
t206 = pkin(9) + qJ(5);
t200 = sin(t206);
t202 = cos(t206);
t277 = rSges(6,1) * t202 - rSges(6,2) * t200;
t132 = -rSges(6,3) * t203 + t201 * t277;
t302 = qJD(5) * t201;
t307 = qJD(3) * t209;
t183 = t211 * t302 + t307;
t379 = t132 * t183;
t306 = qJD(3) * t211;
t184 = t209 * t302 - t306;
t378 = t132 * t184;
t208 = sin(pkin(9));
t210 = cos(pkin(9));
t308 = qJD(3) * t203;
t338 = Icges(4,4) * t203;
t339 = Icges(4,4) * t201;
t377 = -(-Icges(4,1) * t201 - t338) * t308 + (-(-Icges(5,5) * t210 + Icges(5,6) * t208) * t201 - t339 + (-Icges(5,3) - Icges(4,2)) * t203) * qJD(3) * t201;
t186 = rSges(4,1) * t201 + rSges(4,2) * t203;
t376 = t186 * t282;
t362 = -rSges(5,1) * t210 + rSges(5,2) * t208;
t361 = rSges(5,3) * t203 + t201 * t362;
t323 = t209 * t210;
t324 = t208 * t211;
t180 = -t203 * t324 + t323;
t322 = t210 * t211;
t325 = t208 * t209;
t181 = t203 * t322 + t325;
t328 = t201 * t211;
t271 = t208 * (Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t328) - t210 * (Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t328);
t178 = -t203 * t325 - t322;
t179 = t203 * t323 - t324;
t329 = t201 * t209;
t272 = t208 * (Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t329) - t210 * (Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t329);
t263 = -Icges(4,2) * t201 + t338;
t266 = Icges(4,1) * t203 - t339;
t367 = -(Icges(4,6) * t209 + t211 * t263) * t203 - (Icges(4,5) * t209 + t211 * t266) * t201;
t368 = (-Icges(4,6) * t211 + t209 * t263) * t203 + (-Icges(4,5) * t211 + t209 * t266) * t201;
t72 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t329;
t73 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t328;
t375 = (t209 * t271 - t211 * t272) * t201 + (t209 * t73 - t211 * t72) * t203 + t209 * t367 + t211 * t368;
t374 = t377 * t211 + (-t201 * t271 - t203 * t73 - t367) * qJD(3);
t373 = t377 * t209 + (-t201 * t272 - t203 * t72 + t368) * qJD(3);
t304 = qJD(4) * t209;
t191 = t201 * t304;
t185 = t201 * pkin(3) - t203 * qJ(4);
t244 = qJD(3) * t185;
t123 = -t209 * t244 + t191;
t303 = qJD(4) * t211;
t193 = t201 * t303;
t124 = -t211 * t244 + t193;
t372 = -qJD(4) * t201 + t209 * t123 + t211 * t124 + t244 * t363;
t369 = t361 * t282;
t355 = t183 / 0.2e1;
t353 = t184 / 0.2e1;
t259 = Icges(6,5) * t202 - Icges(6,6) * t200;
t125 = -Icges(6,3) * t203 + t201 * t259;
t334 = Icges(6,4) * t202;
t261 = -Icges(6,2) * t200 + t334;
t127 = -Icges(6,6) * t203 + t201 * t261;
t335 = Icges(6,4) * t200;
t264 = Icges(6,1) * t202 - t335;
t129 = -Icges(6,5) * t203 + t201 * t264;
t321 = t211 * t200;
t167 = t202 * t209 - t203 * t321;
t326 = t203 * t211;
t168 = t200 * t209 + t202 * t326;
t336 = Icges(6,4) * t168;
t64 = Icges(6,2) * t167 + Icges(6,6) * t328 + t336;
t149 = Icges(6,4) * t167;
t66 = Icges(6,1) * t168 + Icges(6,5) * t328 + t149;
t273 = -t200 * t64 + t202 * t66;
t327 = t203 * t209;
t165 = -t200 * t327 - t202 * t211;
t166 = t202 * t327 - t321;
t337 = Icges(6,4) * t166;
t63 = Icges(6,2) * t165 + Icges(6,6) * t329 + t337;
t148 = Icges(6,4) * t165;
t65 = Icges(6,1) * t166 + Icges(6,5) * t329 + t148;
t274 = -t200 * t63 + t202 * t65;
t360 = -(-t125 * t211 - t273) * t183 - (-t125 * t209 - t274) * t184;
t212 = -pkin(6) - qJ(4);
t320 = qJ(4) + t212;
t198 = pkin(4) * t210 + pkin(3);
t350 = pkin(3) - t198;
t359 = t201 * t350 - t203 * t320;
t151 = (-Icges(6,2) * t202 - t335) * t201;
t301 = qJD(5) * t203;
t221 = t183 * (-Icges(6,2) * t168 + t149 + t66) + t184 * (-Icges(6,2) * t166 + t148 + t65) - t301 * (t129 + t151);
t292 = t203 * t306;
t294 = t201 * t306;
t92 = -qJD(5) * t168 + t200 * t294;
t93 = qJD(5) * t167 - t202 * t294;
t43 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t292;
t62 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t328;
t249 = t201 * t43 + t308 * t62;
t45 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t292;
t47 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t292;
t10 = t167 * t45 + t168 * t47 + t211 * t249 + t64 * t92 + t66 * t93;
t126 = Icges(6,3) * t201 + t203 * t259;
t150 = (-Icges(6,5) * t200 - Icges(6,6) * t202) * t201;
t57 = qJD(3) * t126 + qJD(5) * t150;
t248 = t125 * t308 + t201 * t57;
t26 = t167 * t64 + t168 * t66 + t328 * t62;
t61 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t329;
t25 = t167 * t63 + t168 * t65 + t328 * t61;
t343 = t209 * t25;
t269 = t211 * t26 + t343;
t40 = t125 * t328 + t127 * t167 + t129 * t168;
t340 = t40 * t201;
t128 = Icges(6,6) * t201 + t203 * t261;
t58 = qJD(3) * t128 + qJD(5) * t151;
t130 = Icges(6,5) * t201 + t203 * t264;
t152 = (-Icges(6,1) * t200 - t334) * t201;
t59 = qJD(3) * t130 + qJD(5) * t152;
t219 = -(t127 * t92 + t129 * t93 + t167 * t58 + t168 * t59 + t211 * t248) * t203 + (t203 * t269 + t340) * qJD(3);
t293 = t203 * t307;
t295 = t201 * t307;
t90 = -qJD(5) * t166 + t200 * t295;
t91 = qJD(5) * t165 - t202 * t295;
t42 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t293;
t250 = t201 * t42 + t308 * t61;
t44 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t293;
t46 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t293;
t9 = t167 * t44 + t168 * t46 + t211 * t250 + t63 * t92 + t65 * t93;
t358 = qJD(5) * t219 / 0.2e1 + t10 * t355 + t9 * t353;
t357 = -t201 * t282 / 0.2e1;
t356 = -t183 / 0.2e1;
t354 = -t184 / 0.2e1;
t352 = -t203 / 0.2e1;
t351 = qJD(3) / 0.2e1;
t24 = t165 * t64 + t166 * t66 + t329 * t62;
t342 = t211 * t24;
t39 = t125 * t329 + t127 * t165 + t129 * t166;
t341 = t39 * t201;
t330 = t125 * t203;
t122 = -t201 * t320 - t203 * t350;
t187 = t203 * pkin(3) + t201 * qJ(4);
t305 = qJD(4) * t203;
t153 = qJD(3) * t187 - t305;
t319 = -t122 * qJD(3) - t153;
t317 = -t185 + t359;
t316 = -t122 - t187;
t135 = rSges(5,3) * t201 - t203 * t362;
t314 = -t135 * qJD(3) - t153;
t313 = -t185 + t361;
t312 = -t135 - t187;
t176 = t187 * t209;
t177 = t187 * t211;
t311 = t209 * t176 + t211 * t177;
t199 = qJD(2) * t209;
t310 = t193 + t199;
t309 = qJD(2) * t211;
t300 = qJD(3) * qJD(4);
t133 = rSges(6,3) * t201 + t203 * t277;
t154 = (-rSges(6,1) * t200 - rSges(6,2) * t202) * t201;
t60 = qJD(3) * t133 + qJD(5) * t154;
t299 = -t60 + t319;
t298 = t123 * t307 + t124 * t306 + t201 * t300;
t297 = -t132 + t317;
t291 = t203 * t300;
t289 = t306 / 0.2e1;
t288 = -t301 / 0.2e1;
t287 = t301 / 0.2e1;
t286 = qJD(5) * t351;
t285 = qJD(3) * t317;
t284 = qJD(3) * t314;
t283 = qJD(3) * t313;
t281 = t191 - t309;
t280 = t201 * t286;
t279 = t203 * t286;
t189 = t209 * t291;
t49 = rSges(6,1) * t93 + rSges(6,2) * t92 + rSges(6,3) * t292;
t68 = rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t328;
t16 = -t49 * t301 - t183 * t60 + t189 + (t319 * t209 + (-t132 * t326 + t201 * t68) * qJD(5)) * qJD(3);
t190 = t211 * t291;
t48 = rSges(6,1) * t91 + rSges(6,2) * t90 + rSges(6,3) * t293;
t67 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t329;
t17 = t48 * t301 + t184 * t60 + t190 + (t319 * t211 + (t132 * t327 - t201 * t67) * qJD(5)) * qJD(3);
t276 = -t16 * t211 + t17 * t209;
t275 = t62 * t183 + t61 * t184;
t23 = t165 * t63 + t166 * t65 + t329 * t61;
t270 = t209 * t23 + t342;
t27 = t201 * t274 - t203 * t61;
t28 = t201 * t273 - t203 * t62;
t268 = t27 * t209 + t28 * t211;
t267 = -t209 * t68 + t211 * t67;
t260 = -Icges(4,5) * t201 - Icges(4,6) * t203;
t258 = -t127 * t200 + t129 * t202;
t255 = t209 * t279;
t254 = t211 * t279;
t251 = t176 * t307 + t177 * t306 + qJD(1) - t305;
t247 = t126 - t258;
t70 = -pkin(4) * t324 + t122 * t209;
t71 = pkin(4) * t325 + t122 * t211;
t22 = t183 * t67 - t184 * t68 + (t209 * t70 + t211 * t71) * qJD(3) + t251;
t246 = t22 * t267;
t241 = qJD(3) * t260;
t240 = t150 * t301 - t183 * (Icges(6,5) * t167 - Icges(6,6) * t168) - t184 * (Icges(6,5) * t165 - Icges(6,6) * t166);
t238 = qJD(3) * t359;
t237 = Icges(5,5) * t203 + (-Icges(5,1) * t210 + Icges(5,4) * t208) * t201;
t235 = Icges(5,6) * t203 + (-Icges(5,4) * t210 + Icges(5,2) * t208) * t201;
t231 = t201 * t240;
t80 = rSges(5,1) * t179 + rSges(5,2) * t178 + rSges(5,3) * t329;
t81 = rSges(5,1) * t181 + rSges(5,2) * t180 + rSges(5,3) * t328;
t34 = (t209 * t80 + t211 * t81) * qJD(3) + t251;
t230 = t34 * t361;
t228 = t22 * (-t201 * t198 - t203 * t212 + t185);
t227 = qJD(3) * t237;
t226 = qJD(3) * t235;
t222 = (Icges(6,1) * t167 - t336 - t64) * t183 + (Icges(6,1) * t165 - t337 - t63) * t184 - (-t127 + t152) * t301;
t220 = -(t127 * t90 + t129 * t91 + t165 * t58 + t166 * t59 + t209 * t248) * t203 + (t203 * t270 + t341) * qJD(3);
t41 = t201 * t258 - t330;
t218 = -((qJD(3) * t258 - t57) * t203 + (qJD(3) * t125 - t200 * t58 + t202 * t59 + (-t127 * t202 - t129 * t200) * qJD(5)) * t201) * t203 + (t41 * t201 + t203 * t268) * qJD(3);
t213 = (-t247 * t301 - t360) * t201;
t194 = t203 * t303;
t192 = t203 * t304;
t171 = t260 * t211;
t170 = t260 * t209;
t158 = t211 * t241;
t157 = t209 * t241;
t145 = -t186 * t307 - t309;
t144 = -t186 * t306 + t199;
t120 = t237 * t211;
t119 = t237 * t209;
t118 = t235 * t211;
t117 = t235 * t209;
t105 = t129 * t211;
t104 = t129 * t209;
t103 = t127 * t211;
t102 = t127 * t209;
t99 = t211 * t227;
t98 = t209 * t227;
t97 = t211 * t226;
t96 = t209 * t226;
t89 = rSges(6,1) * t167 - rSges(6,2) * t168;
t88 = rSges(6,1) * t165 - rSges(6,2) * t166;
t79 = t211 * t238;
t78 = t209 * t238;
t69 = t376 * qJD(3);
t55 = t209 * t283 + t281;
t54 = t211 * t283 + t310;
t51 = t211 * t284 + t190;
t50 = t209 * t284 + t189;
t33 = qJD(3) * t369 + t298;
t32 = t209 * t285 - t301 * t68 + t281 - t379;
t31 = t211 * t285 + t301 * t67 + t310 + t378;
t12 = t183 * t48 - t184 * t49 + (t209 * t78 + t211 * t79 + t267 * t301) * qJD(3) + t298;
t11 = t183 * t28 + t184 * t27 - t301 * t41;
t8 = t165 * t45 + t166 * t47 + t209 * t249 + t64 * t90 + t66 * t91;
t7 = t165 * t44 + t166 * t46 + t209 * t250 + t63 * t90 + t65 * t91;
t6 = t183 * t26 + t184 * t25 - t301 * t40;
t5 = t183 * t24 + t184 * t23 - t301 * t39;
t4 = (qJD(3) * t273 - t43) * t203 + (qJD(3) * t62 - t200 * t45 + t202 * t47 + (-t200 * t66 - t202 * t64) * qJD(5)) * t201;
t3 = (qJD(3) * t274 - t42) * t203 + (qJD(3) * t61 - t200 * t44 + t202 * t46 + (-t200 * t65 - t202 * t63) * qJD(5)) * t201;
t1 = qJD(5) * t220 + t183 * t8 + t184 * t7;
t2 = [-m(4) * t69 + m(5) * t33 + m(6) * t12; m(5) * (t209 * t51 - t211 * t50) + m(6) * t276; (((t103 * t200 - t105 * t202 + t62) * t183 + (t102 * t200 - t104 * t202 + t61) * t184 + t41 * qJD(5)) * t201 + ((t247 * t203 + (t128 * t200 - t130 * t202 - t125) * t201 + t268) * qJD(5) + t360) * t203) * t287 + (t209 * t8 - t211 * t7) * t353 + ((-t103 * t165 - t105 * t166) * t183 + (-t102 * t165 - t104 * t166) * t184 + (t341 + (-t128 * t165 - t130 * t166 + t342) * t203) * qJD(5) + (((t23 - t330) * qJD(5) + t275) * t203 + t213) * t209) * t354 + (t10 * t209 - t211 * t9) * t355 + ((-t103 * t167 - t105 * t168) * t183 + (-t102 * t167 - t104 * t168) * t184 + (t340 + (-t128 * t167 - t130 * t168 + t343) * t203) * qJD(5) + (((t26 - t330) * qJD(5) + t275) * t203 + t213) * t211) * t356 + t209 * t358 + (t209 * t28 - t211 * t27) * t280 - t211 * t1 / 0.2e1 - t11 * t302 / 0.2e1 + (t209 * t26 - t211 * t25) * t254 + (t209 * t24 - t211 * t23) * t255 + ((t158 * t209 + t180 * t97 + t181 * t99) * t209 + (-t180 * t96 - t181 * t98 + t373 * t211 + (-t157 - t374) * t209) * t211) * t307 + ((-t157 * t211 + t178 * t96 + t179 * t98) * t211 + (-t178 * t97 - t179 * t99 + t374 * t209 + (t158 - t373) * t211) * t209) * t306 - ((t118 * t180 + t120 * t181) * t307 + t171 * qJD(3) * t204 + (-t180 * t117 - t181 * t119 - t209 * t170 + t375) * t306) * t307 / 0.2e1 + (-(t117 * t178 + t119 * t179) * t306 + t170 * qJD(3) * t205 + (t178 * t118 + t179 * t120 - t211 * t171 + t375) * t307) * t289 + (-t31 * (t133 * t184 + t194) - t32 * (-t133 * t183 + t192) - (t211 * t228 + t31 * t316) * t306 - (t209 * t228 + t316 * t32) * t307 - ((-t31 * t67 + t32 * t68) * t201 + t246 * t203) * qJD(5) + t12 * t311 + (t17 * t297 + t31 * t299 + t12 * (t68 + t71)) * t211 + (t16 * t297 + t32 * t299 + t12 * (t67 + t70)) * t209 + (t372 + (-t378 + t49 + t79) * t211 + (t379 + t48 + t78) * t209) * t22) * m(6) + (-t54 * t194 - t55 * t192 - ((t211 * t230 + t312 * t54) * t211 + (t209 * t230 + t312 * t55) * t209) * qJD(3) + t33 * t311 + (t313 * t51 + t314 * t54 + t33 * t81) * t211 + (t313 * t50 + t314 * t55 + t33 * t80) * t209 + (t369 + t372) * t34) * m(5) + (-t69 * t363 + (-t144 * t211 - t145 * t209 + t376) * qJD(3) + t144 * t306 + t145 * t307) * (rSges(4,1) * t203 - rSges(4,2) * t201) * m(4) + ((-t3 + t6) * t211 + (t4 + t5) * t209) * t288; 0.2e1 * (t12 * t352 + t22 * t357) * m(6) + 0.2e1 * (t33 * t352 + t34 * t357) * m(5) + 0.2e1 * (m(5) * (qJD(3) * t34 + t209 * t50 + t211 * t51) / 0.2e1 + m(6) * (qJD(3) * t22 + t16 * t209 + t17 * t211) / 0.2e1) * t201; t203 * t6 * t289 + t328 * t358 + (t201 * t269 - t203 * t40) * t254 + ((t10 * t211 + t209 * t9) * t201 + t219) * t355 + t5 * t293 / 0.2e1 + t1 * t329 / 0.2e1 + (t201 * t270 - t203 * t39) * t255 + ((t209 * t7 + t211 * t8) * t201 + t220) * t353 + t201 * t11 * t351 + (qJD(5) * t218 + t183 * t4 + t184 * t3) * t352 + (t201 * t268 - t203 * t41) * t280 + ((t209 * t3 + t211 * t4) * t201 + t218) * t288 + (t167 * t221 + t168 * t222 - t211 * t231) * t356 + (t165 * t221 + t166 * t222 - t209 * t231) * t354 + (t240 * t203 + (-t200 * t221 + t222 * t202) * t201) * t287 + ((-t16 * t68 + t17 * t67 + t31 * t48 - t32 * t49 + (t246 + (t209 * t31 - t211 * t32) * t132) * qJD(3)) * t203 + (t31 * (-qJD(3) * t67 + t209 * t60) + t32 * (qJD(3) * t68 - t211 * t60) + t12 * t267 + t22 * (-t209 * t49 + t211 * t48) + t276 * t132) * t201 - t31 * (t154 * t184 + t301 * t88) - t32 * (-t154 * t183 - t301 * t89) - t22 * (t183 * t88 - t184 * t89)) * m(6);];
tauc = t2(:);
