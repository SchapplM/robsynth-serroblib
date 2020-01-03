% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR9_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:08
% EndTime: 2019-12-31 16:56:24
% DurationCPUTime: 12.58s
% Computational Cost: add. (6718->671), mult. (17791->945), div. (0->0), fcn. (16709->6), ass. (0->325)
t242 = sin(qJ(1));
t244 = cos(qJ(3));
t245 = cos(qJ(1));
t240 = sin(qJ(4));
t243 = cos(qJ(4));
t371 = t242 * t243;
t241 = sin(qJ(3));
t373 = t241 * t245;
t164 = t240 * t373 + t371;
t369 = t243 * t245;
t372 = t242 * t240;
t165 = t241 * t369 - t372;
t158 = Icges(5,4) * t165;
t368 = t244 * t245;
t87 = Icges(5,2) * t164 + Icges(5,6) * t368 - t158;
t157 = Icges(5,4) * t164;
t89 = Icges(5,1) * t165 - Icges(5,5) * t368 - t157;
t298 = -t164 * t87 - t165 * t89;
t162 = -t241 * t372 + t369;
t163 = t240 * t245 + t241 * t371;
t370 = t242 * t244;
t382 = Icges(5,4) * t163;
t85 = Icges(5,2) * t162 - Icges(5,6) * t370 + t382;
t156 = Icges(5,4) * t162;
t88 = Icges(5,1) * t163 - Icges(5,5) * t370 + t156;
t403 = t162 * t85 + t163 * t88;
t82 = Icges(5,5) * t163 + Icges(5,6) * t162 - Icges(5,3) * t370;
t84 = -Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t368;
t453 = t298 + t403 + (-t242 * t82 - t245 * t84) * t244;
t340 = qJD(4) * t244;
t342 = qJD(3) * t245;
t186 = -t242 * t340 + t342;
t341 = qJD(4) * t241;
t219 = qJD(1) + t341;
t27 = t164 * t85 - t165 * t88 + t82 * t368;
t282 = Icges(5,5) * t243 - Icges(5,6) * t240;
t140 = Icges(5,3) * t241 + t244 * t282;
t380 = Icges(5,4) * t243;
t284 = -Icges(5,2) * t240 + t380;
t144 = Icges(5,6) * t241 + t244 * t284;
t381 = Icges(5,4) * t240;
t286 = Icges(5,1) * t243 - t381;
t148 = Icges(5,5) * t241 + t244 * t286;
t52 = t140 * t368 + t144 * t164 - t148 * t165;
t452 = -t186 * t27 - t219 * t52;
t426 = -pkin(1) - pkin(5);
t450 = t426 * t242;
t299 = rSges(5,1) * t243 - rSges(5,2) * t240;
t153 = rSges(5,3) * t241 + t244 * t299;
t339 = qJD(4) * t245;
t344 = qJD(3) * t242;
t185 = t244 * t339 + t344;
t408 = t244 * pkin(3);
t208 = pkin(6) * t241 + t408;
t300 = rSges(5,1) * t165 - rSges(5,2) * t164;
t93 = rSges(5,3) * t368 - t300;
t449 = t153 * t185 + t208 * t344 - t219 * t93;
t296 = t240 * t87 + t243 * t89;
t31 = t241 * t84 - t244 * t296;
t26 = t162 * t87 - t163 * t89 - t370 * t84;
t237 = t245 * rSges(4,3);
t374 = t241 * t242;
t152 = rSges(4,1) * t374 + rSges(4,2) * t370 + t237;
t205 = t245 * pkin(1) + t242 * qJ(2);
t238 = t245 * pkin(5);
t432 = t238 + t205;
t445 = t152 + t432;
t343 = qJD(3) * t244;
t324 = t242 * t343;
t346 = qJD(1) * t245;
t444 = t241 * t346 + t324;
t216 = Icges(4,4) * t370;
t379 = Icges(4,5) * t245;
t149 = Icges(4,1) * t374 + t216 + t379;
t383 = Icges(4,4) * t244;
t287 = Icges(4,1) * t241 + t383;
t150 = -Icges(4,5) * t242 + t245 * t287;
t197 = -Icges(4,2) * t241 + t383;
t171 = t197 * t245;
t255 = t242 * (t150 + t171) - t245 * (-Icges(4,2) * t374 + t149 + t216);
t384 = Icges(4,4) * t241;
t285 = Icges(4,2) * t244 + t384;
t145 = Icges(4,6) * t245 + t242 * t285;
t146 = -Icges(4,6) * t242 + t245 * t285;
t199 = Icges(4,1) * t244 - t384;
t173 = t199 * t242;
t174 = t199 * t245;
t256 = t242 * (t146 - t174) - t245 * (t145 - t173);
t443 = -t256 * t241 + t255 * t244;
t355 = t197 + t287;
t356 = -t285 + t199;
t442 = (t241 * t355 - t244 * t356) * qJD(1);
t439 = 0.2e1 * qJD(3);
t51 = -t140 * t370 + t144 * t162 + t148 * t163;
t438 = t185 * t26 + t51 * t219;
t139 = Icges(5,3) * t244 - t241 * t282;
t281 = t144 * t240 - t148 * t243;
t297 = t240 * t85 - t243 * t88;
t247 = t185 * (-t140 * t245 + t296) + t186 * (t140 * t242 + t297) + t219 * (t139 + t281);
t437 = t247 * t244;
t278 = t146 * t244 + t150 * t241;
t434 = t278 * t245;
t231 = t245 * qJ(2);
t201 = pkin(1) * t242 - t231;
t410 = pkin(5) * t242;
t312 = -t201 - t410;
t309 = -rSges(3,2) * t245 + t242 * rSges(3,3);
t433 = t205 + t309;
t101 = -qJD(3) * t174 + (t242 * t287 + t379) * qJD(1);
t283 = Icges(4,5) * t241 + Icges(4,6) * t244;
t142 = -Icges(4,3) * t242 + t245 * t283;
t350 = qJD(1) * t142;
t75 = t146 * t241 - t150 * t244;
t98 = qJD(1) * t145 - qJD(3) * t171;
t431 = qJD(3) * t75 + t101 * t241 + t244 * t98 + t350;
t308 = qJD(1) * t241 + qJD(4);
t323 = t244 * t342;
t430 = t242 * t308 - t323;
t188 = t285 * qJD(3);
t189 = t287 * qJD(3);
t195 = Icges(4,5) * t244 - Icges(4,6) * t241;
t429 = qJD(1) * t195 + qJD(3) * (t197 * t241 - t199 * t244) + t188 * t244 + t189 * t241;
t102 = qJD(1) * t150 + qJD(3) * t173;
t279 = t145 * t241 - t149 * t244;
t141 = Icges(4,3) * t245 + t242 * t283;
t351 = qJD(1) * t141;
t99 = qJD(1) * t146 + t197 * t344;
t428 = qJD(3) * t279 - t102 * t241 - t244 * t99 + t351;
t169 = (-Icges(5,2) * t243 - t381) * t244;
t253 = t185 * (Icges(5,2) * t165 + t157 - t89) + t186 * (-Icges(5,2) * t163 + t156 + t88) + t219 * (t148 + t169);
t172 = (-Icges(5,1) * t240 - t380) * t244;
t252 = t185 * (Icges(5,1) * t164 + t158 - t87) + t186 * (Icges(5,1) * t162 - t382 - t85) - t219 * (t144 - t172);
t337 = qJD(3) * qJD(4);
t322 = t241 * t337;
t134 = -qJD(1) * t185 + t242 * t322;
t425 = t134 / 0.2e1;
t135 = qJD(1) * t186 - t245 * t322;
t424 = t135 / 0.2e1;
t423 = -t185 / 0.2e1;
t422 = t185 / 0.2e1;
t421 = -t186 / 0.2e1;
t420 = t186 / 0.2e1;
t419 = -t219 / 0.2e1;
t418 = t219 / 0.2e1;
t417 = t241 / 0.2e1;
t416 = t242 / 0.2e1;
t415 = -t245 / 0.2e1;
t413 = rSges(3,2) - pkin(1);
t412 = rSges(5,3) + pkin(6);
t411 = pkin(3) * t241;
t409 = pkin(6) * t244;
t326 = t241 * t344;
t329 = t244 * t346;
t265 = t326 - t329;
t78 = -t219 * t371 + (-t245 * t308 - t324) * t240;
t79 = t308 * t369 + (-t219 * t240 + t243 * t343) * t242;
t41 = Icges(5,5) * t79 + Icges(5,6) * t78 + Icges(5,3) * t265;
t43 = Icges(5,4) * t79 + Icges(5,2) * t78 + Icges(5,6) * t265;
t45 = Icges(5,1) * t79 + Icges(5,4) * t78 + Icges(5,5) * t265;
t7 = (qJD(3) * t297 + t41) * t241 + (qJD(3) * t82 - t240 * t43 + t243 * t45 + (-t240 * t88 - t243 * t85) * qJD(4)) * t244;
t407 = t7 * t186;
t325 = t241 * t342;
t347 = qJD(1) * t244;
t264 = -t242 * t347 - t325;
t274 = t245 * t219;
t76 = -t240 * t430 + t243 * t274;
t77 = t240 * t274 + t243 * t430;
t40 = Icges(5,5) * t77 + Icges(5,6) * t76 + Icges(5,3) * t264;
t42 = Icges(5,4) * t77 + Icges(5,2) * t76 + Icges(5,6) * t264;
t44 = Icges(5,1) * t77 + Icges(5,4) * t76 + Icges(5,5) * t264;
t8 = (qJD(3) * t296 + t40) * t241 + (qJD(3) * t84 - t240 * t42 + t243 * t44 + (t240 * t89 - t243 * t87) * qJD(4)) * t244;
t406 = t8 * t185;
t405 = -qJD(1) / 0.2e1;
t147 = Icges(5,5) * t244 - t241 * t286;
t100 = qJD(3) * t147 + qJD(4) * t172;
t166 = (-Icges(5,5) * t240 - Icges(5,6) * t243) * t244;
t94 = qJD(3) * t139 + qJD(4) * t166;
t143 = Icges(5,6) * t244 - t241 * t284;
t97 = qJD(3) * t143 + qJD(4) * t169;
t18 = (qJD(3) * t281 + t94) * t241 + (qJD(3) * t140 + t100 * t243 - t240 * t97 + (-t144 * t243 - t148 * t240) * qJD(4)) * t244;
t321 = t244 * t337;
t56 = t140 * t241 - t244 * t281;
t404 = t18 * t219 + t56 * t321;
t401 = rSges(3,3) * t245;
t399 = rSges(5,3) * t244;
t235 = t242 * rSges(4,3);
t301 = rSges(4,1) * t241 + rSges(4,2) * t244;
t154 = t301 * t245 - t235;
t204 = rSges(4,1) * t244 - rSges(4,2) * t241;
t182 = t204 * t344;
t228 = qJD(2) * t242;
t67 = t182 + t228 + (t154 + t312) * qJD(1);
t395 = t245 * t67;
t30 = t241 * t82 - t244 * t297;
t394 = t30 * t134;
t393 = t31 * t135;
t334 = rSges(4,1) * t444 + rSges(4,2) * t329;
t345 = qJD(3) * t241;
t107 = (-rSges(4,2) * t345 - rSges(4,3) * qJD(1)) * t242 + t334;
t190 = t301 * qJD(3);
t246 = qJD(1) ^ 2;
t338 = qJD(1) * qJD(2);
t348 = qJD(1) * t242;
t354 = qJ(2) * t346 + t228;
t358 = qJD(1) * (-pkin(1) * t348 + t354) + t242 * t338;
t275 = -t246 * t410 + t358;
t48 = t190 * t342 + (t107 + t182) * qJD(1) + t275;
t392 = t48 * t245;
t177 = t204 * t245;
t106 = -qJD(3) * t177 + (t242 * t301 + t237) * qJD(1);
t229 = qJD(2) * t245;
t159 = qJD(1) * t205 - t229;
t223 = t245 * t338;
t307 = -t238 * t246 + t223;
t328 = t204 * t342;
t49 = -t190 * t344 + (-t106 - t159 + t328) * qJD(1) + t307;
t391 = t49 * t242;
t220 = pkin(3) * t374;
t178 = -pkin(6) * t370 + t220;
t357 = t163 * rSges(5,1) + t162 * rSges(5,2);
t91 = -rSges(5,3) * t370 + t357;
t386 = -t178 - t91;
t221 = pkin(6) * t368;
t180 = pkin(3) * t373 - t221;
t385 = t180 - t93;
t375 = t195 * t245;
t167 = t242 * t195;
t277 = t197 * t244 + t199 * t241;
t81 = t245 * t277 - t167;
t367 = t81 * qJD(1);
t151 = -t241 * t299 + t399;
t176 = (-rSges(5,1) * t240 - rSges(5,2) * t243) * t244;
t105 = qJD(3) * t151 + qJD(4) * t176;
t207 = t409 - t411;
t192 = qJD(3) * t207;
t366 = t105 + t192;
t359 = t153 + t208;
t353 = rSges(3,2) * t348 + rSges(3,3) * t346;
t352 = -qJD(1) * t201 + t228;
t349 = qJD(1) * t283;
t336 = -rSges(4,3) + t426;
t335 = t79 * rSges(5,1) + t78 * rSges(5,2) + rSges(5,3) * t326;
t54 = t245 * t141 + t145 * t370 + t149 * t374;
t55 = -t245 * t142 - t146 * t370 - t150 * t374;
t333 = pkin(3) * t444 + pkin(6) * t326;
t331 = t412 * t244;
t327 = t208 * t342;
t319 = -t347 / 0.2e1;
t317 = -t344 / 0.2e1;
t316 = t344 / 0.2e1;
t315 = t343 / 0.2e1;
t314 = -t342 / 0.2e1;
t305 = qJD(4) * t315;
t304 = rSges(5,1) * t77 + rSges(5,2) * t76;
t303 = -t399 + t411;
t25 = -t370 * t82 + t403;
t295 = t242 * t26 + t245 * t25;
t294 = t242 * t25 - t245 * t26;
t28 = t368 * t84 - t298;
t293 = t242 * t28 + t245 * t27;
t292 = t242 * t27 - t245 * t28;
t291 = t242 * t31 + t245 * t30;
t290 = t242 * t30 - t245 * t31;
t68 = qJD(1) * t445 - t229 - t328;
t289 = t242 * t67 - t245 * t68;
t288 = t242 * t93 + t245 * t91;
t280 = t145 * t244 + t149 * t241;
t267 = (t242 * t55 + t245 * t54) * qJD(3);
t136 = t242 * t141;
t57 = -t280 * t245 + t136;
t58 = -t242 * t142 + t434;
t266 = (t242 * t58 + t245 * t57) * qJD(3);
t69 = (-t152 * t242 - t154 * t245) * qJD(3);
t263 = t140 * t219 + t185 * t84 + t186 * t82;
t260 = (Icges(5,5) * t162 - Icges(5,6) * t163) * t186 + (Icges(5,5) * t164 + Icges(5,6) * t165) * t185 + t166 * t219;
t258 = -qJD(1) * t278 - qJD(3) * t375 + t351;
t257 = qJD(1) * t280 + qJD(3) * t167 + t350;
t254 = t277 * qJD(1) - t283 * qJD(3);
t29 = -t185 * t91 + t186 * t93 + (-t178 * t242 - t180 * t245) * qJD(3);
t38 = t228 + (t180 + t312) * qJD(1) + t449;
t39 = -t327 - t153 * t186 + t219 * t91 - t229 + (t178 + t432) * qJD(1);
t248 = t29 * t288 + (-t242 * t39 - t245 * t38) * t153;
t202 = rSges(3,2) * t242 + t401;
t181 = t208 * t245;
t179 = t208 * t242;
t175 = t204 * t242;
t129 = t153 * t245;
t128 = t153 * t242;
t127 = t148 * t245;
t126 = t148 * t242;
t125 = t144 * t245;
t124 = t144 * t242;
t121 = qJD(1) * t433 - t229;
t120 = t228 + (-t201 + t202) * qJD(1);
t117 = -pkin(6) * t329 + t333;
t116 = t264 * pkin(6) + (t241 * t348 - t323) * pkin(3);
t115 = rSges(5,1) * t164 + rSges(5,2) * t165;
t114 = rSges(5,1) * t162 - rSges(5,2) * t163;
t104 = t223 + (-qJD(1) * t309 - t159) * qJD(1);
t103 = qJD(1) * t353 + t358;
t80 = t242 * t277 + t375;
t71 = t80 * qJD(1);
t47 = -rSges(5,3) * t329 + t335;
t46 = rSges(5,3) * t264 + t304;
t37 = -t242 * t429 + t254 * t245;
t36 = t254 * t242 + t245 * t429;
t35 = t278 * qJD(3) + t101 * t244 - t241 * t98;
t34 = -qJD(3) * t280 + t102 * t244 - t241 * t99;
t24 = t266 - t367;
t23 = t71 + t267;
t16 = t100 * t163 + t140 * t265 + t144 * t78 + t148 * t79 + t162 * t97 - t370 * t94;
t15 = -t100 * t165 + t140 * t264 + t144 * t76 + t148 * t77 + t164 * t97 + t368 * t94;
t14 = t105 * t185 + t135 * t153 - t219 * t46 + (t192 * t242 - t340 * t93) * qJD(3) + (-t116 - t159 + t327) * qJD(1) + t307;
t13 = qJD(1) * t117 - t105 * t186 - t134 * t153 + t219 * t47 + (-t192 * t245 + t208 * t348 + t340 * t91) * qJD(3) + t275;
t12 = t185 * t31 + t186 * t30 + t219 * t56;
t11 = t185 * t28 - t452;
t10 = t186 * t25 + t438;
t9 = t134 * t93 - t135 * t91 - t185 * t47 + t186 * t46 + (t116 * t245 - t117 * t242 + (-t178 * t245 + t180 * t242) * qJD(1)) * qJD(3);
t6 = -t84 * t329 + t162 * t42 + t163 * t44 + t78 * t87 - t79 * t89 + (-t244 * t40 + t345 * t84) * t242;
t5 = -t82 * t329 + t162 * t43 + t163 * t45 + t78 * t85 + t79 * t88 + (-t244 * t41 + t345 * t82) * t242;
t4 = -t84 * t325 + t164 * t42 - t165 * t44 + t76 * t87 - t77 * t89 + (t245 * t40 - t348 * t84) * t244;
t3 = -t82 * t325 + t164 * t43 - t165 * t45 + t76 * t85 + t77 * t88 + (t245 * t41 - t348 * t82) * t244;
t2 = t134 * t25 + t135 * t26 + t16 * t219 + t185 * t6 + t186 * t5 + t321 * t51;
t1 = t134 * t27 + t135 * t28 + t15 * t219 + t185 * t4 + t186 * t3 + t321 * t52;
t17 = [t407 / 0.2e1 + t406 / 0.2e1 + t404 + t393 / 0.2e1 + t394 / 0.2e1 + (t71 + ((-t57 + t136 + t55) * t242 + (t58 - t434 + (t142 - t280) * t242 + t54) * t245) * qJD(3)) * t317 + t16 * t420 + (-qJD(3) * t277 + t188 * t241 - t189 * t244) * qJD(1) + ((t28 + t453) * t186 + t438) * t423 + t52 * t424 + t51 * t425 + (t15 + t10) * t422 + ((-t25 + t453) * t185 + t11 + t452) * t421 + (t24 + t367 + (t242 ^ 2 * t142 + (-t136 + t55 + (t142 + t280) * t245) * t245) * qJD(3)) * t314 + ((-t242 * t331 + t220 + t357 + t432) * t13 + (t245 * t303 - t221 + t231 + t300 + t450) * t14 + (t229 - t304 + (t241 * t412 + t408) * t342 + (t426 * t245 + (-qJ(2) - t303 + t409) * t242) * qJD(1)) * t38 + (t333 + t335 + t354 + t38 - t352 + (-t331 * t245 - t180 + t410 + t450) * qJD(1) - t449) * t39) * m(5) + (t49 * (-t235 + t312) + t67 * t229 + t48 * t445 + t68 * (-rSges(4,2) * t326 + t334 + t354) + (qJD(3) * t204 * t67 + t301 * t49) * t245 + (t336 * t395 + (t67 * (-qJ(2) - t301) + t68 * t336) * t242) * qJD(1) - (t182 - t67 + (t154 - t410) * qJD(1) + t352) * t68) * m(4) + (t104 * (t242 * t413 + t231 + t401) + t120 * t229 + t103 * t433 + t121 * (t353 + t354) + (t120 * t413 * t245 + (t120 * (-rSges(3,3) - qJ(2)) - t121 * pkin(1)) * t242) * qJD(1) - (qJD(1) * t202 - t120 + t352) * t121) * m(3) + (t35 + t36 + t23) * t316 + (qJD(1) * t75 + t34 + t37) * t342 / 0.2e1 + (t245 * t81 + (-t279 + t80) * t242) * qJD(3) * t405; 0.2e1 * (t13 * t415 + t14 * t416) * m(5) + 0.2e1 * (t391 / 0.2e1 - t392 / 0.2e1) * m(4) + 0.2e1 * (t103 * t415 + t104 * t416) * m(3); qJD(1) * (t242 * t35 + t245 * t34 + (t242 * t279 + t75 * t245) * qJD(1)) / 0.2e1 + ((-t241 * t356 - t244 * t355) * qJD(1) + (t241 * t255 + t244 * t256) * qJD(3)) * t405 + ((-t344 * t375 - t349) * t242 + (t442 + (t242 * t167 + t443) * qJD(3)) * t245) * t317 + ((t167 * t342 - t349) * t245 + (-t442 + (-t245 * t375 - t443) * qJD(3)) * t242) * t314 - t12 * t340 / 0.2e1 - t242 * t10 * t341 / 0.2e1 + (-qJD(1) * t290 + t242 * t8 + t245 * t7) * t418 + (((-t124 * t240 + t126 * t243 + t82) * t186 + (t125 * t240 - t127 * t243 + t84) * t185 + (-t143 * t240 + t147 * t243 + t140) * t219 + t56 * qJD(4)) * t244 + (qJD(4) * t290 + t247) * t241) * t419 + (-qJD(1) * t294 + t242 * t6 + t245 * t5) * t420 + ((t124 * t162 + t126 * t163) * t186 + (-t125 * t162 - t127 * t163) * t185 + (t143 * t162 + t147 * t163) * t219 + (t244 * t51 - t26 * t373) * qJD(4) + ((qJD(4) * t25 + t263) * t241 - t437) * t242) * t421 + t11 * t339 * t417 + (-qJD(1) * t292 + t242 * t4 + t245 * t3) * t422 + ((t124 * t164 - t126 * t165) * t186 + (-t125 * t164 + t127 * t165) * t185 + (t143 * t164 - t147 * t165) * t219 + (t244 * t52 + t27 * t374) * qJD(4) + ((-qJD(4) * t28 - t263) * t241 + t437) * t245) * t423 + t293 * t424 + t295 * t425 + t291 * t305 + ((-t13 * t359 - t39 * t366 - t9 * t385 + t29 * (t116 + t46) + (t29 * t386 + t359 * t38) * qJD(1)) * t245 + (t14 * t359 + t38 * t366 + t9 * t386 + t29 * (-t117 - t47) + (t29 * t385 + t359 * t39) * qJD(1)) * t242 - t38 * (qJD(1) * t181 + t129 * t219 + t151 * t185 + t207 * t344) - t39 * (qJD(1) * t179 + t128 * t219 - t151 * t186 - t207 * t342) - t29 * (-t128 * t185 - t129 * t186 - t179 * t344 - t181 * t342) - ((-t38 * t93 + t39 * t91) * t244 + t248 * t241) * qJD(4)) * m(5) + (0.2e1 * t69 * (t245 * t106 - t242 * t107 + (-t152 * t245 + t154 * t242) * qJD(1)) - t289 * t190 + (t391 - t392 + (t242 * t68 + t395) * qJD(1)) * t204 - (t175 * t68 + t177 * t67) * qJD(1) - (t69 * (-t175 * t242 - t177 * t245) - t289 * t301) * qJD(3)) * m(4) + (qJD(1) * t36 + t1 + ((t257 * t242 + t245 * t428) * t245 + (t258 * t242 - t245 * t431) * t242 + (-t57 * t242 + t58 * t245) * qJD(1)) * t439) * t416 + (qJD(1) * t37 + t2 + ((-t242 * t428 + t257 * t245) * t245 + (t242 * t431 + t258 * t245) * t242 + (-t54 * t242 + t55 * t245) * qJD(1)) * t439) * t245 / 0.2e1 - (t267 + t23 + t10) * t348 / 0.2e1 + (t266 + t24 + t11) * t346 / 0.2e1; -t2 * t370 / 0.2e1 + (t241 * t51 - t244 * t294) * t425 + ((qJD(3) * t294 + t16) * t241 + (-qJD(1) * t295 + qJD(3) * t51 - t242 * t5 + t245 * t6) * t244) * t420 + t1 * t368 / 0.2e1 + (t241 * t52 - t244 * t292) * t424 + ((qJD(3) * t292 + t15) * t241 + (-qJD(1) * t293 + qJD(3) * t52 - t242 * t3 + t245 * t4) * t244) * t422 + t12 * t315 + (t393 + t394 + t404 + t406 + t407) * t417 + (t241 * t56 - t244 * t290) * t305 + ((qJD(3) * t290 + t18) * t241 + (-qJD(1) * t291 + qJD(3) * t56 - t242 * t7 + t245 * t8) * t244) * t418 + (t162 * t253 + t163 * t252 - t260 * t370) * t421 + (t253 * t164 - t165 * t252 + t260 * t368) * t423 + (t260 * t241 + (-t240 * t253 + t252 * t243) * t244) * t419 + (t241 * t314 + t242 * t319) * t11 + (t241 * t316 + t245 * t319) * t10 + ((qJD(3) * t248 + t13 * t91 - t14 * t93 - t38 * t46 + t39 * t47) * t241 + (t38 * (-qJD(3) * t93 + t105 * t245) + t39 * (qJD(3) * t91 + t105 * t242) - t9 * t288 + t29 * (-t242 * t46 - t245 * t47 - t346 * t93 + t348 * t91) + (t13 * t242 + t14 * t245 + (-t242 * t38 + t245 * t39) * qJD(1)) * t153) * t244 - t38 * (-t115 * t219 + t176 * t185) - t39 * (t114 * t219 - t176 * t186) - t29 * (-t114 * t185 + t115 * t186)) * m(5);];
tauc = t17(:);
