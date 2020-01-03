% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:31
% DurationCPUTime: 13.68s
% Computational Cost: add. (12926->632), mult. (17293->907), div. (0->0), fcn. (16359->8), ass. (0->316)
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t229 = sin(qJ(4));
t232 = cos(qJ(4));
t228 = qJ(1) + pkin(7);
t224 = sin(t228);
t364 = t232 * t233;
t225 = cos(t228);
t369 = t225 * t229;
t146 = t224 * t364 - t369;
t141 = Icges(5,4) * t146;
t365 = t229 * t233;
t367 = t225 * t232;
t145 = t224 * t365 + t367;
t373 = t224 * t230;
t82 = -Icges(5,2) * t145 + Icges(5,6) * t373 + t141;
t140 = Icges(5,4) * t145;
t86 = -Icges(5,1) * t146 - Icges(5,5) * t373 + t140;
t448 = t229 * t82 + t232 * t86;
t79 = Icges(5,5) * t146 - Icges(5,6) * t145 + Icges(5,3) * t373;
t31 = -t448 * t230 - t233 * t79;
t295 = rSges(5,1) * t232 - rSges(5,2) * t229;
t163 = -rSges(5,3) * t233 + t230 * t295;
t339 = qJD(4) * t230;
t342 = qJD(3) * t225;
t171 = -t224 * t339 + t342;
t415 = t230 * pkin(3);
t207 = -pkin(6) * t233 + t415;
t338 = qJD(4) * t233;
t214 = qJD(1) - t338;
t296 = rSges(5,1) * t146 - rSges(5,2) * t145;
t88 = rSges(5,3) * t373 + t296;
t449 = -t163 * t171 - t207 * t342 - t214 * t88;
t343 = qJD(3) * t224;
t170 = t225 * t339 + t343;
t25 = -t145 * t82 - t146 * t86 + t373 * t79;
t372 = t224 * t232;
t147 = -t225 * t365 + t372;
t374 = t224 * t229;
t148 = t225 * t364 + t374;
t368 = t225 * t230;
t81 = Icges(5,5) * t148 + Icges(5,6) * t147 + Icges(5,3) * t368;
t388 = Icges(5,4) * t148;
t84 = Icges(5,2) * t147 + Icges(5,6) * t368 + t388;
t142 = Icges(5,4) * t147;
t87 = Icges(5,1) * t148 + Icges(5,5) * t368 + t142;
t26 = -t145 * t84 + t146 * t87 + t81 * t373;
t280 = Icges(5,5) * t232 - Icges(5,6) * t229;
t149 = -Icges(5,3) * t233 + t230 * t280;
t386 = Icges(5,4) * t232;
t281 = -Icges(5,2) * t229 + t386;
t151 = -Icges(5,6) * t233 + t230 * t281;
t387 = Icges(5,4) * t229;
t283 = Icges(5,1) * t232 - t387;
t153 = -Icges(5,5) * t233 + t230 * t283;
t50 = -t145 * t151 + t146 * t153 + t149 * t373;
t10 = t170 * t26 - t171 * t25 + t214 * t50;
t27 = t147 * t82 - t148 * t86 + t79 * t368;
t28 = t147 * t84 + t148 * t87 + t81 * t368;
t51 = t147 * t151 + t148 * t153 + t149 * t368;
t11 = t170 * t28 - t171 * t27 + t51 * t214;
t235 = qJD(1) ^ 2;
t443 = t224 * t225;
t218 = t224 * rSges(4,3);
t366 = t225 * t233;
t137 = rSges(4,1) * t366 - rSges(4,2) * t368 + t218;
t180 = t225 * pkin(2) + t224 * pkin(5);
t234 = cos(qJ(1));
t227 = t234 * pkin(1);
t438 = t227 + t180;
t442 = t137 + t438;
t440 = 0.2e1 * qJD(3);
t308 = t225 * rSges(3,1) - rSges(3,2) * t224;
t226 = Icges(4,4) * t233;
t282 = -Icges(4,2) * t230 + t226;
t194 = Icges(4,1) * t230 + t226;
t437 = t227 + t308;
t436 = t224 * t10 + t225 * t11;
t191 = Icges(4,5) * t233 - Icges(4,6) * t230;
t190 = Icges(4,5) * t230 + Icges(4,6) * t233;
t261 = qJD(3) * t190;
t389 = Icges(4,4) * t230;
t195 = Icges(4,1) * t233 - t389;
t135 = Icges(4,5) * t224 + t195 * t225;
t133 = Icges(4,6) * t224 + t225 * t282;
t377 = t133 * t230;
t278 = -t135 * t233 + t377;
t380 = Icges(4,3) * t225;
t433 = -t225 * t261 + (-t191 * t224 + t278 + t380) * qJD(1);
t131 = Icges(4,3) * t224 + t191 * t225;
t203 = Icges(4,4) * t373;
t371 = t224 * t233;
t385 = Icges(4,5) * t225;
t134 = Icges(4,1) * t371 - t203 - t385;
t382 = Icges(4,6) * t225;
t132 = Icges(4,4) * t371 - Icges(4,2) * t373 - t382;
t378 = t132 * t230;
t279 = -t134 * t233 + t378;
t432 = -t224 * t261 + (t131 + t279) * qJD(1);
t130 = Icges(4,5) * t371 - Icges(4,6) * t373 - t380;
t52 = -t225 * t130 - t224 * t279;
t192 = Icges(4,2) * t233 + t389;
t275 = t192 * t230 - t194 * t233;
t431 = t275 * qJD(1) + t191 * qJD(3);
t357 = -Icges(4,2) * t371 + t134 - t203;
t359 = t194 * t224 + t132;
t430 = -t230 * t357 - t233 * t359;
t150 = Icges(5,3) * t230 + t233 * t280;
t276 = -t151 * t229 + t153 * t232;
t285 = -t229 * t84 + t232 * t87;
t429 = t170 * (-t149 * t225 - t285) - t171 * (-t149 * t224 + t448) + t214 * (t150 - t276);
t175 = (-Icges(5,2) * t232 - t387) * t230;
t428 = t170 * (-Icges(5,2) * t148 + t142 + t87) - t171 * (-Icges(5,2) * t146 - t140 - t86) + t214 * (t153 + t175);
t337 = qJD(3) * qJD(4);
t322 = t233 * t337;
t128 = qJD(1) * t170 + t224 * t322;
t427 = t128 / 0.2e1;
t129 = qJD(1) * t171 + t225 * t322;
t426 = t129 / 0.2e1;
t425 = -t170 / 0.2e1;
t424 = t170 / 0.2e1;
t423 = -t171 / 0.2e1;
t422 = t171 / 0.2e1;
t421 = -t214 / 0.2e1;
t420 = t214 / 0.2e1;
t417 = -rSges(5,3) - pkin(6);
t231 = sin(qJ(1));
t416 = pkin(1) * t231;
t414 = t233 * pkin(3);
t345 = qJD(1) * t230;
t331 = t225 * t345;
t340 = qJD(3) * t233;
t259 = t224 * t340 + t331;
t341 = qJD(3) * t230;
t251 = t214 * t232 + t229 * t341;
t344 = qJD(1) * t233;
t305 = -qJD(4) + t344;
t76 = t224 * t251 - t305 * t369;
t250 = t214 * t229 - t232 * t341;
t77 = t224 * t250 + t305 * t367;
t40 = Icges(5,5) * t77 + Icges(5,6) * t76 + Icges(5,3) * t259;
t42 = Icges(5,4) * t77 + Icges(5,2) * t76 + Icges(5,6) * t259;
t44 = Icges(5,1) * t77 + Icges(5,4) * t76 + Icges(5,5) * t259;
t7 = (-qJD(3) * t448 - t40) * t233 + (qJD(3) * t79 - t229 * t42 + t232 * t44 + (t229 * t86 - t232 * t82) * qJD(4)) * t230;
t413 = t7 * t171;
t325 = t225 * t340;
t332 = t224 * t345;
t257 = t325 - t332;
t74 = t225 * t251 + t305 * t374;
t75 = t225 * t250 - t305 * t372;
t39 = Icges(5,5) * t75 + Icges(5,6) * t74 + Icges(5,3) * t257;
t41 = Icges(5,4) * t75 + Icges(5,2) * t74 + Icges(5,6) * t257;
t43 = Icges(5,1) * t75 + Icges(5,4) * t74 + Icges(5,5) * t257;
t8 = (qJD(3) * t285 - t39) * t233 + (qJD(3) * t81 - t229 * t41 + t232 * t43 + (-t229 * t87 - t232 * t84) * qJD(4)) * t230;
t412 = t8 * t170;
t411 = qJD(1) / 0.2e1;
t174 = (-Icges(5,5) * t229 - Icges(5,6) * t232) * t230;
t111 = qJD(3) * t150 + qJD(4) * t174;
t152 = Icges(5,6) * t230 + t233 * t281;
t112 = qJD(3) * t152 + qJD(4) * t175;
t154 = Icges(5,5) * t230 + t233 * t283;
t176 = (-Icges(5,1) * t229 - t386) * t230;
t113 = qJD(3) * t154 + qJD(4) * t176;
t24 = (qJD(3) * t276 - t111) * t233 + (qJD(3) * t149 - t112 * t229 + t113 * t232 + (-t151 * t232 - t153 * t229) * qJD(4)) * t230;
t323 = t230 * t337;
t57 = -t149 * t233 + t230 * t276;
t410 = t24 * t214 + t57 * t323;
t409 = rSges(4,1) * t233;
t407 = rSges(5,3) * t230;
t197 = rSges(4,1) * t230 + rSges(4,2) * t233;
t162 = t197 * t225;
t330 = t197 * t343;
t66 = qJD(1) * t442 - t330;
t405 = t162 * t66;
t350 = rSges(4,2) * t373 + t225 * rSges(4,3);
t136 = rSges(4,1) * t371 - t350;
t222 = t225 * pkin(5);
t179 = pkin(2) * t224 - t222;
t312 = -t179 - t416;
t328 = t197 * t342;
t65 = -t328 + (-t136 + t312) * qJD(1);
t400 = t224 * t65;
t208 = pkin(6) * t230 + t414;
t166 = t208 * t224;
t35 = (-t166 + t312) * qJD(1) + t449;
t398 = t225 * t35;
t397 = t225 * t65;
t396 = t31 * t128;
t32 = t230 * t285 - t233 * t81;
t395 = t32 * t129;
t392 = t166 + t88;
t168 = pkin(3) * t366 + pkin(6) * t368;
t90 = t148 * rSges(5,1) + t147 * rSges(5,2) + rSges(5,3) * t368;
t391 = t168 + t90;
t376 = t190 * t224;
t375 = t190 * t225;
t97 = -t224 * t275 - t375;
t363 = t97 * qJD(1);
t164 = t233 * t295 + t407;
t177 = (-rSges(5,1) * t229 - rSges(5,2) * t232) * t230;
t114 = qJD(3) * t164 + qJD(4) * t177;
t186 = qJD(3) * t208;
t362 = -t114 - t186;
t361 = -t224 * t130 - t134 * t366;
t360 = t224 * t131 + t135 * t366;
t358 = -t194 * t225 - t133;
t356 = -t192 * t225 + t135;
t354 = t163 + t207;
t346 = qJD(1) * t225;
t353 = rSges(4,2) * t332 + rSges(4,3) * t346;
t352 = -t192 + t195;
t351 = t194 + t282;
t348 = qJD(1) * t191;
t347 = qJD(1) * t224;
t336 = t235 * t416;
t335 = t235 * t227;
t334 = t75 * rSges(5,1) + t74 * rSges(5,2) + rSges(5,3) * t325;
t329 = t207 * t343;
t326 = t225 * t341;
t324 = -pkin(2) - t409;
t320 = t346 / 0.2e1;
t319 = -t343 / 0.2e1;
t316 = t342 / 0.2e1;
t315 = t341 / 0.2e1;
t310 = t222 - t416;
t115 = t135 * t371;
t307 = t225 * t131 - t115;
t306 = -t130 + t377;
t213 = pkin(5) * t346;
t304 = qJD(1) * (-pkin(2) * t347 + t213) - t336;
t301 = qJD(4) * t315;
t300 = -rSges(5,1) * t77 - rSges(5,2) * t76;
t178 = rSges(3,1) * t224 + rSges(3,2) * t225;
t297 = -rSges(4,2) * t230 + t409;
t294 = t224 * t26 - t225 * t25;
t293 = t224 * t25 + t225 * t26;
t292 = t224 * t28 - t225 * t27;
t291 = t224 * t27 + t225 * t28;
t290 = t224 * t32 - t225 * t31;
t289 = t224 * t31 + t225 * t32;
t288 = -t224 * t66 - t397;
t287 = -t224 * t90 + t225 * t88;
t70 = t132 * t233 + t134 * t230;
t71 = t133 * t233 + t135 * t230;
t277 = t136 * t224 + t137 * t225;
t272 = t230 * t39 + t340 * t81;
t271 = t230 * t40 + t340 * t79;
t161 = t197 * t224;
t53 = -t133 * t373 - t307;
t266 = (t224 * t53 - t225 * t52) * qJD(3);
t54 = -t132 * t368 - t361;
t55 = -t133 * t368 + t360;
t265 = (t224 * t55 - t225 * t54) * qJD(3);
t264 = t230 * t417 - pkin(2) - t414;
t263 = qJD(3) * t194;
t262 = qJD(3) * t192;
t258 = -t224 * t344 - t326;
t256 = t149 * t214 + t170 * t81 - t171 * t79;
t255 = (-Icges(5,5) * t145 - Icges(5,6) * t146) * t171 - (Icges(5,5) * t147 - Icges(5,6) * t148) * t170 - t174 * t214;
t254 = -t230 * t356 + t233 * t358;
t253 = t230 * t255;
t246 = (-t230 * t351 + t233 * t352) * qJD(1);
t244 = (Icges(5,1) * t147 - t388 - t84) * t170 - (-Icges(5,1) * t145 - t141 - t82) * t171 + (-t151 + t176) * t214;
t100 = -qJD(3) * t161 + (t225 * t297 + t218) * qJD(1);
t99 = rSges(4,1) * t258 - rSges(4,2) * t325 + t353;
t243 = t100 * t224 + t225 * t99 + (t136 * t225 - t137 * t224) * qJD(1);
t29 = t170 * t88 + t171 * t90 + qJD(2) + (t166 * t224 + t168 * t225) * qJD(3);
t36 = -t329 - t163 * t170 + t214 * t90 + (t168 + t438) * qJD(1);
t240 = t29 * t287 + (t224 * t35 - t225 * t36) * t163;
t182 = t282 * qJD(3);
t183 = t195 * qJD(3);
t237 = qJD(1) * t190 - t182 * t230 + t183 * t233 + (-t192 * t233 - t194 * t230) * qJD(3);
t236 = t429 * t230;
t189 = pkin(6) * t325;
t184 = t297 * qJD(3);
t173 = qJD(1) * t179;
t172 = t180 * qJD(1);
t167 = t207 * t225;
t165 = t207 * t224;
t125 = t163 * t225;
t124 = t163 * t224;
t123 = t153 * t225;
t122 = t153 * t224;
t121 = t151 * t225;
t120 = t151 * t224;
t110 = t259 * pkin(6) + (-t224 * t341 + t225 * t344) * pkin(3);
t109 = pkin(3) * t258 - pkin(6) * t332 + t189;
t108 = rSges(5,1) * t147 - rSges(5,2) * t148;
t107 = -rSges(5,1) * t145 - rSges(5,2) * t146;
t98 = -t225 * t275 + t376;
t78 = t98 * qJD(1);
t64 = qJD(3) * t277 + qJD(2);
t49 = -t335 - t184 * t342 + (-t100 - t172 + t330) * qJD(1);
t48 = -t184 * t343 + (t99 - t328) * qJD(1) + t304;
t46 = rSges(5,3) * t259 - t300;
t45 = -rSges(5,3) * t332 + t334;
t38 = t237 * t224 - t225 * t431;
t37 = t224 * t431 + t237 * t225;
t34 = -qJD(3) * t278 + t230 * (-t225 * t263 + (-t195 * t224 + t385) * qJD(1)) + t233 * (-t225 * t262 + (-t224 * t282 + t382) * qJD(1));
t33 = -t279 * qJD(3) + t230 * (qJD(1) * t135 - t224 * t263) + t233 * (qJD(1) * t133 - t224 * t262);
t30 = t243 * qJD(3);
t23 = t78 + t265;
t22 = t266 + t363;
t16 = t111 * t373 - t112 * t145 + t113 * t146 + t149 * t259 + t151 * t76 + t153 * t77;
t15 = t111 * t368 + t112 * t147 + t113 * t148 + t149 * t257 + t151 * t74 + t153 * t75;
t14 = -t335 - t114 * t171 + t128 * t163 - t214 * t46 + (-t186 * t225 - t339 * t88) * qJD(3) + (-t110 - t172 + t329) * qJD(1);
t13 = qJD(1) * t109 - t114 * t170 - t129 * t163 + t214 * t45 + (-t186 * t224 - t207 * t346 + t339 * t90) * qJD(3) + t304;
t12 = t170 * t32 - t171 * t31 + t214 * t57;
t9 = -t128 * t90 + t129 * t88 + t170 * t46 + t171 * t45 + (t109 * t225 + t110 * t224 + (t166 * t225 - t168 * t224) * qJD(1)) * qJD(3);
t6 = -t145 * t41 + t146 * t43 + t224 * t272 + t331 * t81 + t76 * t84 + t77 * t87;
t5 = -t145 * t42 + t146 * t44 + t224 * t271 + t331 * t79 + t76 * t82 - t77 * t86;
t4 = t147 * t41 + t148 * t43 + t225 * t272 - t332 * t81 + t74 * t84 + t75 * t87;
t3 = t147 * t42 + t148 * t44 + t225 * t271 - t332 * t79 + t74 * t82 - t75 * t86;
t2 = t128 * t25 + t129 * t26 + t16 * t214 + t170 * t6 - t171 * t5 + t323 * t50;
t1 = t128 * t27 + t129 * t28 + t15 * t214 + t170 * t4 - t171 * t3 + t323 * t51;
t17 = [m(3) * ((-t178 * t235 - t336) * t437 + (-t335 + (-0.2e1 * t308 - t227 + t437) * t235) * (-t178 - t416)) + (-qJD(3) * t275 + t182 * t233 + t183 * t230) * qJD(1) + t16 * t423 + t15 * t424 + t51 * t426 + t50 * t427 - t413 / 0.2e1 + t412 / 0.2e1 + t410 + (t78 + ((t53 - t115 + (t131 + t378) * t225 + t361) * t225 + t360 * t224) * qJD(3)) * t316 + t395 / 0.2e1 + t396 / 0.2e1 + (t22 - t363 + ((t225 * t306 - t360 + t55) * t225 + (t224 * t306 + t307 + t54) * t224) * qJD(3)) * t319 + (t34 + t37) * t343 / 0.2e1 + (t422 + t423) * t11 + (t14 * (-t296 + t310) + t35 * t300 + t13 * (t438 + t391) + t36 * (-pkin(3) * t326 + t189 + t213 + t334) + (t14 * t264 + t35 * (t233 * t417 + t415) * qJD(3)) * t224 + ((-t231 * t36 - t234 * t35) * pkin(1) + t264 * t398 + (-t35 * pkin(5) + t36 * (-pkin(2) - t208 - t407)) * t224) * qJD(1) - (-t173 - t35 + (-t166 - t416) * qJD(1) + t449) * t36) * m(5) + (-(-t328 - t173 - t65 + (-t136 - t416) * qJD(1)) * t66 + t49 * (t224 * t324 + t310 + t350) + t48 * t442 + t66 * (t213 + t353) + (t197 * t400 - t405) * qJD(3) + ((-t231 * t66 - t234 * t65) * pkin(1) + (-pkin(2) - t297) * t397 + (t65 * (-rSges(4,3) - pkin(5)) + t66 * t324) * t224) * qJD(1)) * m(4) - (t33 + t38 + t23) * t342 / 0.2e1 + ((t70 + t97) * t224 + (t71 + t98) * t225) * qJD(3) * t411; m(4) * t30 + m(5) * t9; t290 * t301 + (((t121 * t229 - t123 * t232 + t81) * t170 - (t120 * t229 - t122 * t232 + t79) * t171 + (-t152 * t229 + t154 * t232 + t149) * t214 + t57 * qJD(4)) * t230 + (qJD(4) * t289 - t429) * t233) * t421 + ((t121 * t145 - t123 * t146) * t170 - (t120 * t145 - t122 * t146) * t171 + (-t145 * t152 + t146 * t154) * t214 + (t230 * t50 + t26 * t366) * qJD(4) + ((qJD(4) * t25 + t256) * t233 + t236) * t224) * t422 + (qJD(1) * t293 + t224 * t6 - t225 * t5) * t423 + (qJD(1) * t291 + t224 * t4 - t225 * t3) * t424 + ((-t121 * t147 - t123 * t148) * t170 - (-t120 * t147 - t122 * t148) * t171 + (t147 * t152 + t148 * t154) * t214 + (t230 * t51 + t27 * t371) * qJD(4) + ((qJD(4) * t28 + t256) * t233 + t236) * t225) * t425 + t292 * t426 + t294 * t427 + (t224 * t34 - t225 * t33 + (t224 * t70 + t71 * t225) * qJD(1)) * t411 + (qJD(1) * t289 + t224 * t8 - t225 * t7) * t420 - t12 * t339 / 0.2e1 + ((-t342 * t376 - t348) * t225 + (t246 + (t254 * t224 + (t375 - t430) * t225) * qJD(3)) * t224) * t316 + ((-t343 * t375 + t348) * t224 + (t246 + (-t430 * t225 + (t376 + t254) * t224) * qJD(3)) * t225) * t319 - qJD(1) * ((t352 * t230 + t351 * t233) * qJD(1) + ((t224 * t356 - t225 * t357) * t233 + (t224 * t358 + t225 * t359) * t230) * qJD(3)) / 0.2e1 - t436 * t338 / 0.2e1 + ((-t14 * t354 + t35 * t362 + t9 * t391 + t29 * (t109 + t45) + (t29 * t392 - t354 * t36) * qJD(1)) * t225 + (-t13 * t354 + t36 * t362 + t9 * t392 + t29 * (t110 + t46) + (-t29 * t391 + t35 * t354) * qJD(1)) * t224 - t35 * (qJD(1) * t165 + t124 * t214 - t164 * t171 - t208 * t342) - t36 * (-qJD(1) * t167 - t125 * t214 - t164 * t170 - t208 * t343) - t29 * (-t124 * t170 - t125 * t171 - t165 * t343 - t167 * t342) - ((-t35 * t88 + t36 * t90) * t230 + t240 * t233) * qJD(4)) * m(5) + (-(t161 * t65 - t405) * qJD(1) - (t64 * (-t161 * t224 - t162 * t225) + t288 * t297) * qJD(3) + t30 * t277 + t64 * t243 + t288 * t184 + (-t48 * t224 - t49 * t225 + (-t225 * t66 + t400) * qJD(1)) * t197) * m(4) + (qJD(1) * t37 + t1 + (-t432 * t443 + t433 * t224 ^ 2 + (t224 * t54 + t55 * t225) * qJD(1)) * t440) * t224 / 0.2e1 - (qJD(1) * t38 + t2 + (t432 * t225 ^ 2 - t433 * t443 + (t224 * t52 + t53 * t225) * qJD(1)) * t440) * t225 / 0.2e1 + (t10 + t266 + t22) * t347 / 0.2e1 + (t11 + t265 + t23) * t320; -t11 * t332 / 0.2e1 + t1 * t368 / 0.2e1 + (t230 * t291 - t233 * t51) * t426 + ((qJD(3) * t291 - t15) * t233 + (-qJD(1) * t292 + qJD(3) * t51 + t224 * t3 + t225 * t4) * t230) * t424 + t230 * t10 * t320 + t2 * t373 / 0.2e1 + (t230 * t293 - t233 * t50) * t427 + ((qJD(3) * t293 - t16) * t233 + (-qJD(1) * t294 + qJD(3) * t50 + t224 * t5 + t225 * t6) * t230) * t423 + t12 * t315 - t233 * (t395 + t396 + t410 + t412 - t413) / 0.2e1 + (t230 * t289 - t233 * t57) * t301 + ((qJD(3) * t289 - t24) * t233 + (-qJD(1) * t290 + qJD(3) * t57 + t224 * t7 + t225 * t8) * t230) * t420 + (t147 * t428 + t244 * t148 - t225 * t253) * t425 + (-t145 * t428 + t146 * t244 - t224 * t253) * t422 + (t255 * t233 + (-t229 * t428 + t232 * t244) * t230) * t421 + t436 * t340 / 0.2e1 + ((qJD(3) * t240 - t13 * t90 + t14 * t88 + t35 * t46 - t36 * t45) * t233 + (t35 * (-qJD(3) * t88 + t114 * t224) + t36 * (qJD(3) * t90 - t114 * t225) + t9 * t287 + t29 * (-t224 * t45 + t225 * t46 - t346 * t90 - t347 * t88) + (-t13 * t225 + t14 * t224 + (t224 * t36 + t398) * qJD(1)) * t163) * t230 - t35 * (-t107 * t214 - t171 * t177) - t36 * (t108 * t214 - t170 * t177) - t29 * (t107 * t170 + t108 * t171)) * m(5);];
tauc = t17(:);
