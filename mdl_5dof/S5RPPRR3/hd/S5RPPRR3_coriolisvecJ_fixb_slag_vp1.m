% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:28
% EndTime: 2020-01-03 11:28:05
% DurationCPUTime: 13.66s
% Computational Cost: add. (13681->564), mult. (9868->720), div. (0->0), fcn. (7604->10), ass. (0->321)
t408 = rSges(4,2) * sin(pkin(9));
t239 = cos(pkin(9));
t412 = rSges(4,1) * t239;
t457 = t408 - t412;
t237 = qJ(1) + pkin(8);
t228 = cos(t237);
t333 = qJD(3) * t228;
t334 = qJD(1) * t228;
t226 = sin(t237);
t335 = qJD(1) * t226;
t344 = pkin(2) * t334 + qJ(3) * t335;
t456 = -t333 + t344;
t235 = pkin(9) + qJ(4);
t229 = qJ(5) + t235;
t220 = sin(t229);
t221 = cos(t229);
t374 = t220 * t228;
t187 = Icges(6,4) * t374;
t372 = t221 * t228;
t391 = Icges(6,5) * t226;
t111 = Icges(6,1) * t372 - t187 + t391;
t236 = qJD(4) + qJD(5);
t169 = t226 * t236;
t367 = t228 * t236;
t213 = Icges(6,4) * t221;
t283 = -Icges(6,2) * t220 + t213;
t440 = Icges(6,1) * t220 + t213;
t447 = t440 + t283;
t393 = Icges(6,4) * t220;
t165 = Icges(6,1) * t221 - t393;
t110 = -Icges(6,5) * t228 + t165 * t226;
t162 = Icges(6,2) * t221 + t393;
t450 = -t162 * t226 + t110;
t247 = qJD(1) * t447 + t169 * (-Icges(6,2) * t372 + t111 - t187) - t367 * t450;
t449 = t162 - t165;
t387 = Icges(6,6) * t226;
t109 = Icges(6,4) * t372 - Icges(6,2) * t374 + t387;
t451 = t440 * t228 + t109;
t108 = -Icges(6,6) * t228 + t226 * t283;
t452 = t440 * t226 + t108;
t427 = qJD(1) * t449 + t169 * t451 - t367 * t452;
t454 = t247 * t220 + t221 * t427;
t225 = sin(t235);
t227 = cos(t235);
t394 = Icges(5,4) * t225;
t176 = Icges(5,1) * t227 - t394;
t119 = -Icges(5,5) * t228 + t176 * t226;
t370 = t225 * t228;
t195 = Icges(5,4) * t370;
t368 = t227 * t228;
t392 = Icges(5,5) * t226;
t120 = Icges(5,1) * t368 - t195 + t392;
t173 = Icges(5,2) * t227 + t394;
t138 = t173 * t226;
t256 = t226 * (-Icges(5,2) * t368 + t120 - t195) - t228 * (t119 - t138);
t216 = Icges(5,4) * t227;
t284 = -Icges(5,2) * t225 + t216;
t117 = -Icges(5,6) * t228 + t226 * t284;
t388 = Icges(5,6) * t226;
t118 = Icges(5,4) * t368 - Icges(5,2) * t370 + t388;
t439 = Icges(5,1) * t225 + t216;
t140 = t439 * t226;
t141 = t439 * t228;
t431 = t226 * (t118 + t141) - t228 * (t117 + t140);
t453 = -t256 * t225 - t227 * t431;
t448 = t457 * t226;
t350 = t439 + t284;
t351 = t173 - t176;
t446 = (t225 * t350 + t227 * t351) * qJD(1);
t445 = 2 * qJD(4);
t418 = pkin(4) * t225;
t161 = Icges(6,5) * t221 - Icges(6,6) * t220;
t106 = -Icges(6,3) * t228 + t161 * t226;
t382 = t109 * t220;
t282 = -t111 * t221 + t382;
t272 = -t106 + t282;
t444 = t367 * t272;
t240 = -pkin(6) - qJ(3);
t205 = t226 * t240;
t384 = qJ(3) * t226;
t222 = t239 * pkin(3) + pkin(2);
t414 = pkin(2) - t222;
t105 = t228 * t414 + t205 + t384;
t181 = pkin(2) * t228 + t384;
t168 = qJD(1) * t181;
t441 = -qJD(1) * t105 + t168;
t342 = t226 * rSges(3,1) + t228 * rSges(3,2);
t403 = rSges(4,3) * t228;
t437 = t403 + t448;
t405 = rSges(6,2) * t221;
t410 = rSges(6,1) * t220;
t166 = t405 + t410;
t132 = t166 * t226;
t133 = t166 * t228;
t406 = rSges(6,2) * t220;
t409 = rSges(6,1) * t221;
t167 = -t406 + t409;
t373 = t221 * t226;
t189 = rSges(6,1) * t373;
t375 = t220 * t226;
t401 = rSges(6,3) * t228;
t112 = -rSges(6,2) * t375 + t189 - t401;
t190 = rSges(6,2) * t374;
t324 = rSges(6,1) * t372;
t113 = rSges(6,3) * t226 - t190 + t324;
t206 = t228 * t240;
t346 = t226 * t222 + t206;
t417 = pkin(4) * t227;
t191 = t222 + t417;
t234 = -pkin(7) + t240;
t354 = t226 * t191 + t228 * t234;
t89 = -t346 + t354;
t153 = t228 * t191;
t90 = t222 * t228 + t226 * t234 - t153 - t205;
t35 = t112 * t169 + t113 * t367 + qJD(2) + (t226 * t89 - t228 * t90) * qJD(4);
t215 = qJD(3) * t226;
t219 = t226 * pkin(2);
t179 = -qJ(3) * t228 + t219;
t232 = sin(qJ(1)) * pkin(1);
t310 = t179 + t232;
t294 = -t179 + t346 + t310;
t330 = qJD(4) * t228;
t296 = t330 * t418;
t39 = t296 + t166 * t367 - t215 + (t112 + t294 + t89) * qJD(1);
t233 = cos(qJ(1)) * pkin(1);
t224 = qJD(1) * t233;
t332 = qJD(4) * t226;
t318 = t225 * t332;
t267 = -pkin(4) * t318 - t333;
t251 = -t169 * t166 + t224 + t267;
t348 = t181 - t105;
t397 = t113 - t90;
t40 = (t348 + t397) * qJD(1) + t251;
t436 = -t39 * (-qJD(1) * t132 + t167 * t367) - t35 * (-t132 * t169 - t133 * t367) - t40 * (-qJD(1) * t133 - t169 * t167);
t172 = Icges(5,5) * t227 - Icges(5,6) * t225;
t115 = -Icges(5,3) * t228 + t172 * t226;
t338 = qJD(1) * t115;
t61 = t117 * t227 + t119 * t225;
t82 = -qJD(4) * t138 + (t228 * t284 + t388) * qJD(1);
t84 = -qJD(4) * t140 + (t176 * t228 + t392) * qJD(1);
t435 = qJD(4) * t61 + t225 * t82 - t227 * t84 - t338;
t157 = t284 * qJD(4);
t158 = t176 * qJD(4);
t171 = Icges(5,5) * t225 + Icges(5,6) * t227;
t434 = -qJD(1) * t171 + qJD(4) * (t173 * t227 + t225 * t439) + t157 * t225 - t158 * t227;
t386 = Icges(5,3) * t226;
t116 = Icges(5,5) * t368 - Icges(5,6) * t370 + t386;
t280 = t118 * t227 + t120 * t225;
t81 = qJD(1) * t117 + t173 * t330;
t83 = qJD(1) * t119 + qJD(4) * t141;
t433 = -qJD(1) * t116 + qJD(4) * t280 - t225 * t81 + t227 * t83;
t160 = Icges(6,5) * t220 + Icges(6,6) * t221;
t298 = t449 * t236;
t299 = t447 * t236;
t430 = -qJD(1) * t160 + t220 * t299 + t221 * t298;
t385 = Icges(6,3) * t226;
t107 = Icges(6,5) * t372 - Icges(6,6) * t374 + t385;
t303 = -qJD(1) * t108 + t111 * t236 - t162 * t367;
t305 = qJD(1) * t110 + t236 * t451;
t429 = -qJD(1) * t107 + t220 * t303 + t221 * t305;
t304 = (t228 * t283 + t387) * qJD(1) + t450 * t236;
t306 = -(t165 * t228 + t391) * qJD(1) + t452 * t236;
t340 = qJD(1) * t106;
t428 = t220 * t304 + t221 * t306 - t340;
t244 = qJD(1) ^ 2;
t154 = qJD(1) * t169;
t426 = t154 / 0.2e1;
t155 = qJD(1) * t367;
t425 = -t155 / 0.2e1;
t424 = -t169 / 0.2e1;
t423 = t169 / 0.2e1;
t422 = t367 / 0.2e1;
t421 = -t367 / 0.2e1;
t420 = -t226 / 0.2e1;
t419 = -t228 / 0.2e1;
t416 = -qJD(1) / 0.2e1;
t415 = qJD(1) / 0.2e1;
t74 = t236 * t133 + (t167 * t226 - t401) * qJD(1);
t341 = t234 - t240;
t76 = t296 + (t341 * t228 + (t191 - t222) * t226) * qJD(1);
t413 = -t74 - t76;
t411 = rSges(5,1) * t227;
t407 = rSges(5,2) * t225;
t404 = rSges(4,3) * t226;
t402 = rSges(5,3) * t228;
t178 = rSges(5,1) * t225 + rSges(5,2) * t227;
t143 = t178 * t226;
t369 = t226 * t227;
t371 = t225 * t226;
t121 = rSges(5,1) * t369 - rSges(5,2) * t371 - t402;
t319 = t178 * t330;
t52 = t319 - t215 + (t121 + t294) * qJD(1);
t400 = t143 * t52;
t307 = -rSges(5,2) * t370 + rSges(5,3) * t226;
t326 = rSges(5,1) * t368;
t122 = t307 + t326;
t271 = -t178 * t332 - t333;
t266 = t224 + t271;
t53 = (t122 + t348) * qJD(1) + t266;
t399 = t228 * t53;
t398 = rSges(6,3) - t234;
t383 = t108 * t220;
t381 = t110 * t221;
t380 = t117 * t225;
t379 = t118 * t225;
t378 = t119 * t227;
t377 = t160 * t226;
t127 = t160 * t228;
t376 = t171 * t226;
t137 = t171 * t228;
t58 = -t162 * t374 + t372 * t440 + t377;
t366 = t58 * qJD(1);
t73 = -t173 * t370 + t368 * t439 + t376;
t365 = t73 * qJD(1);
t208 = qJ(3) * t334;
t343 = t208 + t215;
t142 = pkin(2) * t335 - t343;
t364 = -t208 - (-t226 * t414 + t206) * qJD(1) - t142;
t223 = t244 * t233;
t355 = qJD(1) * t456 + t223;
t349 = rSges(6,3) * t335 + qJD(1) * t324;
t347 = rSges(5,3) * t335 + qJD(1) * t326;
t345 = rSges(4,3) * t335 + t334 * t412;
t336 = qJD(1) * t172;
t331 = qJD(4) * t227;
t329 = (qJD(4) ^ 2) * t417;
t328 = t244 * t232;
t325 = t236 * t410;
t75 = -t226 * t325 + (-t169 * t221 - t220 * t334) * rSges(6,2) + t349;
t327 = t112 * t334 - t113 * t335 + t226 * t75;
t322 = qJD(4) * t418;
t183 = t222 * t334;
t320 = qJD(1) * (-t240 * t335 + t183 - t344) + t355;
t317 = t226 * t331;
t316 = t335 / 0.2e1;
t315 = -t334 / 0.2e1;
t313 = t332 / 0.2e1;
t311 = t330 / 0.2e1;
t309 = t166 + t418;
t302 = -t107 - t383;
t301 = -t107 + t381;
t297 = t224 - t333;
t295 = qJD(1) * t215 - t328;
t182 = rSges(3,1) * t228 - rSges(3,2) * t226;
t288 = -t407 + t411;
t287 = -t226 * t53 + t228 * t52;
t56 = -t109 * t221 - t111 * t220;
t281 = t378 - t380;
t279 = -t120 * t227 + t379;
t278 = t121 * t226 + t122 * t228;
t277 = -t162 * t220 + t221 * t440;
t275 = -t225 * t173 + t227 * t439;
t274 = pkin(2) - t457;
t144 = t178 * t228;
t270 = (-t226 * t39 - t228 * t40) * qJD(1);
t96 = t119 * t369;
t48 = -t115 * t228 - t117 * t371 + t96;
t97 = t120 * t369;
t49 = t116 * t228 + t118 * t371 - t97;
t269 = (-t226 * t49 - t228 * t48) * qJD(4);
t98 = t117 * t370;
t50 = -t115 * t226 - t119 * t368 + t98;
t51 = t116 * t226 - t279 * t228;
t268 = (-t226 * t51 - t228 * t50) * qJD(4);
t263 = -qJD(1) * t161 + t127 * t169 - t367 * t377;
t261 = qJD(1) * t282 - t127 * t236 - t340;
t260 = t236 * t377 + (-t161 * t228 + t381 - t383 - t385) * qJD(1);
t259 = qJD(1) * t279 - qJD(4) * t137 - t338;
t258 = qJD(4) * t376 + (-t172 * t228 + t281 - t386) * qJD(1);
t255 = qJD(1) * t277 - t161 * t236;
t254 = t275 * qJD(1) - t172 * qJD(4);
t10 = t261 * t226 - t228 * t429;
t11 = -t226 * t428 + t260 * t228;
t12 = t226 * t429 + t261 * t228;
t91 = t110 * t373;
t43 = -t106 * t228 - t108 * t375 + t91;
t92 = t111 * t373;
t44 = t107 * t228 + t109 * t375 - t92;
t57 = t226 * t277 - t127;
t54 = t57 * qJD(1);
t20 = -t169 * t44 - t367 * t43 + t54;
t93 = t108 * t374;
t45 = -t106 * t226 - t110 * t372 + t93;
t46 = t107 * t226 - t228 * t282;
t21 = -t169 * t46 - t367 * t45 - t366;
t28 = t255 * t226 + t228 * t430;
t29 = -t226 * t430 + t255 * t228;
t30 = -t220 * t306 + t221 * t304;
t31 = t220 * t305 - t221 * t303;
t55 = t108 * t221 + t110 * t220;
t9 = t260 * t226 + t228 * t428;
t253 = (qJD(1) * t28 - t10 * t169 + t154 * t45 - t155 * t46 - t367 * t9) * t420 + (-t220 * t427 + t221 * t247) * t416 + t20 * t316 + t21 * t315 + (qJD(1) * t29 - t11 * t367 - t12 * t169 + t154 * t43 - t155 * t44) * t419 + (-t226 * t44 - t228 * t43) * t426 + (-t226 * t46 - t228 * t45) * t425 + (-t10 * t226 - t228 * t9 + (t226 * t45 - t228 * t46) * qJD(1)) * t424 + (-t11 * t228 - t12 * t226 + (t226 * t43 - t228 * t44) * qJD(1)) * t421 + (-t226 * t31 - t228 * t30 + (t226 * t55 - t228 * t56) * qJD(1)) * t415 + (t263 * t226 + t228 * t454) * t423 + (-t226 * t454 + t263 * t228) * t422;
t252 = -t236 * t405 - t322 - t325;
t85 = qJD(4) * t144 + (t226 * t288 - t402) * qJD(1);
t86 = -rSges(5,1) * t318 + (-t225 * t334 - t317) * rSges(5,2) + t347;
t249 = t226 * t86 - t228 * t85 + (t121 * t228 - t122 * t226) * qJD(1);
t159 = t288 * qJD(4);
t151 = t191 * t334;
t148 = t167 * t236;
t125 = -t228 * t457 + t404;
t102 = t226 * t112;
t88 = (t125 + t181) * qJD(1) + t297;
t77 = t151 - t183 + (-qJD(1) * t341 - t322) * t226;
t72 = t226 * t275 - t137;
t65 = t72 * qJD(1);
t64 = ((-qJD(1) * t408 - qJD(3)) * t228 + t345) * qJD(1) + t355;
t63 = (qJD(1) * t437 - t142) * qJD(1) + t295;
t59 = qJD(4) * t278 + qJD(2);
t42 = t159 * t330 + (t271 + t86) * qJD(1) + t320;
t41 = -t159 * t332 + (-t85 - t319 + t364) * qJD(1) + t295;
t37 = -t226 * t434 + t254 * t228;
t36 = t254 * t226 + t228 * t434;
t34 = qJD(4) * t279 + t225 * t83 + t227 * t81;
t33 = qJD(4) * t281 + t225 * t84 + t227 * t82;
t32 = t249 * qJD(4);
t25 = t228 * t329 + t148 * t367 - t154 * t166 + (t267 + t75 + t77) * qJD(1) + t320;
t24 = -t226 * t329 - t148 * t169 - t155 * t166 + (-t296 + t364 + t413) * qJD(1) + t295;
t23 = t268 - t365;
t22 = t65 + t269;
t8 = t112 * t155 - t113 * t154 + t169 * t75 - t367 * t74 + (t226 * t77 - t228 * t76 + (t226 * t90 + t228 * t89) * qJD(1)) * qJD(4);
t1 = [t155 * t58 / 0.2e1 + (-t298 * t220 + t299 * t221) * qJD(1) + (t275 * qJD(4) + t157 * t227 + t158 * t225) * qJD(1) + (t54 - (t226 * t302 + t46 + t91) * t367 + (t45 + t92 - t93 + (t106 - t382) * t226) * t169 + (t169 * t301 - t444) * t228) * t423 + m(3) * ((-t244 * t342 - t328) * (t182 + t233) + (-t223 + (-0.2e1 * rSges(3,1) * t334 + 0.2e1 * rSges(3,2) * t335 + qJD(1) * t182) * qJD(1)) * (-t342 - t232)) + t56 * t425 + (t65 + ((t50 + t97 - t98 + (t115 - t379) * t226) * t226 + (-t96 - t51 + (t115 - t279) * t228 + (t378 + t380) * t226) * t228) * qJD(4)) * t313 + (t55 + t57) * t426 + (t366 - (t44 - t93) * t367 + (t43 - t91) * t169 + (-t169 * t272 - t301 * t367) * t228 + (-t169 * t302 + t444) * t226 + t21) * t422 + (t29 + t30) * t421 + (t23 + t365 + ((-t49 + t98 + (t116 - t378) * t228) * t228 + (t48 - t96 + (t116 + t380) * t226) * t226) * qJD(4)) * t311 + (t24 * (t153 - t190 + t233) + t40 * t215 + t25 * (t189 + t232 + t354) + t39 * (t151 + t224 + t349) + (-t25 * rSges(6,3) - t39 * qJD(3) + t24 * t409 + t40 * t252) * t228 + (t24 * t398 - t25 * t406 + t39 * t252) * t226 + (-t40 * t232 + (-t39 * t406 + t398 * t40) * t228 + (t40 * (-t167 - t191) - t39 * t234) * t226) * qJD(1) - (qJD(1) * t397 + t251 - t40 + t441) * t39) * m(6) + (t41 * (-t205 + t233 + (t222 + t411) * t228 + t307) + t53 * t215 + t42 * (t121 + t232 + t346) + t52 * (t183 + t297 + t347) + (-t178 * t399 - t400) * qJD(4) + (-t53 * t232 + (t53 * (rSges(5,3) - t240) - t52 * t407) * t228 + (t53 * (-t222 - t288) - t52 * t240) * t226) * qJD(1) - (qJD(1) * t122 + t266 + t441 - t53) * t52) * m(5) + (t63 * (t233 + t384 + t404) + t64 * (t219 + t232 - t448) + (t63 * t274 + t64 * (-rSges(4,3) - qJ(3))) * t228 + (t343 + (-t226 * t274 - t232 + t403) * qJD(1)) * t88 + (-t168 + t224 - t297 + t345 + t88 + (-t228 * t408 - t125) * qJD(1) + t456) * (-t215 + (t310 - t437) * qJD(1))) * m(4) + (t31 + t28 + t20) * t424 - (t34 + t36 + t22) * t332 / 0.2e1 - (-t280 * qJD(1) + t33 + t37) * t330 / 0.2e1 + (t228 * t73 + (t61 + t72) * t226) * qJD(4) * t415; m(5) * t32 + m(6) * t8; m(4) * (-t226 * t64 - t228 * t63) + m(5) * (-t226 * t42 - t228 * t41) + m(6) * (-t226 * t25 - t228 * t24); ((-t330 * t376 - t336) * t228 + (-t446 + (t228 * t137 + t453) * qJD(4)) * t226) * t311 + ((-t225 * t351 + t227 * t350) * qJD(1) + (-t225 * t431 + t227 * t256) * qJD(4)) * t416 + t253 + ((t137 * t332 - t336) * t226 + (t446 + (-t226 * t376 - t453) * qJD(4)) * t228) * t313 + (-t226 * t34 - t228 * t33 + (t226 * t61 + t228 * t280) * qJD(1)) * t415 + (qJD(1) * t36 + (-(t258 * t226 + t228 * t435) * t228 - (t259 * t226 - t228 * t433) * t226 + (t50 * t226 - t51 * t228) * qJD(1)) * t445) * t420 + (qJD(1) * t37 + (-(-t226 * t435 + t258 * t228) * t228 - (t226 * t433 + t259 * t228) * t226 + (t48 * t226 - t49 * t228) * qJD(1)) * t445) * t419 + (t269 + t22) * t316 + (t268 + t23) * t315 + (t8 * t102 + t35 * t327 + (t8 * t89 + t35 * t77 - t24 * t309 + t40 * (-pkin(4) * t331 - t148) + (-t309 * t39 + t35 * t90) * qJD(1)) * t226 + (t8 * t397 + t35 * t413 + t25 * t309 + t39 * t148 + (-t309 * t40 + t35 * t89) * qJD(1)) * t228 - (-t40 * t317 + (t35 * (-t226 ^ 2 - t228 ^ 2) * qJD(4) + t270) * t225) * pkin(4) + t436) * m(6) + (-(-t144 * t53 - t400) * qJD(1) - (t59 * (-t143 * t226 - t144 * t228) + t287 * t288) * qJD(4) + t32 * t278 + t59 * t249 + t287 * t159 + (-t41 * t226 + t42 * t228 + (-t226 * t52 - t399) * qJD(1)) * t178) * m(5); t253 + (t8 * (t113 * t228 + t102) + t35 * (-t228 * t74 + t327) + (-t226 * t40 + t228 * t39) * t148 + (-t24 * t226 + t25 * t228 + t270) * t166 + t436) * m(6);];
tauc = t1(:);
