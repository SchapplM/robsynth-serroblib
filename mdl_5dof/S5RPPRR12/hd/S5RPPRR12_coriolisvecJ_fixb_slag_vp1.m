% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR12_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:07:16
% DurationCPUTime: 14.64s
% Computational Cost: add. (13116->733), mult. (19243->1002), div. (0->0), fcn. (17681->8), ass. (0->365)
t268 = pkin(8) + qJ(4);
t253 = sin(t268);
t274 = cos(qJ(5));
t275 = cos(qJ(1));
t421 = t274 * t275;
t272 = sin(qJ(5));
t273 = sin(qJ(1));
t423 = t273 * t272;
t195 = -t253 * t423 + t421;
t422 = t273 * t274;
t424 = t272 * t275;
t196 = t253 * t422 + t424;
t254 = cos(t268);
t428 = t254 * t273;
t103 = Icges(6,5) * t196 + Icges(6,6) * t195 - Icges(6,3) * t428;
t197 = t253 * t424 + t422;
t198 = t253 * t421 - t423;
t427 = t254 * t275;
t105 = -Icges(6,5) * t198 + Icges(6,6) * t197 + Icges(6,3) * t427;
t179 = Icges(6,4) * t198;
t108 = Icges(6,2) * t197 + Icges(6,6) * t427 - t179;
t178 = Icges(6,4) * t197;
t110 = Icges(6,1) * t198 - Icges(6,5) * t427 - t178;
t320 = -t108 * t197 - t110 * t198;
t442 = Icges(6,4) * t196;
t106 = Icges(6,2) * t195 - Icges(6,6) * t428 + t442;
t177 = Icges(6,4) * t195;
t109 = Icges(6,1) * t196 - Icges(6,5) * t428 + t177;
t458 = t195 * t106 + t196 * t109;
t513 = t320 + t458 + (-t103 * t273 - t105 * t275) * t254;
t384 = qJD(5) * t273;
t386 = qJD(4) * t275;
t206 = -t254 * t384 + t386;
t234 = qJD(5) * t253 + qJD(1);
t27 = t103 * t427 + t197 * t106 - t198 * t109;
t322 = Icges(6,5) * t274 - Icges(6,6) * t272;
t149 = Icges(6,3) * t253 + t254 * t322;
t440 = Icges(6,4) * t274;
t324 = -Icges(6,2) * t272 + t440;
t151 = Icges(6,6) * t253 + t254 * t324;
t441 = Icges(6,4) * t272;
t326 = Icges(6,1) * t274 - t441;
t153 = Icges(6,5) * t253 + t254 * t326;
t53 = t149 * t427 + t151 * t197 - t153 * t198;
t512 = -t206 * t27 - t234 * t53;
t511 = pkin(1) * t273;
t510 = t275 * pkin(1);
t271 = -pkin(6) - qJ(3);
t419 = qJ(3) + t271;
t269 = sin(pkin(8));
t425 = t269 * t275;
t191 = pkin(3) * t425 + t273 * t419;
t255 = qJD(3) * t275;
t259 = t275 * qJ(2);
t226 = -t259 + t511;
t256 = qJD(2) * t273;
t394 = -qJD(1) * t226 + t256;
t371 = t255 + t394;
t389 = qJD(1) * t275;
t249 = qJ(2) * t389;
t395 = t255 + t256;
t372 = t249 + t395;
t367 = t269 * t389;
t390 = qJD(1) * t273;
t399 = pkin(3) * t367 + t271 * t390;
t508 = -qJD(1) * t191 - t371 + t372 + t399;
t338 = rSges(6,1) * t198 - rSges(6,2) * t197;
t114 = rSges(6,3) * t427 - t338;
t337 = rSges(6,1) * t274 - rSges(6,2) * t272;
t155 = rSges(6,3) * t253 + t254 * t337;
t383 = qJD(5) * t275;
t387 = qJD(4) * t273;
t205 = t254 * t383 + t387;
t464 = t254 * pkin(4);
t216 = pkin(7) * t253 + t464;
t507 = -t114 * t234 + t155 * t205 + t216 * t387;
t506 = qJD(3) * t273;
t319 = t108 * t272 + t110 * t274;
t40 = t105 * t253 - t254 * t319;
t26 = -t105 * t428 + t195 * t108 - t110 * t196;
t362 = t254 * t387;
t502 = t253 * t389 + t362;
t231 = Icges(5,4) * t428;
t430 = t253 * t273;
t439 = Icges(5,5) * t275;
t161 = Icges(5,1) * t430 + t231 + t439;
t443 = Icges(5,4) * t254;
t327 = Icges(5,1) * t253 + t443;
t162 = -Icges(5,5) * t273 + t275 * t327;
t210 = -Icges(5,2) * t253 + t443;
t173 = t210 * t275;
t285 = t273 * (t162 + t173) - t275 * (-Icges(5,2) * t430 + t161 + t231);
t444 = Icges(5,4) * t253;
t325 = Icges(5,2) * t254 + t444;
t159 = Icges(5,6) * t275 + t273 * t325;
t160 = -Icges(5,6) * t273 + t275 * t325;
t212 = Icges(5,1) * t254 - t444;
t175 = t212 * t273;
t176 = t212 * t275;
t286 = t273 * (t160 - t176) - t275 * (t159 - t175);
t501 = -t286 * t253 + t285 * t254;
t401 = t210 + t327;
t402 = -t325 + t212;
t500 = (t253 * t401 - t254 * t402) * qJD(1);
t497 = 0.2e1 * qJD(4);
t52 = -t149 * t428 + t151 * t195 + t153 * t196;
t496 = t205 * t26 + t52 * t234;
t148 = Icges(6,3) * t254 - t253 * t322;
t317 = t151 * t272 - t153 * t274;
t321 = t106 * t272 - t109 * t274;
t277 = t205 * (-t149 * t275 + t319) + t206 * (t149 * t273 + t321) + t234 * (t148 + t317);
t495 = t277 * t254;
t314 = t160 * t254 + t162 * t253;
t492 = t314 * t275;
t229 = t273 * qJ(2) + t510;
t434 = qJ(3) * t275;
t426 = t269 * t273;
t456 = rSges(4,2) * cos(pkin(8));
t489 = -rSges(4,1) * t426 - t275 * rSges(4,3) - t273 * t456;
t491 = t229 + t434 - t489;
t351 = -rSges(3,2) * t275 + t273 * rSges(3,3);
t490 = t229 + t351;
t244 = pkin(3) * t426;
t373 = t244 + t229;
t323 = Icges(5,5) * t253 + Icges(5,6) * t254;
t158 = -Icges(5,3) * t273 + t275 * t323;
t392 = qJD(1) * t158;
t77 = t160 * t253 - t162 * t254;
t95 = qJD(1) * t159 - qJD(4) * t173;
t97 = -qJD(4) * t176 + (t273 * t327 + t439) * qJD(1);
t488 = qJD(4) * t77 + t253 * t97 + t254 * t95 + t392;
t348 = qJD(1) * t253 + qJD(5);
t363 = t254 * t386;
t487 = t273 * t348 - t363;
t200 = t325 * qJD(4);
t201 = t327 * qJD(4);
t208 = Icges(5,5) * t254 - Icges(5,6) * t253;
t486 = qJD(1) * t208 + qJD(4) * (t210 * t253 - t212 * t254) + t200 * t254 + t201 * t253;
t315 = t159 * t253 - t161 * t254;
t157 = Icges(5,3) * t275 + t273 * t323;
t393 = qJD(1) * t157;
t96 = qJD(1) * t160 + t210 * t387;
t98 = qJD(1) * t162 + qJD(4) * t175;
t485 = qJD(4) * t315 - t253 * t98 - t254 * t96 + t393;
t171 = (-Icges(6,2) * t274 - t441) * t254;
t282 = t205 * (Icges(6,2) * t198 - t110 + t178) + t206 * (-Icges(6,2) * t196 + t109 + t177) + t234 * (t153 + t171);
t174 = (-Icges(6,1) * t272 - t440) * t254;
t483 = t205 * (-Icges(6,1) * t197 + t108 - t179) + t206 * (-Icges(6,1) * t195 + t106 + t442) + t234 * (t151 - t174);
t25 = -t103 * t428 + t458;
t10 = t206 * t25 + t496;
t482 = -t10 / 0.2e1;
t380 = qJD(4) * qJD(5);
t361 = t253 * t380;
t141 = -qJD(1) * t205 + t273 * t361;
t481 = t141 / 0.2e1;
t142 = qJD(1) * t206 - t275 * t361;
t480 = t142 / 0.2e1;
t479 = -t205 / 0.2e1;
t478 = t205 / 0.2e1;
t477 = -t206 / 0.2e1;
t476 = t206 / 0.2e1;
t475 = -t234 / 0.2e1;
t474 = t234 / 0.2e1;
t473 = t253 / 0.2e1;
t472 = t273 / 0.2e1;
t471 = -t275 / 0.2e1;
t469 = rSges(3,2) - pkin(1);
t468 = rSges(6,3) + pkin(7);
t467 = pkin(3) * t269;
t466 = pkin(4) * t253;
t465 = pkin(7) * t254;
t364 = t253 * t387;
t368 = t254 * t389;
t296 = t364 - t368;
t90 = -t234 * t422 + (-t275 * t348 - t362) * t272;
t388 = qJD(4) * t254;
t91 = t348 * t421 + (-t234 * t272 + t274 * t388) * t273;
t45 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t296;
t47 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t296;
t49 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t296;
t7 = (qJD(4) * t321 + t45) * t253 + (qJD(4) * t103 - t272 * t47 + t274 * t49 + (-t106 * t274 - t109 * t272) * qJD(5)) * t254;
t463 = t7 * t206;
t295 = -t253 * t386 - t254 * t390;
t311 = t275 * t234;
t88 = -t272 * t487 + t274 * t311;
t89 = t272 * t311 + t274 * t487;
t44 = Icges(6,5) * t89 + Icges(6,6) * t88 + Icges(6,3) * t295;
t46 = Icges(6,4) * t89 + Icges(6,2) * t88 + Icges(6,6) * t295;
t48 = Icges(6,1) * t89 + Icges(6,4) * t88 + Icges(6,5) * t295;
t8 = (qJD(4) * t319 + t44) * t253 + (qJD(4) * t105 - t272 * t46 + t274 * t48 + (-t108 * t274 + t110 * t272) * qJD(5)) * t254;
t462 = t8 * t205;
t461 = -qJD(1) / 0.2e1;
t460 = -pkin(1) + t271;
t168 = (-Icges(6,5) * t272 - Icges(6,6) * t274) * t254;
t85 = qJD(4) * t148 + qJD(5) * t168;
t150 = Icges(6,6) * t254 - t253 * t324;
t86 = qJD(4) * t150 + qJD(5) * t171;
t152 = Icges(6,5) * t254 - t253 * t326;
t87 = qJD(4) * t152 + qJD(5) * t174;
t18 = (qJD(4) * t317 + t85) * t253 + (qJD(4) * t149 - t272 * t86 + t274 * t87 + (-t151 * t274 - t153 * t272) * qJD(5)) * t254;
t360 = t254 * t380;
t55 = t149 * t253 - t254 * t317;
t459 = t18 * t234 + t55 * t360;
t455 = rSges(5,2) * t253;
t454 = rSges(3,3) * t275;
t452 = rSges(6,3) * t254;
t214 = rSges(5,1) * t254 - t455;
t182 = t214 * t275;
t265 = t275 * rSges(5,3);
t339 = rSges(5,1) * t253 + rSges(5,2) * t254;
t101 = -qJD(4) * t182 + (t273 * t339 + t265) * qJD(1);
t203 = t339 * qJD(4);
t382 = qJD(1) * qJD(2);
t248 = t275 * t382;
t381 = qJD(1) * qJD(3);
t433 = qJ(3) * qJD(1) ^ 2;
t310 = -0.2e1 * t273 * t381 - t275 * t433 + t248;
t366 = t214 * t386;
t257 = qJD(2) * t275;
t202 = qJD(1) * t229 - t257;
t241 = t271 * t389;
t405 = t241 - (t244 - t434) * qJD(1) - t202;
t42 = -t203 * t387 + (-t101 + t366 + t405) * qJD(1) + t310;
t449 = t273 * t42;
t376 = qJD(4) * t455;
t308 = -rSges(5,3) * qJD(1) - t376;
t375 = t502 * rSges(5,1) + rSges(5,2) * t368;
t102 = t273 * t308 + t375;
t190 = t214 * t387;
t397 = t249 + t256;
t403 = qJD(1) * (-pkin(1) * t390 + t397) + t273 * t382;
t299 = -t273 * t433 + 0.2e1 * t275 * t381 + t403;
t292 = qJD(1) * (qJ(3) * t390 + t399) + t299;
t41 = t203 * t386 + (t102 + t190) * qJD(1) + t292;
t448 = t275 * t41;
t39 = t103 * t253 - t254 * t321;
t447 = t39 * t141;
t446 = t40 * t142;
t215 = t465 - t466;
t204 = qJD(4) * t215;
t154 = -t253 * t337 + t452;
t181 = (-rSges(6,1) * t272 - rSges(6,2) * t274) * t254;
t92 = qJD(4) * t154 + qJD(5) * t181;
t445 = t204 + t92;
t435 = qJ(3) * t273;
t431 = t208 * t275;
t429 = t253 * t275;
t169 = t273 * t208;
t313 = t210 * t254 + t212 * t253;
t75 = t275 * t313 - t169;
t420 = t75 * qJD(1);
t404 = t196 * rSges(6,1) + t195 * rSges(6,2);
t112 = -rSges(6,3) * t428 + t404;
t235 = pkin(4) * t430;
t186 = -pkin(7) * t428 + t235;
t414 = -t112 - t186;
t236 = pkin(7) * t427;
t188 = pkin(4) * t429 - t236;
t413 = t114 - t188;
t410 = t155 + t216;
t377 = t275 * t456;
t400 = rSges(4,1) * t367 + qJD(1) * t377;
t398 = t241 + t257;
t396 = rSges(3,2) * t390 + rSges(3,3) * t389;
t391 = qJD(1) * t323;
t385 = qJD(5) * t254;
t379 = -rSges(4,3) - pkin(1) - qJ(3);
t378 = t91 * rSges(6,1) + t90 * rSges(6,2) + rSges(6,3) * t364;
t59 = t275 * t157 + t159 * t428 + t161 * t430;
t60 = -t275 * t158 - t160 * t428 - t162 * t430;
t374 = t502 * pkin(4) + pkin(7) * t364;
t163 = rSges(5,1) * t430 + rSges(5,2) * t428 + t265;
t370 = t254 * t468;
t365 = t216 * t386;
t359 = -t390 / 0.2e1;
t357 = t388 / 0.2e1;
t356 = -t387 / 0.2e1;
t355 = t387 / 0.2e1;
t354 = -t386 / 0.2e1;
t350 = -t226 - t435;
t349 = -t257 + t506;
t345 = qJD(5) * t357;
t344 = rSges(6,1) * t89 + rSges(6,2) * t88;
t343 = t191 + t350;
t342 = -t275 * t419 + t373 + t434;
t340 = rSges(4,1) * t269 + t456;
t116 = -pkin(7) * t368 + t374;
t51 = -rSges(6,3) * t368 + t378;
t13 = qJD(1) * t116 - t141 * t155 - t206 * t92 + t234 * t51 + (t112 * t385 - t204 * t275 + t216 * t390) * qJD(4) + t292;
t115 = t295 * pkin(7) + (t253 * t390 - t363) * pkin(4);
t50 = rSges(6,3) * t295 + t344;
t14 = t142 * t155 + t205 * t92 - t234 * t50 + (-t114 * t385 + t204 * t273) * qJD(4) + (-t115 + t365 + t405) * qJD(1) + t310;
t336 = t13 * t273 + t14 * t275;
t335 = t25 * t275 + t26 * t273;
t334 = t25 * t273 - t26 * t275;
t28 = t105 * t427 - t320;
t333 = t27 * t275 + t273 * t28;
t332 = t27 * t273 - t275 * t28;
t331 = t273 * t40 + t275 * t39;
t330 = t273 * t39 - t275 * t40;
t263 = t273 * rSges(5,3);
t164 = t275 * t339 - t263;
t63 = t190 + (t164 + t343) * qJD(1) + t395;
t64 = -t366 + (t163 + t342) * qJD(1) + t349;
t329 = t273 * t63 - t275 * t64;
t318 = t112 * t275 + t114 * t273;
t316 = t159 * t254 + t161 * t253;
t301 = (t273 * t60 + t275 * t59) * qJD(4);
t145 = t273 * t157;
t61 = -t275 * t316 + t145;
t62 = -t158 * t273 + t492;
t300 = (t273 * t62 + t275 * t61) * qJD(4);
t78 = (-t163 * t273 - t164 * t275) * qJD(4);
t298 = -t452 + t466 + t467;
t297 = t339 + t467;
t291 = t103 * t206 + t105 * t205 + t149 * t234;
t290 = (Icges(6,5) * t195 - Icges(6,6) * t196) * t206 + (Icges(6,5) * t197 + Icges(6,6) * t198) * t205 + t168 * t234;
t288 = -qJD(1) * t314 - qJD(4) * t431 + t393;
t287 = qJD(1) * t316 + qJD(4) * t169 + t392;
t284 = t313 * qJD(1) - t323 * qJD(4);
t32 = (t188 + t343) * qJD(1) + t395 + t507;
t33 = -t365 + t112 * t234 - t155 * t206 + (t186 + t342) * qJD(1) + t349;
t36 = -t112 * t205 + t114 * t206 + (-t186 * t273 - t188 * t275) * qJD(4);
t278 = t36 * t318 + (-t273 * t33 - t275 * t32) * t155;
t227 = rSges(3,2) * t273 + t454;
t189 = t216 * t275;
t187 = t216 * t273;
t180 = t214 * t273;
t167 = -rSges(4,3) * t273 + t275 * t340;
t144 = qJD(1) * t490 - t257;
t143 = t256 + (-t226 + t227) * qJD(1);
t134 = t155 * t275;
t133 = t155 * t273;
t132 = t153 * t275;
t131 = t153 * t273;
t130 = t151 * t275;
t129 = t151 * t273;
t126 = rSges(6,1) * t197 + rSges(6,2) * t198;
t125 = rSges(6,1) * t195 - rSges(6,2) * t196;
t118 = t248 + (-qJD(1) * t351 - t202) * qJD(1);
t117 = qJD(1) * t396 + t403;
t100 = qJD(1) * t491 + t349;
t99 = (t167 + t350) * qJD(1) + t395;
t74 = t273 * t313 + t431;
t73 = (qJD(1) * t489 - t202) * qJD(1) + t310;
t72 = qJD(1) * (-rSges(4,3) * t390 + t400) + t299;
t71 = t74 * qJD(1);
t38 = qJD(4) * t314 - t253 * t95 + t254 * t97;
t37 = -qJD(4) * t316 - t253 * t96 + t254 * t98;
t35 = -t273 * t486 + t284 * t275;
t34 = t284 * t273 + t275 * t486;
t24 = t300 - t420;
t23 = t71 + t301;
t16 = t149 * t296 + t151 * t90 + t153 * t91 + t195 * t86 + t196 * t87 - t428 * t85;
t15 = t149 * t295 + t151 * t88 + t153 * t89 + t197 * t86 - t198 * t87 + t427 * t85;
t12 = t205 * t40 + t206 * t39 + t234 * t55;
t11 = t205 * t28 - t512;
t9 = -t112 * t142 + t114 * t141 - t205 * t51 + t206 * t50 + (t115 * t275 - t116 * t273 + (-t186 * t275 + t188 * t273) * qJD(1)) * qJD(4);
t6 = t105 * t296 + t108 * t90 - t110 * t91 + t195 * t46 + t196 * t48 - t428 * t44;
t5 = t103 * t296 + t106 * t90 + t109 * t91 + t195 * t47 + t196 * t49 - t428 * t45;
t4 = t105 * t295 + t108 * t88 - t110 * t89 + t197 * t46 - t198 * t48 + t427 * t44;
t3 = t103 * t295 + t106 * t88 + t109 * t89 + t197 * t47 - t198 * t49 + t427 * t45;
t2 = t141 * t25 + t142 * t26 + t16 * t234 + t205 * t6 + t206 * t5 + t360 * t52;
t1 = t141 * t27 + t142 * t28 + t15 * t234 + t205 * t4 + t206 * t3 + t360 * t53;
t17 = [t459 + (t71 + ((-t61 + t145 + t60) * t273 + (t62 - t492 + (t158 - t316) * t273 + t59) * t275) * qJD(4)) * t356 + t447 / 0.2e1 + t462 / 0.2e1 + t463 / 0.2e1 + t16 * t476 + ((t28 + t513) * t206 + t496) * t479 + t53 * t480 + t52 * t481 + (-qJD(4) * t313 + t200 * t253 - t201 * t254) * qJD(1) + t446 / 0.2e1 + (t15 + t10) * t478 + ((-t25 + t513) * t205 + t11 + t512) * t477 + (t420 + (t158 * t273 ^ 2 + (-t145 + t60 + (t158 + t316) * t275) * t275) * qJD(4) + t24) * t354 + (t14 * (-t236 + t259 + t338) + t13 * (t235 + t373 + t404) + (-t13 * t370 + t14 * t460) * t273 + (-t13 * t271 + t14 * t298) * t275 + (-t344 + t398 - t506 + (t253 * t468 + t464) * t386 + (-t510 + (-qJ(2) - t298 + t465) * t273) * qJD(1)) * t32 + (t374 + t378 + t32 + (-t370 * t275 - t188 + t435 - t511) * qJD(1) - t507 + t508) * t33) * m(6) + ((-t271 * t275 + t163 + t373) * t41 + (t273 * t460 + t275 * t297 + t259 - t263) * t42 + (t398 + (rSges(5,1) * t388 - pkin(1) * qJD(1) + t308) * t275 + (-qJD(3) + (-qJ(2) - t297) * qJD(1)) * t273) * t63 + (t375 + (-t376 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t273 - t190 + t63 - (t164 - t435) * qJD(1) + t508) * t64) * m(5) + (t73 * (rSges(4,1) * t425 + t259 + t377) + t99 * t257 + t72 * t491 + t100 * (t372 + t400) + (-t99 * qJD(3) + t73 * t379) * t273 + (t99 * t379 * t275 + (t99 * (-qJ(2) - t340) + t100 * t379) * t273) * qJD(1) - (-t99 + (t167 - t435) * qJD(1) + t371) * t100) * m(4) + (t118 * (t273 * t469 + t259 + t454) + t143 * t257 + t117 * t490 + t144 * (t396 + t397) + (t143 * t469 * t275 + (t143 * (-rSges(3,3) - qJ(2)) - t144 * pkin(1)) * t273) * qJD(1) - (qJD(1) * t227 - t143 + t394) * t144) * m(3) + (t38 + t34 + t23) * t355 + (qJD(1) * t77 + t35 + t37) * t386 / 0.2e1 + (t275 * t75 + (-t315 + t74) * t273) * qJD(4) * t461; 0.2e1 * (t13 * t471 + t14 * t472) * m(6) + 0.2e1 * (t449 / 0.2e1 - t448 / 0.2e1) * m(5) + 0.2e1 * (t471 * t72 + t472 * t73) * m(4) + 0.2e1 * (t117 * t471 + t118 * t472) * m(3); m(4) * (t273 * t72 + t275 * t73) + m(5) * (t273 * t41 + t275 * t42) + m(6) * t336; t11 * t383 * t473 + t253 * t384 * t482 + t331 * t345 + ((-t387 * t431 - t391) * t273 + (t500 + (t273 * t169 + t501) * qJD(4)) * t275) * t356 + ((t169 * t386 - t391) * t275 + (-t500 + (-t275 * t431 - t501) * qJD(4)) * t273) * t354 + qJD(1) * (t273 * t38 + t275 * t37 + (t273 * t315 + t77 * t275) * qJD(1)) / 0.2e1 - t12 * t385 / 0.2e1 + ((-t253 * t402 - t254 * t401) * qJD(1) + (t253 * t285 + t254 * t286) * qJD(4)) * t461 + (-qJD(1) * t330 + t273 * t8 + t275 * t7) * t474 + (((-t129 * t272 + t131 * t274 + t103) * t206 + (t130 * t272 - t132 * t274 + t105) * t205 + (-t150 * t272 + t152 * t274 + t149) * t234 + t55 * qJD(5)) * t254 + (qJD(5) * t330 + t277) * t253) * t475 + (-qJD(1) * t334 + t273 * t6 + t275 * t5) * t476 + ((t129 * t195 + t131 * t196) * t206 + (-t130 * t195 - t132 * t196) * t205 + (t150 * t195 + t152 * t196) * t234 + (t254 * t52 - t26 * t429) * qJD(5) + ((qJD(5) * t25 + t291) * t253 - t495) * t273) * t477 + (-qJD(1) * t332 + t273 * t4 + t275 * t3) * t478 + ((t129 * t197 - t131 * t198) * t206 + (-t130 * t197 + t132 * t198) * t205 + (t150 * t197 - t152 * t198) * t234 + (t254 * t53 + t27 * t430) * qJD(5) + ((-qJD(5) * t28 - t291) * t253 + t495) * t275) * t479 + t333 * t480 + t335 * t481 + ((-t13 * t410 - t33 * t445 + t9 * t413 + t36 * (t115 + t50) + (t32 * t410 + t36 * t414) * qJD(1)) * t275 + (t14 * t410 + t32 * t445 + t9 * t414 + t36 * (-t116 - t51) + (t33 * t410 - t36 * t413) * qJD(1)) * t273 - t32 * (qJD(1) * t189 + t134 * t234 + t154 * t205 + t215 * t387) - t33 * (qJD(1) * t187 + t133 * t234 - t154 * t206 - t215 * t386) - t36 * (-t133 * t205 - t134 * t206 - t187 * t387 - t189 * t386) - ((t112 * t33 - t114 * t32) * t254 + t278 * t253) * qJD(5)) * m(6) + (-(t180 * t64 + t182 * t63) * qJD(1) - (t78 * (-t180 * t273 - t182 * t275) - t329 * t339) * qJD(4) + 0.2e1 * t78 * (t101 * t275 - t102 * t273 + (-t163 * t275 + t164 * t273) * qJD(1)) - t329 * t203 + (t449 - t448 + (t273 * t64 + t275 * t63) * qJD(1)) * t214) * m(5) + (qJD(1) * t34 + t1 + ((t287 * t273 + t275 * t485) * t275 + (t288 * t273 - t275 * t488) * t273 + (-t61 * t273 + t62 * t275) * qJD(1)) * t497) * t472 + (qJD(1) * t35 + t2 + ((-t273 * t485 + t287 * t275) * t275 + (t273 * t488 + t288 * t275) * t273 + (-t59 * t273 + t60 * t275) * qJD(1)) * t497) * t275 / 0.2e1 + (t301 + t23 + t10) * t359 + (t300 + t24 + t11) * t389 / 0.2e1; t368 * t482 + t253 * t10 * t355 - t2 * t428 / 0.2e1 + (t253 * t52 - t254 * t334) * t481 + ((qJD(4) * t334 + t16) * t253 + (-qJD(1) * t335 + qJD(4) * t52 - t273 * t5 + t275 * t6) * t254) * t476 + t1 * t427 / 0.2e1 + (t253 * t53 - t254 * t332) * t480 + ((qJD(4) * t332 + t15) * t253 + (-qJD(1) * t333 + qJD(4) * t53 - t273 * t3 + t275 * t4) * t254) * t478 + t12 * t357 + (t446 + t447 + t459 + t462 + t463) * t473 + (t253 * t55 - t254 * t330) * t345 + ((qJD(4) * t330 + t18) * t253 + (-qJD(1) * t331 + qJD(4) * t55 - t273 * t7 + t275 * t8) * t254) * t474 + (t195 * t282 - t196 * t483 - t290 * t428) * t477 + (t282 * t197 + t198 * t483 + t290 * t427) * t479 + (t290 * t253 + (-t272 * t282 - t274 * t483) * t254) * t475 + (t253 * t354 + t254 * t359) * t11 + ((qJD(4) * t278 + t13 * t112 - t14 * t114 - t32 * t50 + t33 * t51) * t253 + (t32 * (-qJD(4) * t114 + t275 * t92) + t33 * (qJD(4) * t112 + t273 * t92) - t9 * t318 + t36 * (t112 * t390 - t114 * t389 - t273 * t50 - t275 * t51) + ((-t273 * t32 + t275 * t33) * qJD(1) + t336) * t155) * t254 - t32 * (-t126 * t234 + t181 * t205) - t33 * (t125 * t234 - t181 * t206) - t36 * (-t125 * t205 + t126 * t206)) * m(6);];
tauc = t17(:);
