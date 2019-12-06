% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:19
% EndTime: 2019-12-05 15:28:35
% DurationCPUTime: 10.90s
% Computational Cost: add. (9127->394), mult. (7869->509), div. (0->0), fcn. (6041->6), ass. (0->230)
t464 = Icges(6,4) + Icges(5,5);
t221 = pkin(7) + qJ(2);
t217 = sin(t221);
t219 = cos(t221);
t220 = pkin(8) + qJ(4);
t216 = sin(t220);
t204 = Icges(6,5) * t216;
t218 = cos(t220);
t268 = Icges(6,1) * t218 + t204;
t463 = -t217 * t268 + t464 * t219;
t339 = t216 * t217;
t177 = Icges(5,4) * t339;
t337 = t217 * t218;
t453 = Icges(5,1) * t337 - t177 - t463;
t266 = Icges(6,3) * t218 - t204;
t356 = Icges(5,4) * t216;
t462 = Icges(5,2) * t218 + t266 + t356;
t353 = Icges(6,5) * t218;
t155 = Icges(6,1) * t216 - t353;
t207 = Icges(5,4) * t218;
t461 = Icges(5,1) * t216 + t155 + t207;
t150 = Icges(5,5) * t218 - Icges(5,6) * t216;
t100 = Icges(5,3) * t217 + t150 * t219;
t152 = Icges(6,4) * t218 + Icges(6,6) * t216;
t102 = Icges(6,2) * t217 + t152 * t219;
t460 = t100 + t102;
t106 = Icges(6,4) * t217 + t219 * t268;
t158 = Icges(5,1) * t218 - t356;
t108 = Icges(5,5) * t217 + t158 * t219;
t452 = t106 + t108;
t148 = Icges(6,3) * t216 + t353;
t267 = -Icges(5,2) * t216 + t207;
t459 = t148 - t267;
t451 = (Icges(5,6) - Icges(6,6)) * t218 + t464 * t216;
t458 = t158 + t268;
t350 = Icges(5,6) * t219;
t103 = Icges(5,4) * t337 - Icges(5,2) * t339 - t350;
t345 = t103 * t216;
t97 = -Icges(6,6) * t219 + t148 * t217;
t428 = -t216 * t97 - t453 * t218 + t345;
t347 = Icges(5,3) * t219;
t99 = Icges(5,5) * t337 - Icges(5,6) * t339 - t347;
t457 = -t217 * t428 - t219 * t99;
t456 = t103 - t97;
t104 = Icges(5,6) * t217 + t219 * t267;
t336 = t218 * t219;
t176 = Icges(6,5) * t336;
t338 = t216 * t219;
t349 = Icges(6,6) * t217;
t98 = Icges(6,3) * t338 + t176 + t349;
t455 = t104 - t98;
t450 = t462 * qJD(4);
t449 = t461 * qJD(4);
t441 = -t462 * t216 + t461 * t218;
t448 = t460 * t217 + t452 * t336 + t98 * t338;
t101 = -Icges(6,2) * t219 + t152 * t217;
t92 = t217 * t101;
t447 = -t217 * t99 - t453 * t336 - t97 * t338 - t92;
t446 = t459 * qJD(4);
t445 = t458 * qJD(4);
t444 = -t150 - t152;
t393 = t451 * t219;
t394 = t451 * t217;
t335 = t219 * t101;
t421 = -t335 + t457;
t420 = -t103 * t338 - t447;
t419 = -t104 * t338 + t448;
t401 = t453 * t216 + t456 * t218;
t400 = t452 * t216 + t455 * t218;
t440 = t441 * t217 - t393;
t439 = t441 * t219 + t394;
t438 = t450 * t219 + (t217 * t267 - t350 - t97) * qJD(2);
t437 = t450 * t217 + (t148 * t219 - t104 + t349) * qJD(2);
t436 = -t449 * t219 + (-t158 * t217 + t463) * qJD(2);
t435 = -t452 * qJD(2) + t449 * t217;
t278 = t102 * t219 - t106 * t337 - t98 * t339;
t82 = t108 * t337;
t282 = t219 * t100 - t82;
t37 = -t104 * t339 - t282;
t434 = -t278 + t37;
t433 = t451 * qJD(4);
t344 = t104 * t216;
t432 = t216 * t98 + t452 * t218 - t344;
t431 = t217 * (-t219 * t462 + t452) - t219 * (-Icges(5,2) * t337 - t266 * t217 - t177 + t453);
t417 = rSges(6,3) + qJ(5);
t430 = t445 * t218 + t446 * t216 + (-t216 * t461 - t218 * t462) * qJD(4) + t451 * qJD(2);
t429 = t456 * t219 + (-Icges(6,1) * t338 + t155 * t219 + t176 - t455) * t217;
t427 = t460 * qJD(2);
t426 = t461 - t459;
t425 = -t462 + t458;
t422 = t441 * qJD(2) + qJD(4) * t444;
t376 = rSges(6,1) + pkin(4);
t418 = t439 * qJD(2);
t159 = pkin(4) * t216 - qJ(5) * t218;
t160 = rSges(6,1) * t216 - rSges(6,3) * t218;
t318 = t159 + t160;
t165 = rSges(6,1) * t218 + rSges(6,3) * t216;
t416 = pkin(4) * t218 + qJ(5) * t216 + t165;
t415 = -t400 * qJD(4) + t438 * t216 + t436 * t218 + t427;
t386 = qJD(2) * t101;
t414 = -qJD(2) * t99 + t401 * qJD(4) - t437 * t216 + t435 * t218 - t386;
t413 = (t419 * t217 - t420 * t219) * qJD(4);
t412 = (t434 * t217 - t421 * t219) * qJD(4);
t411 = t440 * qJD(2);
t410 = t428 * qJD(2) - t217 * t433 + t427;
t409 = -t386 - t433 * t219 + (-t150 * t217 + t347 - t432) * qJD(2);
t408 = 0.2e1 * qJD(4);
t407 = t411 + t412;
t406 = t413 + t418;
t405 = t428 * qJD(4) + t435 * t216 + t437 * t218;
t404 = t432 * qJD(4) + t436 * t216 - t438 * t218;
t403 = -t217 * t422 + t219 * t430;
t402 = t217 * t430 + t219 * t422;
t199 = t219 * qJ(3);
t162 = pkin(2) * t217 - t199;
t146 = qJD(2) * t162;
t224 = -pkin(6) - qJ(3);
t191 = t219 * t224;
t223 = cos(pkin(8));
t215 = pkin(3) * t223 + pkin(2);
t313 = -t217 * t215 - t191;
t95 = t162 + t313;
t399 = qJD(2) * t95 - t146;
t301 = qJD(5) * t218;
t210 = t217 * rSges(6,2);
t327 = t336 * t376 + t417 * t338 + t210;
t213 = t219 * rSges(6,2);
t328 = t217 * t416 - t213;
t31 = -t301 + qJD(1) + (t217 * t328 + t219 * t327) * qJD(4);
t398 = qJD(4) * t31;
t397 = t216 * t376;
t279 = t318 * qJD(4);
t302 = qJD(5) * t216;
t246 = -t279 + t302;
t396 = t246 * t217;
t208 = t217 * rSges(5,3);
t112 = rSges(5,1) * t336 - rSges(5,2) * t338 + t208;
t179 = t219 * t215;
t280 = -t217 * t224 + t179;
t395 = t112 + t280;
t198 = t217 * qJ(3);
t167 = t219 * pkin(2) + t198;
t366 = rSges(4,2) * sin(pkin(8));
t368 = rSges(4,1) * t223;
t259 = t217 * rSges(4,3) + (-t366 + t368) * t219;
t392 = t167 + t259;
t303 = qJD(4) * t219;
t291 = t218 * t303;
t305 = qJD(2) * t219;
t391 = rSges(6,2) * t305 + t291 * t417;
t390 = -t216 * t431 + t429 * t218;
t389 = (-t216 * t426 + t218 * t425) * qJD(2);
t388 = t335 + t448;
t387 = t444 * qJD(2);
t130 = t160 * t217;
t304 = qJD(4) * t217;
t371 = (pkin(4) * t305 + qJ(5) * t304) * t218 + (qJ(5) * t305 + (-pkin(4) * qJD(4) + qJD(5)) * t217) * t216 - qJD(4) * t130 + (t165 * t219 + t210) * qJD(2);
t172 = t219 * t302;
t306 = qJD(2) * t217;
t244 = -t216 * t303 - t218 * t306;
t294 = t216 * t306;
t372 = t244 * t376 - t417 * t294 + t172 + t391;
t1 = (t302 + t372 * t219 + t371 * t217 + (-t217 * t327 + t219 * t328) * qJD(2)) * qJD(4);
t379 = m(6) * t1;
t378 = t217 / 0.2e1;
t377 = -t219 / 0.2e1;
t374 = qJD(2) / 0.2e1;
t373 = pkin(2) - t215;
t367 = rSges(5,1) * t218;
t161 = rSges(5,1) * t216 + rSges(5,2) * t218;
t135 = t161 * t219;
t197 = qJD(3) * t219;
t293 = t161 * t304;
t43 = qJD(2) * t395 - t197 - t293;
t365 = t135 * t43;
t110 = rSges(5,1) * t337 - rSges(5,2) * t339 - t219 * rSges(5,3);
t196 = qJD(3) * t217;
t292 = t161 * t303;
t275 = t196 - t292;
t359 = t95 - t162;
t42 = (-t110 + t359) * qJD(2) + t275;
t364 = t217 * t42;
t128 = qJD(2) * t167 - t197;
t186 = t224 * t306;
t361 = -t128 + t186 - (-t219 * t373 - t198) * qJD(2);
t298 = qJD(2) * qJD(3);
t192 = qJ(3) * t305;
t309 = t192 + t196;
t326 = qJD(2) * (-pkin(2) * t306 + t309) + t217 * t298;
t325 = -qJD(4) * t416 + t301;
t324 = -t159 * t217 - t130;
t323 = t318 * t219;
t315 = rSges(5,2) * t294 + rSges(5,3) * t305;
t314 = t172 + t196;
t187 = t217 * t366;
t312 = rSges(4,3) * t305 + qJD(2) * t187;
t311 = t186 + t197;
t310 = t219 * rSges(4,3) + t187;
t297 = t217 * t368;
t295 = qJD(2) * (-t192 + (t217 * t373 - t191) * qJD(2)) + t326;
t290 = -pkin(2) - t368;
t287 = -t304 / 0.2e1;
t284 = t303 / 0.2e1;
t281 = -t99 + t344;
t273 = -rSges(5,2) * t216 + t367;
t272 = -t217 * t43 - t219 * t42;
t263 = t110 * t217 + t112 * t219;
t258 = qJD(4) * (t301 + t325);
t131 = t161 * t217;
t245 = -t219 * t279 + t314;
t239 = -t215 - t416;
t72 = rSges(5,1) * t244 - rSges(5,2) * t291 + t315;
t74 = -qJD(4) * t131 + (t219 * t273 + t208) * qJD(2);
t238 = t217 * t74 + t219 * t72 + (t110 * t219 - t112 * t217) * qJD(2);
t190 = t219 * t298;
t144 = t273 * qJD(4);
t113 = t297 - t310;
t78 = qJD(2) * t392 - t197;
t77 = t196 + (-t113 - t162) * qJD(2);
t56 = t190 + (-qJD(2) * t259 - t128) * qJD(2);
t55 = qJD(2) * (-qJD(2) * t297 + t312) + t326;
t44 = qJD(4) * t263 + qJD(1);
t30 = -t197 + t396 + (t327 + t280) * qJD(2);
t29 = (-t328 + t359) * qJD(2) + t245;
t26 = -t144 * t303 + t190 + (-t74 + t293 + t361) * qJD(2);
t25 = -t144 * t304 + (t72 - t292) * qJD(2) + t295;
t16 = t238 * qJD(4);
t11 = t190 + t219 * t258 + (t361 - t371 - t396) * qJD(2);
t10 = t217 * t258 + (t219 * t246 + t372) * qJD(2) + t295;
t2 = [m(5) * t16 + t379; (((t37 - t82 + (t100 + t345) * t219 + t447) * t219 + (t388 + t421 - t457) * t217) * qJD(4) + t418) * t284 + (t441 * qJD(4) + t445 * t216 - t446 * t218) * qJD(2) + (t11 * (t213 + t313) + t29 * t311 + t10 * (t179 + t327) + t30 * (t314 + t391) + (-t30 * qJD(4) * t397 + (-t30 * t224 + t239 * t29) * qJD(2)) * t219 + (-t10 * t224 - t11 * t376 * t218 + (-t29 * qJD(5) - t11 * t417) * t216 + t29 * (-t218 * t417 + t397) * qJD(4) + (-t29 * rSges(6,2) + t239 * t30) * qJD(2)) * t217 - (-qJD(2) * t328 + t245 - t29 + t399) * t30) * m(6) + (t26 * (-t110 + t313) + t42 * t311 + t25 * t395 + t43 * (t196 + t315) + (t161 * t364 - t365) * qJD(4) + ((-t42 * rSges(5,3) + t43 * (-t215 - t367)) * t217 + (t42 * (-t215 - t273) - t43 * t224) * t219) * qJD(2) - (-qJD(2) * t110 + t275 + t399 - t42) * t43) * m(5) + (t56 * (t217 * t290 + t199 + t310) + t77 * t197 + t55 * t392 + t78 * (t309 + t312) + (t77 * (t290 + t366) * t219 + (t77 * (-rSges(4,3) - qJ(3)) + t78 * t290) * t217) * qJD(2) - (-qJD(2) * t113 - t146 + t196 - t77) * t78) * m(4) + (((t219 * t281 - t388 + t419) * t219 + (t217 * t281 + t278 + t282 + t420 - t92) * t217) * qJD(4) + t407 - t411) * t287 + (t403 + t404) * t304 / 0.2e1 - (t402 - t405 + t406) * t303 / 0.2e1 + ((t401 + t440) * t217 + (t400 + t439) * t219) * qJD(4) * t374; 0.2e1 * (t10 * t377 + t11 * t378) * m(6) + 0.2e1 * (t25 * t377 + t26 * t378) * m(5) + 0.2e1 * (t377 * t55 + t378 * t56) * m(4); -((t429 * t216 + t218 * t431) * qJD(4) + (t425 * t216 + t426 * t218) * qJD(2)) * qJD(2) / 0.2e1 + (t405 * t219 + t404 * t217 + (t217 * t401 + t219 * t400) * qJD(2)) * t374 + ((-t304 * t393 - t387) * t217 + ((t217 * t394 + t390) * qJD(4) + t389) * t219) * t287 + ((-t303 * t394 + t387) * t219 + ((t219 * t393 + t390) * qJD(4) + t389) * t217) * t284 + ((-t11 * t318 + t29 * t325 + t1 * t327 + t31 * t372 + (-t30 * t318 + t31 * t328) * qJD(2)) * t219 + (-t10 * t318 + t30 * t325 + t1 * t328 + t31 * t371 + (t29 * t318 - t31 * t327) * qJD(2)) * t217 - (t216 * t31 + (t217 * t30 + t219 * t29) * t218) * qJD(5) - (-t29 * t324 - t30 * t323) * qJD(2) - ((-t29 * t416 - t31 * t323) * t219 + (-t30 * t416 + t31 * t324) * t217) * qJD(4)) * m(6) + (-(t131 * t42 - t365) * qJD(2) - (t44 * (-t131 * t217 - t135 * t219) + t272 * t273) * qJD(4) + t16 * t263 + t44 * t238 + t272 * t144 + (-t25 * t217 - t26 * t219 + (-t219 * t43 + t364) * qJD(2)) * t161) * m(5) + (t403 * qJD(2) + ((t419 * qJD(2) + t414 * t219) * t219 + (t409 * t217 + t420 * qJD(2) + (-t410 + t415) * t219) * t217) * t408) * t378 + (t402 * qJD(2) + ((t434 * qJD(2) + t410 * t219) * t219 + (t415 * t217 + t421 * qJD(2) + (-t409 + t414) * t219) * t217) * t408) * t377 + (t407 + t412) * t306 / 0.2e1 + (t406 + t413) * t305 / 0.2e1; -t218 * t379 + 0.2e1 * (m(6) * (t10 * t217 + t11 * t219 + t398) / 0.2e1 - m(6) * (t217 ^ 2 + t219 ^ 2) * t398 / 0.2e1) * t216;];
tauc = t2(:);
