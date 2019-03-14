% Calculate time derivative of joint inertia matrix for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:49
% EndTime: 2019-03-08 20:27:24
% DurationCPUTime: 23.09s
% Computational Cost: add. (131389->1153), mult. (331141->1594), div. (0->0), fcn. (396847->14), ass. (0->495)
t565 = sin(pkin(12));
t566 = cos(pkin(12));
t576 = sin(qJ(2));
t578 = cos(qJ(2));
t457 = t565 * t576 - t566 * t578;
t428 = t457 * qJD(2);
t445 = sin(pkin(11));
t447 = cos(pkin(11));
t567 = cos(pkin(6));
t480 = t567 * t565;
t481 = t567 * t566;
t414 = t480 * t578 + t481 * t576;
t452 = qJD(2) * t414;
t370 = t445 * t428 - t447 * t452;
t453 = -t480 * t576 + t481 * t578;
t408 = t453 * qJD(2);
t458 = t565 * t578 + t566 * t576;
t429 = t458 * qJD(2);
t371 = t408 * t447 - t429 * t445;
t486 = t567 * t578;
t424 = -t445 * t576 + t447 * t486;
t415 = t424 * qJD(2);
t485 = t567 * t576;
t463 = -t445 * t578 - t447 * t485;
t416 = t463 * qJD(2);
t607 = -Icges(3,5) * t415 - Icges(4,5) * t371 - Icges(3,6) * t416 - Icges(4,6) * t370;
t446 = sin(pkin(6));
t413 = t458 * t446;
t406 = qJD(2) * t413;
t533 = qJD(2) * t446;
t407 = t457 * t533;
t606 = -Icges(4,5) * t407 - Icges(4,6) * t406 + (Icges(3,5) * t578 - Icges(3,6) * t576) * t533;
t372 = t428 * t447 + t445 * t452;
t373 = -t408 * t445 - t429 * t447;
t426 = -t445 * t486 - t447 * t576;
t417 = t426 * qJD(2);
t464 = t445 * t485 - t447 * t578;
t418 = t464 * qJD(2);
t605 = Icges(3,5) * t417 + Icges(4,5) * t373 + Icges(3,6) * t418 + Icges(4,6) * t372;
t444 = qJ(5) + qJ(6);
t441 = sin(t444);
t442 = cos(t444);
t394 = -t414 * t445 - t447 * t457;
t449 = sin(qJ(4));
t556 = t446 * t449;
t577 = cos(qJ(4));
t358 = t394 * t577 + t445 * t556;
t443 = qJD(5) + qJD(6);
t495 = -t358 * t443 - t372;
t516 = t446 * t577;
t470 = -t394 * t449 + t445 * t516;
t308 = qJD(4) * t470 + t373 * t577;
t393 = -t445 * t453 - t447 * t458;
t498 = t393 * t443 - t308;
t218 = t441 * t498 + t442 * t495;
t219 = t441 * t495 - t442 * t498;
t307 = qJD(4) * t358 + t373 * t449;
t164 = rSges(7,1) * t219 + rSges(7,2) * t218 + rSges(7,3) * t307;
t311 = -t358 * t441 - t393 * t442;
t312 = t358 * t442 - t393 * t441;
t230 = rSges(7,1) * t312 + rSges(7,2) * t311 - rSges(7,3) * t470;
t392 = t414 * t447 - t445 * t457;
t356 = t392 * t577 - t447 * t556;
t305 = qJD(4) * t356 + t371 * t449;
t469 = -t392 * t449 - t447 * t516;
t496 = -t356 * t443 - t370;
t306 = qJD(4) * t469 + t371 * t577;
t391 = -t445 * t458 + t447 * t453;
t499 = t391 * t443 - t306;
t216 = t441 * t499 + t442 * t496;
t217 = t441 * t496 - t442 * t499;
t163 = rSges(7,1) * t217 + rSges(7,2) * t216 + rSges(7,3) * t305;
t309 = -t356 * t441 - t391 * t442;
t310 = t356 * t442 - t391 * t441;
t229 = rSges(7,1) * t310 + rSges(7,2) * t309 - rSges(7,3) * t469;
t555 = -t163 * t470 + t229 * t307;
t96 = t164 * t469 - t230 * t305 + t555;
t231 = Icges(5,5) * t306 - Icges(5,6) * t305 - Icges(5,3) * t370;
t233 = Icges(5,4) * t306 - Icges(5,2) * t305 - Icges(5,6) * t370;
t235 = Icges(5,1) * t306 - Icges(5,4) * t305 - Icges(5,5) * t370;
t277 = Icges(5,5) * t356 + Icges(5,6) * t469 - Icges(5,3) * t391;
t279 = Icges(5,4) * t356 + Icges(5,2) * t469 - Icges(5,6) * t391;
t281 = Icges(5,1) * t356 + Icges(5,4) * t469 - Icges(5,5) * t391;
t107 = -t231 * t393 + t233 * t470 + t235 * t358 - t277 * t372 - t279 * t307 + t281 * t308;
t232 = Icges(5,5) * t308 - Icges(5,6) * t307 - Icges(5,3) * t372;
t234 = Icges(5,4) * t308 - Icges(5,2) * t307 - Icges(5,6) * t372;
t236 = Icges(5,1) * t308 - Icges(5,4) * t307 - Icges(5,5) * t372;
t278 = Icges(5,5) * t358 + Icges(5,6) * t470 - Icges(5,3) * t393;
t280 = Icges(5,4) * t358 + Icges(5,2) * t470 - Icges(5,6) * t393;
t282 = Icges(5,1) * t358 + Icges(5,4) * t470 - Icges(5,5) * t393;
t108 = -t232 * t393 + t234 * t470 + t236 * t358 - t278 * t372 - t280 * t307 + t282 * t308;
t401 = t413 * t577 + t449 * t567;
t351 = qJD(4) * t401 - t407 * t449;
t467 = -t413 * t449 + t567 * t577;
t352 = qJD(4) * t467 - t407 * t577;
t292 = Icges(5,5) * t352 - Icges(5,6) * t351 + Icges(5,3) * t406;
t293 = Icges(5,4) * t352 - Icges(5,2) * t351 + Icges(5,6) * t406;
t294 = Icges(5,1) * t352 - Icges(5,4) * t351 + Icges(5,5) * t406;
t412 = t457 * t446;
t344 = Icges(5,5) * t401 + Icges(5,6) * t467 + Icges(5,3) * t412;
t345 = Icges(5,4) * t401 + Icges(5,2) * t467 + Icges(5,6) * t412;
t346 = Icges(5,1) * t401 + Icges(5,4) * t467 + Icges(5,5) * t412;
t123 = -t292 * t393 + t293 * t470 + t294 * t358 - t307 * t345 + t308 * t346 - t344 * t372;
t157 = Icges(7,5) * t217 + Icges(7,6) * t216 + Icges(7,3) * t305;
t159 = Icges(7,4) * t217 + Icges(7,2) * t216 + Icges(7,6) * t305;
t161 = Icges(7,1) * t217 + Icges(7,4) * t216 + Icges(7,5) * t305;
t223 = Icges(7,5) * t310 + Icges(7,6) * t309 - Icges(7,3) * t469;
t225 = Icges(7,4) * t310 + Icges(7,2) * t309 - Icges(7,6) * t469;
t227 = Icges(7,1) * t310 + Icges(7,4) * t309 - Icges(7,5) * t469;
t44 = -t157 * t470 + t159 * t311 + t161 * t312 + t218 * t225 + t219 * t227 + t223 * t307;
t158 = Icges(7,5) * t219 + Icges(7,6) * t218 + Icges(7,3) * t307;
t160 = Icges(7,4) * t219 + Icges(7,2) * t218 + Icges(7,6) * t307;
t162 = Icges(7,1) * t219 + Icges(7,4) * t218 + Icges(7,5) * t307;
t224 = Icges(7,5) * t312 + Icges(7,6) * t311 - Icges(7,3) * t470;
t226 = Icges(7,4) * t312 + Icges(7,2) * t311 - Icges(7,6) * t470;
t228 = Icges(7,1) * t312 + Icges(7,4) * t311 - Icges(7,5) * t470;
t45 = -t158 * t470 + t160 * t311 + t162 * t312 + t218 * t226 + t219 * t228 + t224 * t307;
t494 = -t401 * t443 + t406;
t497 = t412 * t443 + t352;
t270 = -t441 * t497 + t442 * t494;
t271 = t441 * t494 + t442 * t497;
t198 = Icges(7,5) * t271 + Icges(7,6) * t270 + Icges(7,3) * t351;
t199 = Icges(7,4) * t271 + Icges(7,2) * t270 + Icges(7,6) * t351;
t200 = Icges(7,1) * t271 + Icges(7,4) * t270 + Icges(7,5) * t351;
t349 = -t401 * t441 + t412 * t442;
t350 = t401 * t442 + t412 * t441;
t272 = Icges(7,5) * t350 + Icges(7,6) * t349 - Icges(7,3) * t467;
t273 = Icges(7,4) * t350 + Icges(7,2) * t349 - Icges(7,6) * t467;
t274 = Icges(7,1) * t350 + Icges(7,4) * t349 - Icges(7,5) * t467;
t92 = -t198 * t470 + t199 * t311 + t200 * t312 + t218 * t273 + t219 * t274 + t272 * t307;
t28 = t92 * t567 + (-t44 * t447 + t445 * t45) * t446;
t450 = cos(qJ(5));
t448 = sin(qJ(5));
t561 = t391 * t448;
t319 = t356 * t450 - t561;
t245 = -qJD(5) * t319 - t306 * t448 - t370 * t450;
t318 = -t356 * t448 - t391 * t450;
t461 = qJD(5) * t318 - t370 * t448;
t246 = t306 * t450 + t461;
t166 = Icges(6,5) * t246 + Icges(6,6) * t245 + Icges(6,3) * t305;
t168 = Icges(6,4) * t246 + Icges(6,2) * t245 + Icges(6,6) * t305;
t170 = Icges(6,1) * t246 + Icges(6,4) * t245 + Icges(6,5) * t305;
t239 = Icges(6,5) * t319 + Icges(6,6) * t318 - Icges(6,3) * t469;
t241 = Icges(6,4) * t319 + Icges(6,2) * t318 - Icges(6,6) * t469;
t243 = Icges(6,1) * t319 + Icges(6,4) * t318 - Icges(6,5) * t469;
t560 = t393 * t448;
t321 = t358 * t450 - t560;
t247 = -qJD(5) * t321 - t308 * t448 - t372 * t450;
t320 = -t358 * t448 - t393 * t450;
t460 = qJD(5) * t320 - t372 * t448;
t248 = t308 * t450 + t460;
t50 = -t166 * t470 + t168 * t320 + t170 * t321 + t239 * t307 + t241 * t247 + t243 * t248;
t167 = Icges(6,5) * t248 + Icges(6,6) * t247 + Icges(6,3) * t307;
t169 = Icges(6,4) * t248 + Icges(6,2) * t247 + Icges(6,6) * t307;
t171 = Icges(6,1) * t248 + Icges(6,4) * t247 + Icges(6,5) * t307;
t240 = Icges(6,5) * t321 + Icges(6,6) * t320 - Icges(6,3) * t470;
t242 = Icges(6,4) * t321 + Icges(6,2) * t320 - Icges(6,6) * t470;
t244 = Icges(6,1) * t321 + Icges(6,4) * t320 - Icges(6,5) * t470;
t51 = -t167 * t470 + t169 * t320 + t171 * t321 + t240 * t307 + t242 * t247 + t244 * t248;
t559 = t412 * t448;
t354 = t401 * t450 + t559;
t290 = -qJD(5) * t354 - t352 * t448 + t406 * t450;
t353 = -t401 * t448 + t412 * t450;
t459 = qJD(5) * t353 + t406 * t448;
t291 = t352 * t450 + t459;
t205 = Icges(6,5) * t291 + Icges(6,6) * t290 + Icges(6,3) * t351;
t206 = Icges(6,4) * t291 + Icges(6,2) * t290 + Icges(6,6) * t351;
t207 = Icges(6,1) * t291 + Icges(6,4) * t290 + Icges(6,5) * t351;
t286 = Icges(6,5) * t354 + Icges(6,6) * t353 - Icges(6,3) * t467;
t287 = Icges(6,4) * t354 + Icges(6,2) * t353 - Icges(6,6) * t467;
t288 = Icges(6,1) * t354 + Icges(6,4) * t353 - Icges(6,5) * t467;
t95 = -t205 * t470 + t206 * t320 + t207 * t321 + t247 * t287 + t248 * t288 + t286 * t307;
t30 = t95 * t567 + (t445 * t51 - t447 * t50) * t446;
t604 = t28 + t30 + t123 * t567 + (-t107 * t447 + t108 * t445) * t446;
t115 = t231 * t412 + t233 * t467 + t235 * t401 + t277 * t406 - t279 * t351 + t281 * t352;
t116 = t232 * t412 + t234 * t467 + t236 * t401 + t278 * t406 - t280 * t351 + t282 * t352;
t141 = t292 * t412 + t293 * t467 + t294 * t401 + t344 * t406 - t345 * t351 + t346 * t352;
t52 = -t157 * t467 + t159 * t349 + t161 * t350 + t223 * t351 + t225 * t270 + t227 * t271;
t53 = -t158 * t467 + t160 * t349 + t162 * t350 + t224 * t351 + t226 * t270 + t228 * t271;
t99 = -t198 * t467 + t199 * t349 + t200 * t350 + t270 * t273 + t271 * t274 + t272 * t351;
t32 = t99 * t567 + (t445 * t53 - t447 * t52) * t446;
t102 = -t205 * t467 + t206 * t353 + t207 * t354 + t286 * t351 + t287 * t290 + t288 * t291;
t57 = -t166 * t467 + t168 * t353 + t170 * t354 + t239 * t351 + t241 * t290 + t243 * t291;
t58 = -t167 * t467 + t169 * t353 + t171 * t354 + t240 * t351 + t242 * t290 + t244 * t291;
t33 = t102 * t567 + (t445 * t58 - t447 * t57) * t446;
t603 = t32 + t33 + t141 * t567 + (-t115 * t447 + t116 * t445) * t446;
t105 = -t231 * t391 + t233 * t469 + t235 * t356 - t277 * t370 - t279 * t305 + t281 * t306;
t106 = -t232 * t391 + t234 * t469 + t236 * t356 - t278 * t370 - t280 * t305 + t282 * t306;
t122 = -t292 * t391 + t293 * t469 + t294 * t356 - t305 * t345 + t306 * t346 - t344 * t370;
t42 = -t157 * t469 + t159 * t309 + t161 * t310 + t216 * t225 + t217 * t227 + t223 * t305;
t43 = -t158 * t469 + t160 * t309 + t162 * t310 + t216 * t226 + t217 * t228 + t224 * t305;
t91 = -t198 * t469 + t199 * t309 + t200 * t310 + t216 * t273 + t217 * t274 + t272 * t305;
t27 = t91 * t567 + (-t42 * t447 + t43 * t445) * t446;
t48 = -t166 * t469 + t168 * t318 + t170 * t319 + t239 * t305 + t241 * t245 + t243 * t246;
t49 = -t167 * t469 + t169 * t318 + t171 * t319 + t240 * t305 + t242 * t245 + t244 * t246;
t94 = -t205 * t469 + t206 * t318 + t207 * t319 + t245 * t287 + t246 * t288 + t286 * t305;
t29 = t94 * t567 + (t445 * t49 - t447 * t48) * t446;
t602 = t122 * t567 + (-t105 * t447 + t106 * t445) * t446 + t29 + t27;
t601 = t446 ^ 2;
t569 = pkin(5) * t450;
t183 = pkin(5) * t461 + pkin(10) * t305 + t306 * t569;
t553 = t163 + t183;
t184 = pkin(5) * t460 + pkin(10) * t307 + t308 * t569;
t552 = t164 + t184;
t221 = -pkin(5) * t561 - pkin(10) * t469 + t356 * t569;
t549 = t221 + t229;
t222 = -pkin(5) * t560 - pkin(10) * t470 + t358 * t569;
t548 = t222 + t230;
t261 = pkin(4) * t306 + pkin(9) * t305;
t302 = pkin(4) * t356 - pkin(9) * t469;
t545 = -t261 * t393 - t302 * t372;
t431 = pkin(2) * qJD(2) * t486 - qJD(3) * t446;
t510 = qJD(2) * t576;
t493 = pkin(2) * t510;
t402 = t431 * t447 - t445 * t493;
t403 = -t431 * t445 - t447 * t493;
t557 = t446 * t447;
t558 = t445 * t446;
t534 = t402 * t558 + t403 * t557;
t591 = m(7) / 0.2e1;
t592 = m(6) / 0.2e1;
t599 = t592 + t591;
t124 = -t223 * t469 + t225 * t309 + t227 * t310;
t125 = -t224 * t469 + t226 * t309 + t228 * t310;
t147 = -t272 * t469 + t273 * t309 + t274 * t310;
t11 = -t124 * t370 - t125 * t372 + t147 * t406 - t391 * t42 - t393 * t43 + t412 * t91;
t129 = -t239 * t469 + t241 * t318 + t243 * t319;
t130 = -t240 * t469 + t242 * t318 + t244 * t319;
t154 = -t286 * t469 + t287 * t318 + t288 * t319;
t15 = -t129 * t370 - t130 * t372 + t154 * t406 - t391 * t48 - t393 * t49 + t412 * t94;
t178 = -t277 * t391 + t279 * t469 + t281 * t356;
t179 = -t278 * t391 + t280 * t469 + t282 * t356;
t195 = -t344 * t391 + t345 * t469 + t346 * t356;
t598 = -t105 * t391 - t106 * t393 + t122 * t412 - t178 * t370 - t179 * t372 + t195 * t406 + t11 + t15;
t126 = -t223 * t470 + t225 * t311 + t227 * t312;
t127 = -t224 * t470 + t226 * t311 + t228 * t312;
t148 = -t272 * t470 + t273 * t311 + t274 * t312;
t12 = -t126 * t370 - t127 * t372 + t148 * t406 - t391 * t44 - t393 * t45 + t412 * t92;
t131 = -t239 * t470 + t241 * t320 + t243 * t321;
t132 = -t240 * t470 + t242 * t320 + t244 * t321;
t155 = -t286 * t470 + t287 * t320 + t288 * t321;
t16 = -t131 * t370 - t132 * t372 + t155 * t406 - t391 * t50 - t393 * t51 + t412 * t95;
t180 = -t277 * t393 + t279 * t470 + t281 * t358;
t181 = -t278 * t393 + t280 * t470 + t282 * t358;
t196 = -t344 * t393 + t345 * t470 + t346 * t358;
t597 = -t107 * t391 - t108 * t393 + t123 * t412 - t180 * t370 - t181 * t372 + t196 * t406 + t12 + t16;
t186 = t277 * t412 + t279 * t467 + t281 * t401;
t187 = t278 * t412 + t280 * t467 + t282 * t401;
t136 = -t223 * t467 + t225 * t349 + t227 * t350;
t137 = -t224 * t467 + t226 * t349 + t228 * t350;
t177 = -t272 * t467 + t273 * t349 + t274 * t350;
t22 = -t136 * t370 - t137 * t372 + t177 * t406 - t391 * t52 - t393 * t53 + t412 * t99;
t138 = -t239 * t467 + t241 * t353 + t243 * t354;
t139 = -t240 * t467 + t242 * t353 + t244 * t354;
t185 = -t286 * t467 + t287 * t353 + t288 * t354;
t24 = t102 * t412 - t138 * t370 - t139 * t372 + t185 * t406 - t391 * t57 - t393 * t58;
t209 = t344 * t412 + t345 * t467 + t346 * t401;
t564 = t209 * t406;
t596 = -t115 * t391 - t116 * t393 + t141 * t412 - t186 * t370 - t187 * t372 + t22 + t24 + t564;
t237 = rSges(5,1) * t306 - rSges(5,2) * t305 - rSges(5,3) * t370;
t172 = rSges(6,1) * t246 + rSges(6,2) * t245 + rSges(6,3) * t305;
t575 = m(6) * t172;
t595 = -m(5) * t237 - m(7) * t553 - t575;
t594 = m(4) / 0.2e1;
t593 = m(5) / 0.2e1;
t590 = t305 / 0.2e1;
t589 = t307 / 0.2e1;
t588 = t351 / 0.2e1;
t587 = -t469 / 0.2e1;
t586 = -t470 / 0.2e1;
t585 = -t370 / 0.2e1;
t584 = -t372 / 0.2e1;
t583 = -t391 / 0.2e1;
t582 = -t393 / 0.2e1;
t581 = -t467 / 0.2e1;
t580 = t406 / 0.2e1;
t579 = t412 / 0.2e1;
t173 = rSges(6,1) * t248 + rSges(6,2) * t247 + rSges(6,3) * t307;
t574 = m(6) * t173;
t249 = rSges(6,1) * t319 + rSges(6,2) * t318 - rSges(6,3) * t469;
t573 = m(6) * t249;
t250 = rSges(6,1) * t321 + rSges(6,2) * t320 - rSges(6,3) * t470;
t572 = m(6) * t250;
t571 = pkin(2) * t446;
t570 = pkin(2) * t578;
t554 = -t164 * t467 + t230 * t351;
t203 = rSges(7,1) * t271 + rSges(7,2) * t270 + rSges(7,3) * t351;
t275 = rSges(7,1) * t350 + rSges(7,2) * t349 - rSges(7,3) * t467;
t551 = -t203 * t469 + t275 * t305;
t210 = pkin(5) * t459 + pkin(10) * t351 + t352 * t569;
t550 = -t203 - t210;
t547 = -t249 - t302;
t303 = pkin(4) * t358 - pkin(9) * t470;
t546 = t250 + t303;
t262 = t308 * pkin(4) + t307 * pkin(9);
t544 = t262 * t412 + t303 * t406;
t298 = pkin(4) * t352 + t351 * pkin(9);
t348 = pkin(4) * t401 - pkin(9) * t467;
t543 = -t298 * t391 - t348 * t370;
t276 = pkin(5) * t559 - pkin(10) * t467 + t401 * t569;
t542 = -t275 - t276;
t289 = rSges(6,1) * t354 + rSges(6,2) * t353 - rSges(6,3) * t467;
t541 = t289 + t348;
t332 = pkin(3) * t373 - pkin(8) * t372;
t399 = t567 * t403;
t540 = t332 * t567 + t399;
t343 = pkin(3) * t394 - pkin(8) * t393;
t509 = pkin(2) * t485 - qJ(3) * t446;
t381 = -t445 * t509 + t447 * t570;
t369 = t567 * t381;
t539 = t343 * t567 + t369;
t380 = t445 * t570 + t447 * t509;
t538 = t380 * t558 + t381 * t557;
t511 = qJD(2) * t578;
t430 = qJD(3) * t567 + t511 * t571;
t537 = pkin(3) * t407 - pkin(8) * t406 - t430;
t434 = qJ(3) * t567 + t571 * t576;
t536 = -pkin(3) * t413 - pkin(8) * t412 - t434;
t535 = 0.2e1 * t534;
t530 = 0.2e1 * t567;
t529 = 0.2e1 * t446;
t527 = t567 / 0.2e1;
t526 = -t302 - t549;
t525 = t303 + t548;
t524 = t262 * t567 + t540;
t523 = t348 - t542;
t522 = -t298 + t537;
t521 = t303 * t567 + t539;
t520 = -t348 + t536;
t519 = m(6) * t567;
t518 = m(7) * t567;
t515 = t578 * Icges(3,4);
t514 = t576 * Icges(3,4);
t513 = t558 / 0.2e1;
t512 = -t557 / 0.2e1;
t508 = t567 * t380;
t507 = t567 * t402;
t506 = 0.2e1 * m(5);
t504 = 0.2e1 * m(6);
t502 = 0.2e1 * m(7);
t501 = t446 * (-rSges(4,1) * t413 + rSges(4,2) * t412 - rSges(4,3) * t567 - t434);
t500 = (rSges(4,1) * t407 + rSges(4,2) * t406 - t430) * t446;
t331 = pkin(3) * t371 - pkin(8) * t370;
t316 = t331 * t558;
t317 = t332 * t557;
t491 = 0.2e1 * t316 + 0.2e1 * t317 + t535;
t490 = t316 + t317 + t534;
t342 = pkin(3) * t392 - pkin(8) * t391;
t489 = t342 * t558 + t343 * t557 + t538;
t295 = rSges(5,1) * t352 - rSges(5,2) * t351 + rSges(5,3) * t406;
t488 = (-t295 + t537) * t446;
t347 = rSges(5,1) * t401 + rSges(5,2) * t467 + rSges(5,3) * t412;
t487 = (-t347 + t536) * t446;
t20 = t136 * t305 + t137 * t307 + t177 * t351 - t467 * t99 - t469 * t52 - t470 * t53;
t65 = -t124 * t469 - t125 * t470 - t147 * t467;
t66 = -t126 * t469 - t127 * t470 - t148 * t467;
t7 = t124 * t305 + t125 * t307 + t147 * t351 - t42 * t469 - t43 * t470 - t467 * t91;
t8 = t126 * t305 + t127 * t307 + t148 * t351 - t44 * t469 - t45 * t470 - t467 * t92;
t83 = -t136 * t469 - t137 * t470 - t177 * t467;
t479 = -t20 * t467 + t305 * t65 + t307 * t66 + t351 * t83 - t469 * t7 - t470 * t8;
t208 = rSges(6,1) * t291 + rSges(6,2) * t290 + rSges(6,3) * t351;
t478 = (-t208 + t522) * t446;
t477 = (-t289 + t520) * t446;
t258 = t261 * t558;
t259 = t262 * t557;
t475 = t258 + t259 + t490;
t474 = t302 * t558 + t303 * t557 + t489;
t473 = 0.2e1 * t96 * t591;
t472 = (t522 + t550) * t446;
t471 = (t520 + t542) * t446;
t73 = t147 * t567 + (-t124 * t447 + t125 * t445) * t446;
t74 = t148 * t567 + (-t126 * t447 + t127 * t445) * t446;
t87 = t177 * t567 + (-t136 * t447 + t137 * t445) * t446;
t468 = t20 * t527 + t27 * t587 + t28 * t586 + t32 * t581 + t512 * t7 + t513 * t8 + t588 * t87 + t589 * t74 + t590 * t73;
t466 = -t331 * t567 - t507;
t465 = -t342 * t567 - t508;
t69 = -t124 * t391 - t125 * t393 + t147 * t412;
t70 = -t126 * t391 - t127 * t393 + t148 * t412;
t84 = -t136 * t391 - t137 * t393 + t177 * t412;
t462 = t11 * t587 + t12 * t586 + t20 * t579 + t22 * t581 + t580 * t83 + t582 * t8 + t583 * t7 + t584 * t66 + t585 * t65 + t588 * t84 + t589 * t70 + t590 * t69;
t238 = rSges(5,1) * t308 - rSges(5,2) * t307 - rSges(5,3) * t372;
t456 = m(5) * t238 + m(7) * t552 + t574;
t455 = -t261 * t567 + t466;
t454 = -t302 * t567 + t465;
t422 = (rSges(3,1) * t578 - rSges(3,2) * t576) * t533;
t421 = (Icges(3,1) * t578 - t514) * t533;
t420 = (-Icges(3,2) * t576 + t515) * t533;
t411 = t567 * rSges(3,3) + (rSges(3,1) * t576 + rSges(3,2) * t578) * t446;
t410 = Icges(3,5) * t567 + (Icges(3,1) * t576 + t515) * t446;
t409 = Icges(3,6) * t567 + (Icges(3,2) * t578 + t514) * t446;
t390 = rSges(3,1) * t417 + rSges(3,2) * t418;
t388 = rSges(3,1) * t415 + rSges(3,2) * t416;
t387 = Icges(3,1) * t417 + Icges(3,4) * t418;
t386 = Icges(3,1) * t415 + Icges(3,4) * t416;
t385 = Icges(3,4) * t417 + Icges(3,2) * t418;
t384 = Icges(3,4) * t415 + Icges(3,2) * t416;
t379 = -rSges(3,1) * t464 + rSges(3,2) * t426 + rSges(3,3) * t558;
t378 = -rSges(3,1) * t463 + rSges(3,2) * t424 - rSges(3,3) * t557;
t377 = -Icges(3,1) * t464 + Icges(3,4) * t426 + Icges(3,5) * t558;
t376 = -Icges(3,1) * t463 + Icges(3,4) * t424 - Icges(3,5) * t557;
t375 = -Icges(3,4) * t464 + Icges(3,2) * t426 + Icges(3,6) * t558;
t374 = -Icges(3,4) * t463 + Icges(3,2) * t424 - Icges(3,6) * t557;
t366 = Icges(4,1) * t413 - Icges(4,4) * t412 + Icges(4,5) * t567;
t365 = Icges(4,4) * t413 - Icges(4,2) * t412 + Icges(4,6) * t567;
t363 = -Icges(4,1) * t407 - Icges(4,4) * t406;
t362 = -Icges(4,4) * t407 - Icges(4,2) * t406;
t338 = rSges(4,1) * t394 + rSges(4,2) * t393 + rSges(4,3) * t558;
t337 = rSges(4,1) * t392 + rSges(4,2) * t391 - rSges(4,3) * t557;
t336 = Icges(4,1) * t394 + Icges(4,4) * t393 + Icges(4,5) * t558;
t335 = Icges(4,1) * t392 + Icges(4,4) * t391 - Icges(4,5) * t557;
t334 = Icges(4,4) * t394 + Icges(4,2) * t393 + Icges(4,6) * t558;
t333 = Icges(4,4) * t392 + Icges(4,2) * t391 - Icges(4,6) * t557;
t330 = rSges(4,1) * t373 + rSges(4,2) * t372;
t329 = rSges(4,1) * t371 + rSges(4,2) * t370;
t328 = Icges(4,1) * t373 + Icges(4,4) * t372;
t327 = Icges(4,1) * t371 + Icges(4,4) * t370;
t326 = Icges(4,4) * t373 + Icges(4,2) * t372;
t325 = Icges(4,4) * t371 + Icges(4,2) * t370;
t313 = t391 * t348;
t297 = t412 * t303;
t285 = t393 * t302;
t284 = rSges(5,1) * t358 + rSges(5,2) * t470 - rSges(5,3) * t393;
t283 = rSges(5,1) * t356 + rSges(5,2) * t469 - rSges(5,3) * t391;
t265 = -t329 * t567 + t447 * t500 - t507;
t264 = t330 * t567 + t445 * t500 + t399;
t263 = t469 * t275;
t220 = (t329 * t445 + t330 * t447) * t446 + t534;
t215 = t284 * t412 + t347 * t393;
t214 = -t283 * t412 - t347 * t391;
t213 = t467 * t230;
t212 = t470 * t229;
t197 = -t283 * t393 + t284 * t391;
t193 = -t283 * t567 + t447 * t487 + t465;
t192 = t284 * t567 + t445 * t487 + t539;
t191 = -t250 * t467 + t289 * t470;
t190 = t249 * t467 - t289 * t469;
t189 = t275 * t470 - t213;
t188 = t229 * t467 - t263;
t182 = (t283 * t445 + t284 * t447) * t446 + t489;
t176 = -t237 * t567 + t447 * t488 + t466;
t175 = t238 * t567 + t445 * t488 + t540;
t174 = -t249 * t470 + t250 * t469;
t165 = t230 * t469 - t212;
t153 = t250 * t412 + t393 * t541 + t297;
t152 = -t289 * t391 + t412 * t547 - t313;
t146 = t238 * t412 + t284 * t406 + t295 * t393 + t347 * t372;
t145 = -t237 * t412 - t283 * t406 - t295 * t391 - t347 * t370;
t144 = -t249 * t567 + t447 * t477 + t454;
t143 = t250 * t567 + t445 * t477 + t521;
t142 = (t237 * t445 + t238 * t447) * t446 + t490;
t140 = -t249 * t393 + t391 * t546 - t285;
t135 = -t237 * t393 + t238 * t391 - t283 * t372 + t284 * t370;
t134 = -t222 * t467 - t470 * t542 - t213;
t133 = -t276 * t469 + t467 * t549 - t263;
t128 = (t249 * t445 + t250 * t447) * t446 + t474;
t121 = t393 * t523 + t412 * t548 + t297;
t120 = t391 * t542 + t412 * t526 - t313;
t119 = t447 * t471 - t549 * t567 + t454;
t118 = t445 * t471 + t548 * t567 + t521;
t117 = -t221 * t470 + t469 * t548 - t212;
t114 = -t172 * t567 + t447 * t478 + t455;
t113 = t173 * t567 + t445 * t478 + t524;
t112 = t391 * t525 - t393 * t549 - t285;
t111 = (t445 * t549 + t447 * t548) * t446 + t474;
t110 = -t173 * t467 + t208 * t470 + t250 * t351 - t289 * t307;
t109 = t172 * t467 - t208 * t469 - t249 * t351 + t289 * t305;
t104 = t203 * t470 - t275 * t307 + t554;
t103 = t163 * t467 - t229 * t351 + t551;
t101 = (t172 * t445 + t173 * t447) * t446 + t475;
t100 = -t172 * t470 + t173 * t469 + t249 * t307 - t250 * t305;
t98 = t173 * t412 + t250 * t406 + (t208 + t298) * t393 + t541 * t372 + t544;
t97 = -t208 * t391 - t289 * t370 + (-t172 - t261) * t412 + t547 * t406 + t543;
t93 = t185 * t567 + (-t138 * t447 + t139 * t445) * t446;
t90 = t447 * t472 - t553 * t567 + t455;
t89 = t445 * t472 + t552 * t567 + t524;
t88 = -t138 * t391 - t139 * t393 + t185 * t412;
t86 = -t138 * t469 - t139 * t470 - t185 * t467;
t79 = -t172 * t393 - t249 * t372 + (t173 + t262) * t391 + t546 * t370 + t545;
t78 = t155 * t567 + (-t131 * t447 + t132 * t445) * t446;
t77 = t154 * t567 + (-t129 * t447 + t130 * t445) * t446;
t76 = -t131 * t391 - t132 * t393 + t155 * t412;
t75 = -t129 * t391 - t130 * t393 + t154 * t412;
t72 = -t131 * t469 - t132 * t470 - t155 * t467;
t71 = -t129 * t469 - t130 * t470 - t154 * t467;
t56 = (t445 * t553 + t447 * t552) * t446 + t475;
t55 = -t184 * t467 + t222 * t351 + t307 * t542 - t470 * t550 + t554;
t54 = -t210 * t469 + t276 * t305 - t351 * t549 + t467 * t553 + t551;
t47 = t552 * t412 + t548 * t406 + (t298 - t550) * t393 + t523 * t372 + t544;
t46 = t550 * t391 + t542 * t370 + (-t261 - t553) * t412 + t526 * t406 + t543;
t40 = -t183 * t470 + t221 * t307 - t305 * t548 + t469 * t552 + t555;
t37 = -t553 * t393 - t549 * t372 + (t262 + t552) * t391 + t525 * t370 + t545;
t23 = -t102 * t467 + t138 * t305 + t139 * t307 + t185 * t351 - t469 * t57 - t470 * t58;
t14 = t131 * t305 + t132 * t307 + t155 * t351 - t467 * t95 - t469 * t50 - t470 * t51;
t13 = t129 * t305 + t130 * t307 + t154 * t351 - t467 * t94 - t469 * t48 - t470 * t49;
t1 = [0; t535 * t594 + t491 * t593 + (m(3) * t390 + m(4) * t330 + t456) * t557 + (m(3) * t388 + m(4) * t329 - t595) * t558 + t599 * (0.2e1 * t258 + 0.2e1 * t259 + t491); (t111 * t56 + t118 * t89 + t119 * t90) * t502 + (t101 * t128 + t113 * t143 + t114 * t144) * t504 + (t142 * t182 + t175 * t192 + t176 * t193) * t506 + (((-t375 * t510 + t377 * t511 + t385 * t578 + t387 * t576) * t445 - (-t374 * t510 + t376 * t511 + t384 * t578 + t386 * t576) * t447) * t601 + (-t412 * t362 + t413 * t363 - t406 * t365 - t407 * t366 + t606 * t567) * t567 + ((-t409 * t510 + t410 * t511 + t420 * t578 + t421 * t576) * t567 + (t412 * t325 - t413 * t327 + t406 * t333 + t407 * t335 + t567 * t607) * t447 + (-t412 * t326 + t413 * t328 - t406 * t334 - t407 * t336 + t567 * t605) * t445) * t446 + t603) * t567 + ((t326 * t393 + t328 * t394 + t334 * t372 + t336 * t373 + t375 * t418 + t377 * t417 + t385 * t426 - t387 * t464 + t558 * t605) * t558 + (t362 * t393 + t363 * t394 + t365 * t372 + t366 * t373 + t409 * t418 + t410 * t417 + t420 * t426 - t421 * t464 + t558 * t606) * t567 + t604) * t558 + ((t325 * t391 + t327 * t392 + t333 * t370 + t335 * t371 + t374 * t416 + t376 * t415 + t384 * t424 - t386 * t463 + t557 * t607) * t557 + (-t362 * t391 - t363 * t392 - t365 * t370 - t366 * t371 - t409 * t416 - t410 * t415 - t420 * t424 + t421 * t463 + t557 * t606) * t567 + (-t374 * t418 - t376 * t417 - t384 * t426 + t386 * t464 - t325 * t393 - t327 * t394 - t333 * t372 - t335 * t373 - t326 * t391 - t328 * t392 - t334 * t370 - t336 * t371 - t375 * t416 - t377 * t415 - t385 * t424 + t387 * t463 + t605 * t557 + t607 * t558) * t558 - t602) * t557 + 0.2e1 * m(4) * ((-t337 * t567 + t447 * t501 - t508) * t265 + (t338 * t567 + t445 * t501 + t369) * t264 + ((t337 * t445 + t338 * t447) * t446 + t538) * t220) + 0.2e1 * m(3) * ((-t378 * t567 - t411 * t557) * (-t388 * t567 - t422 * t557) + (t379 * t567 - t411 * t558) * (t390 * t567 - t422 * t558) + (t378 * t445 + t379 * t447) * t601 * (t388 * t445 + t390 * t447)); 0; (t56 * t530 + (t445 * t90 - t447 * t89) * t529) * t591 + (t101 * t530 + (-t113 * t447 + t114 * t445) * t529) * t592 + (t142 * t530 + (-t175 * t447 + t176 * t445) * t529) * t593 + (t220 * t530 + (-t264 * t447 + t265 * t445) * t529) * t594; 0; t595 * t393 + t456 * t391 + (-m(5) * t283 - m(7) * t549 - t573) * t372 + (m(5) * t284 + m(7) * t548 + t572) * t370 + 0.2e1 * t599 * (t391 * t262 + t370 * t303 + t545); (t111 * t37 + t112 * t56 + t118 * t47 + t119 * t46 + t120 * t90 + t121 * t89) * m(7) + (t101 * t140 + t113 * t153 + t114 * t152 + t128 * t79 + t143 * t98 + t144 * t97) * m(6) + m(5) * (t135 * t182 + t142 * t197 + t145 * t193 + t146 * t192 + t175 * t215 + t176 * t214) + (t195 * t567 + (-t178 * t447 + t179 * t445) * t446 + t77 + t73) * t585 + (t196 * t567 + (-t180 * t447 + t181 * t445) * t446 + t78 + t74) * t584 + t602 * t583 + t604 * t582 + (t209 * t567 + (-t186 * t447 + t187 * t445) * t446 + t93 + t87) * t580 + t603 * t579 + t596 * t527 + t597 * t513 + t598 * t512; m(5) * t135 * t567 + t79 * t519 + t37 * t518 + ((-m(5) * t146 - m(6) * t98 - m(7) * t47) * t447 + (m(5) * t145 + m(6) * t97 + m(7) * t46) * t445) * t446; (t112 * t37 + t120 * t46 + t121 * t47) * t502 + (t140 * t79 + t152 * t97 + t153 * t98) * t504 + (t135 * t197 + t145 * t214 + t146 * t215) * t506 + (t84 + t88) * t406 + (t564 + t596) * t412 + (-t187 * t406 - t597) * t393 + (-t186 * t406 - t598) * t391 + (t180 * t391 + t181 * t393 - t196 * t412 - t70 - t76) * t372 + (t178 * t391 + t179 * t393 - t195 * t412 - t69 - t75) * t370; t473 - (m(7) * t183 + t575) * t470 - (-m(7) * t184 - t574) * t469 + (m(7) * t221 + t573) * t307 + (-m(7) * t222 - t572) * t305; (t100 * t128 + t101 * t174 + t109 * t144 + t110 * t143 + t113 * t191 + t114 * t190) * m(6) + (t111 * t40 + t117 * t56 + t118 * t55 + t119 * t54 + t133 * t90 + t134 * t89) * m(7) + (t445 * t14 / 0.2e1 - t447 * t13 / 0.2e1) * t446 + t33 * t581 + t30 * t586 + t29 * t587 + t93 * t588 + t78 * t589 + t77 * t590 + t468 + t23 * t527; t100 * t519 + t40 * t518 + ((-m(6) * t110 - m(7) * t55) * t447 + (m(6) * t109 + m(7) * t54) * t445) * t446; (t100 * t140 + t109 * t152 + t110 * t153 + t174 * t79 + t190 * t97 + t191 * t98) * m(6) + (t112 * t40 + t117 * t37 + t120 * t54 + t121 * t55 + t133 * t46 + t134 * t47) * m(7) + t23 * t579 + t86 * t580 + t24 * t581 + t14 * t582 + t13 * t583 + t72 * t584 + t71 * t585 + t16 * t586 + t15 * t587 + t88 * t588 + t76 * t589 + t75 * t590 + t462; (t117 * t40 + t133 * t54 + t134 * t55) * t502 + t307 * t72 - t470 * t14 + t305 * t71 - t469 * t13 + t351 * t86 - t467 * t23 + (t100 * t174 + t109 * t190 + t110 * t191) * t504 + t479; t473; (t103 * t119 + t104 * t118 + t111 * t96 + t165 * t56 + t188 * t90 + t189 * t89) * m(7) + t468; (t96 * t530 + (t103 * t445 - t104 * t447) * t529) * t591; (t103 * t120 + t104 * t121 + t112 * t96 + t165 * t37 + t188 * t46 + t189 * t47) * m(7) + t462; (t103 * t133 + t104 * t134 + t117 * t96 + t165 * t40 + t188 * t54 + t189 * t55) * m(7) + t479; (t103 * t188 + t104 * t189 + t165 * t96) * t502 + t479;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;