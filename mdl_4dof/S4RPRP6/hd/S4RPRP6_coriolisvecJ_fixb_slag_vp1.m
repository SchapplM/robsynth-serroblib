% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:10
% DurationCPUTime: 9.87s
% Computational Cost: add. (2769->328), mult. (6991->435), div. (0->0), fcn. (5354->4), ass. (0->186)
t455 = Icges(4,4) + Icges(5,4);
t452 = Icges(4,1) + Icges(5,1);
t450 = Icges(4,2) + Icges(5,2);
t451 = Icges(4,5) + Icges(5,5);
t449 = Icges(4,6) + Icges(5,6);
t208 = cos(qJ(3));
t454 = t455 * t208;
t206 = sin(qJ(3));
t453 = t455 * t206;
t442 = t450 * t208 + t453;
t432 = t452 * t206 + t454;
t448 = Icges(4,3) + Icges(5,3);
t207 = sin(qJ(1));
t321 = t207 * t208;
t447 = t455 * t321;
t444 = t451 * t206 + t449 * t208;
t209 = cos(qJ(1));
t446 = t451 * t209;
t440 = t207 * t442 + t449 * t209;
t439 = -t207 * t449 + t442 * t209;
t322 = t206 * t207;
t437 = t322 * t452 + t446 + t447;
t436 = -t207 * t451 + t432 * t209;
t429 = -t206 * t450 + t454;
t445 = t208 * t452 - t453;
t441 = t444 * t207 + t448 * t209;
t443 = -t206 * t449 + t208 * t451;
t375 = -t448 * t207 + t444 * t209;
t410 = t206 * t436 + t208 * t439;
t438 = t410 * t209;
t435 = t429 * t209;
t434 = t445 * t207;
t433 = t445 * t209;
t427 = t206 * t445 + t208 * t429;
t414 = t437 * t206 + t208 * t440;
t376 = t441 * t207;
t399 = t443 * t207;
t431 = t442 * qJD(3);
t430 = t432 * qJD(3);
t391 = t443 * t209;
t425 = t209 * t441 + t440 * t321 + t437 * t322;
t379 = -t209 * t375 - t321 * t439 - t322 * t436;
t424 = -t209 * t414 + t376;
t401 = -t375 * t207 + t438;
t378 = t206 * t440 - t208 * t437;
t423 = t207 * t427 + t391;
t422 = t209 * t427 - t399;
t421 = qJD(1) * t440 - qJD(3) * t435;
t299 = qJD(3) * t207;
t420 = -qJD(1) * t439 - t299 * t429;
t419 = -t433 * qJD(3) + (t207 * t432 + t446) * qJD(1);
t418 = qJD(1) * t436 + qJD(3) * t434;
t192 = t209 * qJ(2);
t156 = pkin(1) * t207 - t192;
t189 = qJD(2) * t207;
t301 = qJD(1) * t209;
t305 = qJ(2) * t301 + t189;
t417 = qJD(1) * t156 + t305;
t416 = (t322 * t450 - t437 - t447) * t209 + (t435 + t436) * t207;
t415 = (t434 - t440) * t209 + (-t433 + t439) * t207;
t377 = t206 * t439 - t208 * t436;
t413 = -t429 - t432;
t412 = t445 - t442;
t411 = t427 * qJD(1) - qJD(3) * t444;
t409 = -t189 + t417;
t408 = t375 * qJD(1);
t407 = t441 * qJD(1);
t406 = t431 * t208 + t430 * t206 + (t429 * t206 - t208 * t445) * qJD(3) + t443 * qJD(1);
t405 = (t206 * t413 + t208 * t412) * qJD(1);
t404 = pkin(3) * t208;
t400 = t423 * qJD(1);
t398 = qJD(1) * t414 + t399 * qJD(3) + t408;
t397 = -qJD(1) * t410 - t391 * qJD(3) + t407;
t396 = (t207 * t401 + t209 * t424) * qJD(3);
t395 = (t379 * t207 + t209 * t425) * qJD(3);
t260 = rSges(5,1) * t206 + rSges(5,2) * t208;
t353 = pkin(3) * t206;
t234 = t260 + t353;
t284 = t208 * t299;
t288 = t206 * t301;
t394 = t284 + t288;
t393 = t444 * qJD(1);
t392 = t422 * qJD(1);
t389 = qJD(3) * t377 + t206 * t419 + t208 * t421 + t408;
t388 = qJD(3) * t378 - t206 * t418 + t208 * t420 + t407;
t387 = -t206 * t415 + t208 * t416;
t386 = 0.2e1 * qJD(3);
t385 = t395 + t400;
t384 = -t392 + t396;
t383 = -qJD(3) * t414 + t206 * t420 + t208 * t418;
t382 = t410 * qJD(3) - t206 * t421 + t208 * t419;
t381 = t411 * t207 + t406 * t209;
t380 = -t406 * t207 + t411 * t209;
t161 = t209 * pkin(1) + t207 * qJ(2);
t272 = -rSges(3,2) * t209 + t207 * rSges(3,3);
t370 = t161 + t272;
t182 = pkin(3) * t322;
t205 = -qJ(4) - pkin(5);
t347 = pkin(5) + t205;
t369 = -t347 * t209 + t182;
t287 = t208 * t301;
t302 = qJD(1) * t207;
t188 = qJD(4) * t209;
t307 = pkin(3) * t284 + t188;
t367 = -rSges(5,1) * t394 - rSges(5,2) * t287 - pkin(3) * t288 - t205 * t302 - t307;
t358 = t207 / 0.2e1;
t357 = -t209 / 0.2e1;
t355 = rSges(3,2) - pkin(1);
t354 = -rSges(5,3) - pkin(1);
t352 = pkin(3) * t209;
t351 = pkin(5) * t207;
t350 = pkin(5) * qJD(1) ^ 2;
t349 = -qJD(1) / 0.2e1;
t159 = rSges(5,1) * t208 - rSges(5,2) * t206;
t126 = t159 * t209;
t202 = t209 * rSges(5,3);
t297 = qJD(4) * t207;
t298 = qJD(3) * t209;
t271 = -t298 * t404 + t297;
t346 = -qJD(3) * t126 + t271 + (t207 * t260 + t202 + t369) * qJD(1);
t300 = qJD(3) * t206;
t345 = -(-rSges(5,2) * t300 - rSges(5,3) * qJD(1)) * t207 - pkin(5) * t302 + t367;
t343 = rSges(3,3) * t209;
t203 = t209 * rSges(4,3);
t160 = rSges(4,1) * t208 - rSges(4,2) * t206;
t129 = t160 * t299;
t200 = t207 * rSges(4,3);
t261 = rSges(4,1) * t206 + rSges(4,2) * t208;
t107 = t261 * t209 - t200;
t275 = t107 - t351;
t41 = t129 + t189 + (-t156 + t275) * qJD(1);
t342 = t209 * t41;
t139 = t261 * qJD(3);
t294 = qJD(1) * qJD(2);
t312 = qJD(1) * (-pkin(1) * t302 + t305) + t207 * t294;
t241 = -t207 * t350 + t312;
t289 = rSges(4,1) * t394 + rSges(4,2) * t287;
t71 = (-rSges(4,2) * t300 - rSges(4,3) * qJD(1)) * t207 + t289;
t23 = t139 * t298 + (t71 + t129) * qJD(1) + t241;
t341 = t23 * t209;
t190 = qJD(2) * t209;
t111 = qJD(1) * t161 - t190;
t184 = t209 * t294;
t268 = -t209 * t350 + t184;
t286 = t160 * t298;
t127 = t160 * t209;
t69 = -qJD(3) * t127 + (t207 * t261 + t203) * qJD(1);
t24 = -t139 * t299 + (-t111 - t69 + t286) * qJD(1) + t268;
t340 = t24 * t207;
t199 = t207 * rSges(5,3);
t313 = t206 * t352 + t207 * t347 + t209 * t260 - t199;
t267 = t313 - t351;
t270 = t159 * t299 + t189 + t307;
t29 = (-t156 + t267) * qJD(1) + t270;
t339 = t29 * t159;
t104 = rSges(5,1) * t322 + rSges(5,2) * t321 + t202;
t314 = -t104 - t369;
t304 = rSges(3,2) * t302 + rSges(3,3) * t301;
t293 = -rSges(4,3) - pkin(1) - pkin(5);
t105 = rSges(4,1) * t322 + rSges(4,2) * t321 + t203;
t285 = t206 * t299;
t280 = -t299 / 0.2e1;
t278 = -t298 / 0.2e1;
t277 = t298 / 0.2e1;
t276 = t159 + t404;
t274 = pkin(5) * t209 + t161;
t263 = -t190 + t271;
t42 = -t286 - t190 + (t105 + t274) * qJD(1);
t259 = t207 * t41 - t209 * t42;
t258 = qJD(3) * t276;
t138 = t260 * qJD(3);
t240 = (t353 * qJD(3) + t138) * qJD(3);
t43 = (-t105 * t207 - t107 * t209) * qJD(3);
t219 = t207 * t314 - t209 * t313;
t157 = rSges(3,2) * t207 + t343;
t125 = t160 * t207;
t124 = t159 * t207;
t77 = qJD(1) * t370 - t190;
t76 = t189 + (-t156 + t157) * qJD(1);
t67 = t184 + (-qJD(1) * t272 - t111) * qJD(1);
t66 = qJD(1) * t304 + t312;
t30 = -t159 * t298 + (t274 - t314) * qJD(1) + t263;
t28 = t219 * qJD(3);
t14 = -t240 * t207 + (t209 * t258 - t111 - t297 - t346) * qJD(1) + t268;
t13 = t240 * t209 + (t207 * t258 + t188 - t345) * qJD(1) + t241;
t1 = [(((t376 + t379 - t424) * t207 + (-t438 + (-t414 + t375) * t207 + t401 + t425) * t209) * qJD(3) + t400) * t280 + (t14 * (t205 * t207 - t156 - t199) - t29 * t263 + t13 * (t182 + t104 + t161) + (qJD(3) * t339 - t13 * t205 + t14 * t234) * t209 + (-rSges(5,2) * t285 - t270 + t29 - t367 + t417) * t30) * m(5) + (t24 * (-t156 - t200 - t351) + t41 * t190 + t23 * (t105 + t161) + (qJD(3) * t160 * t41 + t23 * pkin(5) + t24 * t261) * t209 + (-rSges(4,2) * t285 - t129 + t289 + t409 + t41) * t42) * m(4) + (t67 * (t207 * t355 + t192 + t343) + t76 * t190 + t66 * t370 + (t304 + t76 + t409) * t77) * m(3) + ((((t414 + t375) * t209 - t376 + t379) * t209 + t375 * t207 ^ 2) * qJD(3) + t384 + t392) * t278 + (t380 + t383) * t277 + (t381 + t382 + t385) * t299 / 0.2e1 + ((-t378 + t423) * t207 + t422 * t209) * qJD(3) * t349 + (t377 * t277 - t430 * t208 + t431 * t206 - t427 * qJD(3) + (-t267 * t30 + t29 * (t205 + t354) * t209 + (t29 * (-qJ(2) - t234) + t30 * t354) * t207) * m(5) + (t293 * t342 + (t41 * (-qJ(2) - t261) + t42 * t293) * t207 - t275 * t42) * m(4) + (t76 * t355 * t209 + (t76 * (-rSges(3,3) - qJ(2)) - t77 * pkin(1)) * t207 - t157 * t77) * m(3)) * qJD(1); 0.2e1 * (t13 * t357 + t14 * t358) * m(5) + 0.2e1 * (t340 / 0.2e1 - t341 / 0.2e1) * m(4) + 0.2e1 * (t66 * t357 + t358 * t67) * m(3); ((t206 * t416 + t208 * t415) * qJD(3) + (-t206 * t412 + t208 * t413) * qJD(1)) * t349 + (t383 * t209 + t382 * t207 + (t378 * t207 + t377 * t209) * qJD(1)) * qJD(1) / 0.2e1 + ((-t299 * t391 - t393) * t207 + ((t207 * t399 + t387) * qJD(3) - t405) * t209) * t280 + ((t298 * t399 - t393) * t209 + ((-t209 * t391 - t387) * qJD(3) + t405) * t207) * t278 + ((-t13 * t276 + t30 * t138 + t28 * t346) * t209 + (-t29 * t138 + t14 * t276 + t28 * t345) * t207 + (-t29 * t126 - t30 * t124 + (t28 * t314 + t339) * t209 + (t30 * t159 + t28 * t313) * t207) * qJD(1) + (-(t30 * t260 + (-t208 * t352 - t126) * t28) * t209 - (-t29 * t234 + (-pkin(3) * t321 - t124) * t28) * t207 - t29 * t182 + (t346 * t209 + t345 * t207 + (t207 * t313 + t209 * t314) * qJD(1)) * t219) * qJD(3)) * m(5) + (0.2e1 * t43 * (-t207 * t71 + t209 * t69 + (-t105 * t209 + t107 * t207) * qJD(1)) - t259 * t139 + (t340 - t341 + (t207 * t42 + t342) * qJD(1)) * t160 - (t125 * t42 + t127 * t41) * qJD(1) - (t43 * (-t125 * t207 - t127 * t209) - t259 * t261) * qJD(3)) * m(4) + (t381 * qJD(1) + ((t401 * qJD(1) + t388 * t209) * t209 + (t397 * t207 - t424 * qJD(1) + (-t389 + t398) * t209) * t207) * t386) * t358 + (t380 * qJD(1) + ((t379 * qJD(1) + t398 * t209) * t209 + (t389 * t207 - t425 * qJD(1) + (-t388 + t397) * t209) * t207) * t386) * t209 / 0.2e1 - (t385 + t395) * t302 / 0.2e1 + (t384 + t396) * t301 / 0.2e1; m(5) * (t13 * t207 + t14 * t209);];
tauc = t1(:);
