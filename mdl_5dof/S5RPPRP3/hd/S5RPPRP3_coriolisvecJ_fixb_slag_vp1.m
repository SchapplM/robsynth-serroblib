% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:59
% DurationCPUTime: 12.74s
% Computational Cost: add. (6104->403), mult. (7383->501), div. (0->0), fcn. (5554->6), ass. (0->237)
t214 = qJ(1) + pkin(7);
t211 = sin(t214);
t212 = cos(t214);
t218 = cos(qJ(4));
t216 = sin(qJ(4));
t356 = Icges(5,4) * t216;
t264 = Icges(5,2) * t218 + t356;
t101 = -Icges(5,6) * t211 + t212 * t264;
t354 = Icges(6,4) * t216;
t263 = Icges(6,2) * t218 + t354;
t99 = -Icges(6,6) * t211 + t212 * t263;
t450 = t101 + t99;
t457 = Icges(5,3) + Icges(6,3);
t353 = Icges(6,4) * t218;
t265 = Icges(6,1) * t216 + t353;
t103 = -Icges(6,5) * t211 + t212 * t265;
t355 = Icges(5,4) * t218;
t266 = Icges(5,1) * t216 + t355;
t105 = -Icges(5,5) * t211 + t212 * t266;
t449 = t103 + t105;
t261 = Icges(6,5) * t216 + Icges(6,6) * t218;
t262 = Icges(5,5) * t216 + Icges(5,6) * t218;
t455 = t261 + t262;
t100 = Icges(5,6) * t212 + t211 * t264;
t98 = Icges(6,6) * t212 + t211 * t263;
t451 = t100 + t98;
t343 = t211 * t218;
t181 = Icges(6,4) * t343;
t344 = t211 * t216;
t351 = Icges(6,5) * t212;
t102 = Icges(6,1) * t344 + t181 + t351;
t182 = Icges(5,4) * t343;
t352 = Icges(5,5) * t212;
t104 = Icges(5,1) * t344 + t182 + t352;
t452 = t102 + t104;
t167 = -Icges(6,2) * t216 + t353;
t169 = -Icges(5,2) * t216 + t355;
t445 = t167 + t169;
t171 = Icges(6,1) * t218 - t354;
t173 = Icges(5,1) * t218 - t356;
t456 = t171 + t173;
t453 = t455 * t211 + t457 * t212;
t432 = t449 * t216 + t450 * t218;
t454 = t265 + t266;
t415 = -t457 * t211 + t455 * t212;
t444 = t456 * t216 + t445 * t218;
t412 = t452 * t216 + t451 * t218;
t416 = t453 * t211;
t163 = Icges(6,5) * t218 - Icges(6,6) * t216;
t116 = t211 * t163;
t165 = Icges(5,5) * t218 - Icges(5,6) * t216;
t118 = t211 * t165;
t448 = t116 + t118;
t447 = (t263 + t264) * qJD(4);
t446 = t454 * qJD(4);
t345 = t165 * t212;
t346 = t163 * t212;
t443 = -t345 - t346;
t442 = t432 * t212;
t441 = t453 * t212 + t451 * t343 + t452 * t344;
t419 = -t415 * t212 - t450 * t343 - t449 * t344;
t440 = -t412 * t212 + t416;
t417 = -t415 * t211 + t442;
t439 = t444 * t211 - t443;
t438 = t444 * t212 - t448;
t121 = t167 * t212;
t123 = t169 * t212;
t437 = (-t121 - t123) * qJD(4) + t451 * qJD(1);
t318 = qJD(4) * t211;
t436 = -t450 * qJD(1) - t445 * t318;
t125 = t171 * t212;
t127 = t173 * t212;
t435 = (-t125 - t127) * qJD(4) + (t454 * t211 + t351 + t352) * qJD(1);
t124 = t171 * t211;
t126 = t173 * t211;
t434 = (t124 + t126) * qJD(4) + t449 * qJD(1);
t433 = t444 * qJD(1) - t455 * qJD(4);
t431 = t447 * t218 + t446 * t216 + (t445 * t216 - t456 * t218) * qJD(4) + (t163 + t165) * qJD(1);
t430 = t439 * qJD(1);
t429 = (t417 * t211 + t440 * t212) * qJD(4);
t428 = (t419 * t211 + t441 * t212) * qJD(4);
t427 = t438 * qJD(1);
t426 = t428 + t430;
t425 = -t427 + t429;
t424 = -t412 * qJD(4) + t436 * t216 + t434 * t218;
t423 = t432 * qJD(4) - t437 * t216 + t435 * t218;
t422 = t433 * t211 + t431 * t212;
t421 = -t431 * t211 + t433 * t212;
t414 = t450 * t216 - t449 * t218;
t413 = t451 * t216 - t452 * t218;
t288 = -rSges(4,2) * t212 + t211 * rSges(4,3);
t219 = cos(qJ(1));
t213 = t219 * pkin(1);
t396 = t212 * pkin(2) + t211 * qJ(3);
t305 = t213 + t396;
t407 = t288 + t305;
t411 = t415 * qJD(1);
t410 = t453 * qJD(1);
t221 = qJD(1) ^ 2;
t409 = t412 * qJD(1) + t448 * qJD(4) + t411;
t408 = -t432 * qJD(1) + t443 * qJD(4) + t410;
t274 = rSges(6,1) * t216 + rSges(6,2) * t218;
t373 = pkin(4) * t216;
t242 = t274 + t373;
t300 = t218 * t318;
t319 = qJD(1) * t212;
t304 = t216 * t319;
t406 = t300 + t304;
t405 = t413 * qJD(4) - t434 * t216 + t436 * t218 + t410;
t404 = t414 * qJD(4) + t435 * t216 + t437 * t218 + t411;
t328 = t167 + t265;
t329 = -t263 + t171;
t403 = (t216 * t328 - t218 * t329) * qJD(1);
t326 = t169 + t266;
t327 = -t264 + t173;
t402 = (t216 * t326 - t218 * t327) * qJD(1);
t400 = 0.2e1 * qJD(4);
t187 = pkin(4) * t344;
t215 = -qJ(5) - pkin(6);
t368 = pkin(6) + t215;
t397 = -t368 * t212 + t187;
t289 = t212 * rSges(3,1) - rSges(3,2) * t211;
t395 = t213 + t289;
t303 = t218 * t319;
t320 = qJD(1) * t211;
t193 = qJD(5) * t212;
t331 = pkin(4) * t300 + t193;
t394 = -rSges(6,1) * t406 - rSges(6,2) * t303 - pkin(4) * t304 - t215 * t320 - t331;
t217 = sin(qJ(1));
t375 = pkin(1) * t217;
t393 = rSges(4,3) * t212 - t375;
t335 = t105 + t123;
t339 = t101 - t127;
t383 = t216 * t339 - t218 * t335;
t336 = -Icges(5,2) * t344 + t104 + t182;
t340 = t100 - t126;
t382 = t216 * t340 - t218 * t336;
t337 = t103 + t121;
t357 = -t125 + t99;
t381 = t216 * t357 - t218 * t337;
t338 = -Icges(6,2) * t344 + t102 + t181;
t358 = t124 - t98;
t380 = -t216 * t358 - t218 * t338;
t379 = t211 / 0.2e1;
t378 = -t212 / 0.2e1;
t376 = -rSges(6,3) - pkin(2);
t374 = pkin(2) * t211;
t372 = pkin(4) * t218;
t371 = pkin(6) * t211;
t370 = -qJD(1) / 0.2e1;
t178 = rSges(6,1) * t218 - rSges(6,2) * t216;
t130 = t178 * t212;
t207 = t212 * rSges(6,3);
t311 = t212 * t372;
t315 = qJD(5) * t211;
t287 = -qJD(4) * t311 + t315;
t367 = -qJD(4) * t130 + t287 + (t211 * t274 + t207 + t397) * qJD(1);
t316 = qJD(4) * t216;
t366 = -(-rSges(6,2) * t316 - rSges(6,3) * qJD(1)) * t211 - pkin(6) * t320 + t394;
t208 = t212 * rSges(5,3);
t179 = rSges(5,1) * t218 - rSges(5,2) * t216;
t133 = t179 * t318;
t197 = t212 * qJ(3);
t137 = -t197 + t374;
t194 = qJD(3) * t211;
t205 = t211 * rSges(5,3);
t275 = rSges(5,1) * t216 + rSges(5,2) * t218;
t109 = t275 * t212 - t205;
t280 = -t371 - t375;
t251 = t109 + t280;
t43 = t133 + t194 + (-t137 + t251) * qJD(1);
t362 = t212 * t43;
t148 = t275 * qJD(4);
t314 = qJD(1) * qJD(3);
t325 = qJ(3) * t319 + t194;
t332 = qJD(1) * (-pkin(2) * t320 + t325) + t211 * t314;
t227 = t221 * t280 + t332;
t317 = qJD(4) * t212;
t306 = rSges(5,1) * t406 + rSges(5,2) * t303;
t75 = (-rSges(5,2) * t316 - rSges(5,3) * qJD(1)) * t211 + t306;
t25 = t148 * t317 + (t75 + t133) * qJD(1) + t227;
t361 = t25 * t212;
t195 = qJD(3) * t212;
t113 = qJD(1) * t396 - t195;
t189 = t212 * t314;
t279 = -pkin(6) * t212 - t213;
t236 = t221 * t279 + t189;
t302 = t179 * t317;
t131 = t179 * t212;
t73 = -qJD(4) * t131 + (t211 * t275 + t208) * qJD(1);
t26 = -t148 * t318 + (-t113 - t73 + t302) * qJD(1) + t236;
t360 = t26 * t211;
t204 = t211 * rSges(6,3);
t333 = t211 * t368 + t242 * t212 - t204;
t243 = t280 + t333;
t286 = t178 * t318 + t194 + t331;
t31 = (-t137 + t243) * qJD(1) + t286;
t359 = t31 * t178;
t106 = rSges(6,1) * t344 + rSges(6,2) * t343 + t207;
t334 = -t106 - t397;
t324 = rSges(4,2) * t320 + rSges(4,3) * t319;
t134 = qJD(1) * t137;
t323 = t194 - t134;
t322 = qJD(1) * t261;
t321 = qJD(1) * t262;
t313 = -rSges(5,3) - pkin(2) - pkin(6);
t310 = t221 * t375;
t309 = t221 * t213;
t107 = rSges(5,1) * t344 + rSges(5,2) * t343 + t208;
t301 = t211 * t316;
t296 = -t318 / 0.2e1;
t294 = -t317 / 0.2e1;
t291 = t178 + t372;
t281 = -t374 - t375;
t278 = -t195 + t287;
t139 = rSges(3,1) * t211 + rSges(3,2) * t212;
t250 = t396 - t279;
t44 = -t302 - t195 + (t107 + t250) * qJD(1);
t271 = t211 * t43 - t212 * t44;
t270 = qJD(4) * t291;
t257 = -t107 * t211 - t109 * t212;
t252 = t197 + t281;
t147 = t274 * qJD(4);
t249 = (t373 * qJD(4) + t147) * qJD(4);
t226 = -t211 * t75 + t212 * t73 + (-t107 * t212 + t109 * t211) * qJD(1);
t129 = t179 * t211;
t128 = t178 * t211;
t51 = -t309 + t189 + (-qJD(1) * t288 - t113) * qJD(1);
t50 = qJD(1) * t324 - t310 + t332;
t45 = qJD(4) * t257 + qJD(2);
t32 = -t178 * t317 + (t250 - t334) * qJD(1) + t278;
t30 = qJD(2) + (t211 * t334 - t212 * t333) * qJD(4);
t16 = t226 * qJD(4);
t15 = -t249 * t211 + (t212 * t270 - t113 - t315 - t367) * qJD(1) + t236;
t14 = t249 * t212 + (t211 * t270 + t193 - t366) * qJD(1) + t227;
t1 = (t367 * t212 + t366 * t211 + (t211 * t333 + t212 * t334) * qJD(1)) * qJD(4);
t2 = [m(3) * ((-t139 * t221 - t310) * t395 + (-t309 + (-0.2e1 * t289 - t213 + t395) * t221) * (-t139 - t375)) + (((t417 + t441 - t442) * t212 + ((-t412 + t415) * t212 + t416 - t440 + t419) * t211) * qJD(4) + t430) * t296 + (-t444 * qJD(4) + t447 * t216 - t446 * t218) * qJD(1) + (t15 * (t211 * t215 - t204 + t252) - t31 * t278 + t14 * (t187 + t305 + t106) + t32 * (-rSges(6,2) * t301 + t325 - t394) + (qJD(4) * t359 - t14 * t215 + t15 * t242) * t212 + ((-t217 * t32 - t219 * t31) * pkin(1) + t31 * (t215 + t376) * t212 + (t31 * (-qJ(3) - t242) + t32 * t376) * t211) * qJD(1) - (qJD(1) * t243 - t134 + t286 - t31) * t32) * m(6) + (t26 * (-t205 + t252 - t371) + t43 * t195 + t25 * (t305 + t107) + t44 * (-rSges(5,2) * t301 + t306 + t325) + (qJD(4) * t179 * t43 + t25 * pkin(6) + t26 * t275) * t212 + ((-t217 * t44 - t219 * t43) * pkin(1) + t313 * t362 + (t43 * (-qJ(3) - t275) + t44 * t313) * t211) * qJD(1) - (qJD(1) * t251 + t133 + t323 - t43) * t44) * m(5) + (t51 * (t197 + (rSges(4,2) - pkin(2)) * t211 + t393) + t50 * t407 + (-t323 + t324 + t325 + (-rSges(4,2) * t211 + t281 - t393) * qJD(1)) * (qJD(1) * t407 - t195)) * m(4) + ((t415 * t211 ^ 2 + ((t412 + t415) * t212 - t416 + t419) * t212) * qJD(4) + t425 + t427) * t294 + (t422 + t423 + t426) * t318 / 0.2e1 + (t414 * qJD(1) + t421 + t424) * t317 / 0.2e1 + (t438 * t212 + (-t413 + t439) * t211) * qJD(4) * t370; m(5) * t16 + m(6) * t1; 0.2e1 * (t14 * t378 + t15 * t379) * m(6) + 0.2e1 * (t360 / 0.2e1 - t361 / 0.2e1) * m(5) + 0.2e1 * (t378 * t50 + t379 * t51) * m(4); ((((-t340 + t358) * t212 + (t339 + t357) * t211) * t218 + ((-t336 - t338) * t212 + (t335 + t337) * t211) * t216) * qJD(4) + ((-t326 - t328) * t218 + (-t327 - t329) * t216) * qJD(1)) * t370 + (t424 * t212 + t423 * t211 + (t413 * t211 + t414 * t212) * qJD(1)) * qJD(1) / 0.2e1 + ((-t318 * t346 - t322) * t211 + (t403 + (t380 * t212 + (t116 - t381) * t211) * qJD(4)) * t212 + (-t318 * t345 - t321) * t211 + (t402 + (t382 * t212 + (t118 - t383) * t211) * qJD(4)) * t212) * t296 + ((t116 * t317 - t322) * t212 + (-t403 + (t381 * t211 + (-t346 - t380) * t212) * qJD(4)) * t211 + (t118 * t317 - t321) * t212 + (-t402 + (t383 * t211 + (-t345 - t382) * t212) * qJD(4)) * t211) * t294 + (-((t32 * t274 + (-t130 - t311) * t30) * t212 + (-t31 * t242 + (-pkin(4) * t343 - t128) * t30) * t211) * qJD(4) + (-t1 * t333 - t14 * t291 + t32 * t147 + t30 * t367) * t212 + (t15 * t291 + t31 * (-pkin(4) * t316 - t147) + t1 * t334 + t30 * t366) * t211 + (-t31 * t130 - t32 * t128 + (t30 * t334 + t359) * t212 + (t32 * t178 + t30 * t333) * t211) * qJD(1)) * m(6) + (t16 * t257 + t45 * t226 - t271 * t148 + (t360 - t361 + (t211 * t44 + t362) * qJD(1)) * t179 - (t129 * t44 + t131 * t43) * qJD(1) - (t45 * (-t129 * t211 - t131 * t212) - t271 * t275) * qJD(4)) * m(5) + (t422 * qJD(1) + ((t417 * qJD(1) + t405 * t212) * t212 + (t408 * t211 - t440 * qJD(1) + (-t404 + t409) * t212) * t211) * t400) * t379 + (t421 * qJD(1) + ((t419 * qJD(1) + t409 * t212) * t212 + (t404 * t211 - t441 * qJD(1) + (-t405 + t408) * t212) * t211) * t400) * t212 / 0.2e1 - (t426 + t428) * t320 / 0.2e1 + (t425 + t429) * t319 / 0.2e1; m(6) * (t14 * t211 + t15 * t212);];
tauc = t2(:);
