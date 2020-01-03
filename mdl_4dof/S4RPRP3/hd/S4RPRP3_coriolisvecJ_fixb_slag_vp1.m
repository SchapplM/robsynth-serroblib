% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:44
% DurationCPUTime: 10.85s
% Computational Cost: add. (5495->372), mult. (6631->469), div. (0->0), fcn. (5130->6), ass. (0->222)
t429 = Icges(4,3) + Icges(5,3);
t192 = qJ(1) + pkin(6);
t187 = sin(t192);
t196 = cos(qJ(3));
t298 = t187 * t196;
t194 = sin(qJ(3));
t299 = t187 * t194;
t188 = cos(t192);
t308 = Icges(5,6) * t188;
t87 = Icges(5,4) * t298 - Icges(5,2) * t299 - t308;
t309 = Icges(4,6) * t188;
t89 = Icges(4,4) * t298 - Icges(4,2) * t299 - t309;
t419 = t87 + t89;
t139 = Icges(5,5) * t196 - Icges(5,6) * t194;
t141 = Icges(4,5) * t196 - Icges(4,6) * t194;
t414 = t139 + t141;
t160 = Icges(5,4) * t299;
t312 = Icges(5,5) * t188;
t91 = Icges(5,1) * t298 - t160 - t312;
t161 = Icges(4,4) * t299;
t313 = Icges(4,5) * t188;
t93 = Icges(4,1) * t298 - t161 - t313;
t426 = t91 + t93;
t314 = Icges(5,4) * t194;
t147 = Icges(5,1) * t196 - t314;
t92 = Icges(5,5) * t187 + t147 * t188;
t315 = Icges(4,4) * t194;
t149 = Icges(4,1) * t196 - t315;
t94 = Icges(4,5) * t187 + t149 * t188;
t417 = t92 + t94;
t142 = Icges(5,2) * t196 + t314;
t144 = Icges(4,2) * t196 + t315;
t428 = t142 + t144;
t189 = Icges(5,4) * t196;
t146 = Icges(5,1) * t194 + t189;
t190 = Icges(4,4) * t196;
t148 = Icges(4,1) * t194 + t190;
t424 = t146 + t148;
t427 = t429 * t188;
t394 = t427 + (Icges(4,6) + Icges(5,6)) * t299 + (-Icges(4,5) - Icges(5,5)) * t298;
t393 = t429 * t187 + t414 * t188;
t423 = t147 + t149;
t233 = -Icges(5,2) * t194 + t189;
t234 = -Icges(4,2) * t194 + t190;
t422 = t233 + t234;
t421 = t419 * t194;
t420 = t417 * t298;
t88 = Icges(5,6) * t187 + t188 * t233;
t90 = Icges(4,6) * t187 + t188 * t234;
t418 = t88 + t90;
t410 = t428 * t194 - t424 * t196;
t397 = -t426 * t196 + t421;
t416 = t422 * qJD(3);
t415 = t423 * qJD(3);
t138 = Icges(5,5) * t194 + Icges(5,6) * t196;
t140 = Icges(4,5) * t194 + Icges(4,6) * t196;
t413 = t140 + t138;
t412 = t428 * qJD(3);
t411 = t424 * qJD(3);
t409 = t393 * t188 - t420;
t296 = t188 * t196;
t408 = t394 * t187 - t426 * t296;
t376 = t393 * t187 + t417 * t296;
t407 = -t397 * t187 + t394 * t188;
t383 = -t418 * t299 - t409;
t297 = t188 * t194;
t382 = -t419 * t297 - t408;
t381 = -t418 * t297 + t376;
t300 = t140 * t188;
t302 = t138 * t188;
t406 = -t410 * t187 - t300 - t302;
t301 = t140 * t187;
t303 = t138 * t187;
t405 = -t410 * t188 + t301 + t303;
t404 = t418 * t194;
t403 = -t412 * t188 + (-t422 * t187 + t308 + t309) * qJD(1);
t402 = t418 * qJD(1) - t412 * t187;
t401 = -t411 * t188 + (-t423 * t187 + t312 + t313) * qJD(1);
t400 = -t417 * qJD(1) + t411 * t187;
t399 = t415 * t196 - t416 * t194 + (-t424 * t194 - t196 * t428) * qJD(3) + t413 * qJD(1);
t398 = t417 * t196 - t404;
t396 = t410 * qJD(1) + t414 * qJD(3);
t395 = t405 * qJD(1);
t392 = (t381 * t187 - t382 * t188) * qJD(3);
t391 = (t383 * t187 - t407 * t188) * qJD(3);
t390 = t406 * qJD(1);
t389 = t390 + t391;
t388 = t392 + t395;
t387 = qJD(3) * t397 + t194 * t400 - t196 * t402;
t386 = t398 * qJD(3) + t401 * t194 + t196 * t403;
t385 = t187 * t396 + t188 * t399;
t384 = t187 * t399 - t188 * t396;
t380 = t426 * t194 + t419 * t196;
t379 = t417 * t194 + t418 * t196;
t378 = t413 * qJD(3);
t377 = t394 + t404;
t375 = t393 * qJD(1);
t199 = qJD(1) ^ 2;
t374 = -t379 * qJD(3) - t194 * t403 + t401 * t196 + t375;
t373 = t394 * qJD(1) + t380 * qJD(3) + t194 * t402 + t196 * t400;
t372 = qJD(1) * t397 - t378 * t187 + t375;
t371 = -t378 * t188 + (-t414 * t187 - t398 + t427) * qJD(1);
t370 = 0.2e1 * qJD(3);
t178 = t187 * rSges(5,3);
t352 = pkin(3) * t196;
t186 = pkin(2) + t352;
t193 = -qJ(4) - pkin(5);
t369 = rSges(5,1) * t296 - rSges(5,2) * t297 + t188 * t186 - t187 * t193 + t178;
t122 = t188 * pkin(2) + t187 * pkin(5);
t197 = cos(qJ(1));
t191 = t197 * pkin(1);
t261 = t122 + t191;
t179 = t187 * rSges(4,3);
t98 = rSges(4,1) * t296 - rSges(4,2) * t297 + t179;
t368 = t261 + t98;
t166 = t188 * t193;
t367 = rSges(5,1) * t298 - rSges(5,2) * t299 - t188 * rSges(5,3) + t187 * t186 + t166;
t170 = qJD(4) * t187;
t328 = t196 * rSges(5,2);
t150 = rSges(5,1) * t194 + t328;
t260 = pkin(3) * t194 + t150;
t279 = qJD(3) * t188;
t215 = -t260 * t279 + t170;
t282 = qJD(1) * t187;
t274 = t194 * t282;
t281 = qJD(1) * t188;
t366 = rSges(5,2) * t274 + rSges(5,3) * t281 + t170;
t253 = t188 * rSges(3,1) - rSges(3,2) * t187;
t365 = t191 + t253;
t319 = -t148 * t187 - t89;
t323 = -Icges(4,2) * t298 - t161 + t93;
t357 = -t194 * t323 + t196 * t319;
t321 = -t146 * t187 - t87;
t325 = -Icges(5,2) * t298 - t160 + t91;
t356 = -t194 * t325 + t196 * t321;
t355 = t187 / 0.2e1;
t354 = -t188 / 0.2e1;
t195 = sin(qJ(1));
t353 = pkin(1) * t195;
t350 = qJD(1) / 0.2e1;
t349 = pkin(2) - t186;
t169 = pkin(5) * t281;
t271 = t194 * t279;
t216 = -t196 * t282 - t271;
t278 = qJD(3) * t196;
t270 = t188 * t278;
t348 = -pkin(3) * t271 - t169 + (t187 * t349 - t166) * qJD(1) + rSges(5,1) * t216 - rSges(5,2) * t270 + t366;
t113 = t150 * t187;
t339 = rSges(5,1) * t196;
t153 = -rSges(5,2) * t194 + t339;
t171 = qJD(4) * t188;
t277 = pkin(3) * t299;
t290 = qJD(3) * t277 + t171;
t347 = -qJD(3) * t113 - t290 + ((-pkin(5) - t193) * t187 + t178 + (-t349 + t153) * t188) * qJD(1);
t184 = t188 * pkin(5);
t121 = pkin(2) * t187 - t184;
t342 = -t121 + t367;
t341 = -t122 + t369;
t340 = rSges(4,1) * t196;
t151 = rSges(4,1) * t194 + rSges(4,2) * t196;
t116 = t151 * t188;
t280 = qJD(3) * t187;
t273 = t151 * t280;
t42 = qJD(1) * t368 - t273;
t337 = t116 * t42;
t285 = rSges(4,2) * t299 + t188 * rSges(4,3);
t96 = rSges(4,1) * t298 - t285;
t268 = -t96 - t353;
t272 = t151 * t279;
t41 = -t272 + (-t121 + t268) * qJD(1);
t336 = t187 * t41;
t335 = t188 * t41;
t251 = -t342 - t353;
t28 = (-t121 + t251) * qJD(1) + t215;
t326 = t28 * t150;
t324 = -t142 * t188 + t92;
t322 = -t144 * t188 + t94;
t320 = -t146 * t188 - t88;
t318 = -t148 * t188 - t90;
t291 = rSges(4,2) * t274 + rSges(4,3) * t281;
t289 = -t142 + t147;
t288 = t146 + t233;
t287 = -t144 + t149;
t286 = t148 + t234;
t284 = qJD(1) * t139;
t283 = qJD(1) * t141;
t276 = t199 * t353;
t275 = t199 * t191;
t269 = -pkin(2) - t340;
t265 = -t280 / 0.2e1;
t262 = t279 / 0.2e1;
t259 = -t153 - t352;
t250 = qJD(1) * (-pkin(2) * t282 + t169) - t276;
t249 = -pkin(3) * t297 - t150 * t188;
t129 = t153 * qJD(3);
t246 = -pkin(3) * t278 - t129;
t120 = rSges(3,1) * t187 + rSges(3,2) * t188;
t243 = -rSges(4,2) * t194 + t340;
t242 = -t187 * t42 - t335;
t241 = t187 * t96 + t188 * t98;
t229 = (-t352 * qJD(3) - t129) * qJD(3);
t114 = t151 * t187;
t214 = -t194 * t324 + t196 * t320;
t213 = -t194 * t322 + t196 * t318;
t212 = (-t194 * t288 + t196 * t289) * qJD(1);
t211 = (-t194 * t286 + t196 * t287) * qJD(1);
t68 = rSges(4,1) * t216 - rSges(4,2) * t270 + t291;
t70 = -qJD(3) * t114 + (t188 * t243 + t179) * qJD(1);
t210 = t187 * t70 + t188 * t68 + (-t187 * t98 + t188 * t96) * qJD(1);
t130 = t243 * qJD(3);
t119 = qJD(1) * t121;
t118 = t122 * qJD(1);
t40 = qJD(3) * t241 + qJD(2);
t31 = -t275 - t130 * t279 + (-t118 - t70 + t273) * qJD(1);
t30 = -t130 * t280 + (t68 - t272) * qJD(1) + t250;
t29 = -t150 * t280 + (t261 + t341) * qJD(1) - t290;
t25 = qJD(2) + (t187 * t342 + t188 * t341) * qJD(3);
t16 = t210 * qJD(3);
t15 = -t275 + t229 * t188 + (t260 * t280 - t118 + t171 - t347) * qJD(1);
t14 = t229 * t187 + (t348 + t215) * qJD(1) + t250;
t1 = (t348 * t188 + t347 * t187 + (-t187 * t341 + t188 * t342) * qJD(1)) * qJD(3);
t2 = [m(3) * ((-t120 * t199 - t276) * t365 + (-t275 + (-0.2e1 * t253 - t191 + t365) * t199) * (-t120 - t353)) + ((t376 * t187 + ((t393 + t421) * t188 + t383 + t408 - t420) * t188) * qJD(3) + t395) * t262 + (-t410 * qJD(3) + t415 * t194 + t416 * t196) * qJD(1) + (t15 * (-t353 - t367) + t28 * t290 + t14 * (t191 + t369) + t29 * t366 + (t187 * t326 + t29 * (-t328 + (-rSges(5,1) - pkin(3)) * t194) * t188) * qJD(3) + ((-t195 * t29 - t197 * t28) * pkin(1) + (t28 * (-t153 - t186) - t29 * t193) * t188 + (t28 * (-rSges(5,3) + t193) + t29 * (-t186 - t339)) * t187) * qJD(1) - (qJD(1) * t251 - t119 + t215 - t28) * t29) * m(5) + (-(qJD(1) * t268 - t119 - t272 - t41) * t42 + t31 * (t187 * t269 + t184 + t285 - t353) + t30 * t368 + t42 * (t169 + t291) + (t151 * t336 - t337) * qJD(3) + ((-t195 * t42 - t197 * t41) * pkin(1) + (-pkin(2) - t243) * t335 + (t41 * (-rSges(4,3) - pkin(5)) + t42 * t269) * t187) * qJD(1)) * m(4) + (((t188 * t377 - t376 + t381) * t188 + (t187 * t377 + t382 + t409) * t187) * qJD(3) + t389 - t390) * t265 + (t385 + t386) * t280 / 0.2e1 - (t384 - t387 + t388) * t279 / 0.2e1 + ((t380 + t406) * t187 + (t379 + t405) * t188) * qJD(3) * t350; m(4) * t16 + m(5) * t1; -((((-t323 - t325) * t188 + (t322 + t324) * t187) * t196 + ((-t319 - t321) * t188 + (t318 + t320) * t187) * t194) * qJD(3) + ((t286 + t288) * t196 + (t287 + t289) * t194) * qJD(1)) * qJD(1) / 0.2e1 + (t387 * t188 + t386 * t187 + (t187 * t380 + t188 * t379) * qJD(1)) * t350 + ((-t280 * t300 + t283) * t187 + (t211 + (-t357 * t188 + (t301 + t213) * t187) * qJD(3)) * t188 + (-t280 * t302 + t284) * t187 + (t212 + (-t356 * t188 + (t303 + t214) * t187) * qJD(3)) * t188) * t265 + ((-t279 * t301 - t283) * t188 + (t211 + (t213 * t187 + (t300 - t357) * t188) * qJD(3)) * t187 + (-t279 * t303 - t284) * t188 + (t212 + (t214 * t187 + (t302 - t356) * t188) * qJD(3)) * t187) * t262 + (t16 * t241 + t40 * t210 + t242 * t130 + (-t30 * t187 - t31 * t188 + (-t188 * t42 + t336) * qJD(1)) * t151 - (t114 * t41 - t337) * qJD(1) - (t40 * (-t114 * t187 - t116 * t188) + t242 * t243) * qJD(3)) * m(4) + (t385 * qJD(1) + ((t381 * qJD(1) + t373 * t188) * t188 + (t371 * t187 + t382 * qJD(1) + (-t372 + t374) * t188) * t187) * t370) * t355 + (t384 * qJD(1) + ((t383 * qJD(1) + t372 * t188) * t188 + (t374 * t187 + t407 * qJD(1) + (-t371 + t373) * t188) * t187) * t370) * t354 + ((-t14 * t260 + t29 * t246 + t1 * t342 + t25 * t347 + (-t25 * t341 + t326) * qJD(1)) * t187 + (-t15 * t260 + t28 * t246 + t1 * t341 + t25 * t348 + (t25 * t342 - t260 * t29) * qJD(1)) * t188 - (t28 * t113 + t249 * t29) * qJD(1) - ((t249 * t25 + t259 * t28) * t188 + (t29 * t259 + (-t113 - t277) * t25) * t187) * qJD(3)) * m(5) + (t389 + t391) * t282 / 0.2e1 + (t388 + t392) * t281 / 0.2e1; 0.2e1 * (t14 * t354 + t15 * t355) * m(5);];
tauc = t2(:);
