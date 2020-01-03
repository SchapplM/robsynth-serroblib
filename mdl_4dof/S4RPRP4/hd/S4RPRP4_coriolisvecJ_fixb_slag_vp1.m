% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:50
% DurationCPUTime: 10.69s
% Computational Cost: add. (5529->379), mult. (6979->485), div. (0->0), fcn. (5407->6), ass. (0->222)
t195 = cos(qJ(3));
t193 = sin(qJ(3));
t317 = Icges(4,4) * t193;
t145 = Icges(4,2) * t195 + t317;
t189 = Icges(5,5) * t193;
t231 = Icges(5,3) * t195 - t189;
t408 = t145 + t231;
t314 = Icges(5,5) * t195;
t147 = Icges(5,1) * t193 - t314;
t190 = Icges(4,4) * t195;
t149 = Icges(4,1) * t193 + t190;
t407 = t147 + t149;
t142 = Icges(4,5) * t195 - Icges(4,6) * t193;
t192 = qJ(1) + pkin(6);
t187 = sin(t192);
t188 = cos(t192);
t85 = Icges(4,3) * t187 + t142 * t188;
t144 = Icges(5,4) * t195 + Icges(5,6) * t193;
t87 = Icges(5,2) * t187 + t144 * t188;
t413 = t85 + t87;
t233 = Icges(5,1) * t195 + t189;
t91 = Icges(5,4) * t187 + t188 * t233;
t150 = Icges(4,1) * t195 - t317;
t93 = Icges(4,5) * t187 + t150 * t188;
t401 = t91 + t93;
t406 = -t408 * t193 + t407 * t195;
t297 = t188 * t195;
t298 = t188 * t193;
t163 = Icges(5,5) * t297;
t310 = Icges(5,6) * t187;
t83 = Icges(5,3) * t298 + t163 + t310;
t412 = t413 * t187 + t401 * t297 + t83 * t298;
t140 = Icges(5,3) * t193 + t314;
t232 = -Icges(4,2) * t193 + t190;
t411 = (-t140 + t232) * qJD(3);
t410 = (t150 + t233) * qJD(3);
t141 = Icges(4,5) * t193 + Icges(4,6) * t195;
t143 = Icges(5,4) * t193 - Icges(5,6) * t195;
t409 = t141 + t143;
t299 = t187 * t195;
t300 = t187 * t193;
t311 = Icges(4,6) * t188;
t88 = Icges(4,4) * t299 - Icges(4,2) * t300 - t311;
t332 = t193 * t88;
t164 = Icges(4,4) * t300;
t315 = Icges(4,5) * t188;
t92 = Icges(4,1) * t299 - t164 - t315;
t237 = -t195 * t92 + t332;
t86 = -Icges(5,2) * t188 + t144 * t187;
t333 = t188 * t86;
t82 = -Icges(5,6) * t188 + t140 * t187;
t90 = -Icges(5,4) * t188 + t187 * t233;
t239 = t193 * t82 + t195 * t90;
t365 = t187 * t239;
t30 = -t333 + t365;
t308 = Icges(4,3) * t188;
t84 = Icges(4,5) * t299 - Icges(4,6) * t300 - t308;
t405 = -t187 * t237 - t188 * t84 + t30;
t89 = Icges(4,6) * t187 + t188 * t232;
t383 = -t298 * t89 + t412;
t301 = t143 * t188;
t303 = t141 * t188;
t404 = t406 * t187 - t301 - t303;
t302 = t143 * t187;
t304 = t141 * t187;
t403 = t406 * t188 + t302 + t304;
t249 = t188 * t87 - t91 * t299 - t83 * t300;
t71 = t93 * t299;
t255 = t188 * t85 - t71;
t33 = -t300 * t89 - t255;
t402 = -t249 + t33;
t400 = t410 * t195 - t411 * t193 + (-t407 * t193 - t408 * t195) * qJD(3) + t409 * qJD(1);
t399 = t408 * qJD(3);
t398 = t407 * qJD(3);
t397 = t237 - t239;
t331 = t193 * t89;
t396 = t193 * t83 + t401 * t195 - t331;
t395 = (-t142 - t144) * qJD(3) + t406 * qJD(1);
t394 = t403 * qJD(1);
t342 = -t187 * t84 - t92 * t297;
t36 = -t298 * t88 - t342;
t80 = t187 * t86;
t34 = t90 * t297 + t82 * t298 + t80;
t377 = t188 * t34;
t393 = (t383 * t187 - t188 * t36 - t377) * qJD(3);
t392 = (t402 * t187 - t405 * t188) * qJD(3);
t391 = t404 * qJD(1);
t390 = t391 + t392;
t389 = t393 + t394;
t388 = (t399 * t187 + (t140 * t188 + t310 - t89) * qJD(1)) * t195 + (-t401 * qJD(1) + t398 * t187) * t193 + t397 * qJD(3);
t387 = (-t399 * t188 + (-t187 * t232 + t311 + t82) * qJD(1)) * t195 + (-t398 * t188 + (-t150 * t187 + t315 - t90) * qJD(1)) * t193 + t396 * qJD(3);
t386 = -t395 * t187 + t400 * t188;
t385 = t400 * t187 + t395 * t188;
t384 = t34 + t36;
t382 = (-t82 + t88) * t195 + (t90 + t92) * t193;
t381 = (-t83 + t89) * t195 + t401 * t193;
t380 = t409 * qJD(3);
t379 = t333 + t412;
t376 = rSges(5,3) + qJ(4);
t378 = (-t380 * t187 + (t397 + t413) * qJD(1)) * t188;
t278 = qJD(1) * t188;
t197 = qJD(1) ^ 2;
t348 = rSges(5,1) + pkin(3);
t126 = t188 * pkin(2) + t187 * pkin(5);
t196 = cos(qJ(1));
t191 = t196 * pkin(1);
t363 = t191 + t126;
t179 = t187 * rSges(4,3);
t97 = rSges(4,1) * t297 - rSges(4,2) * t298 + t179;
t375 = t363 + t97;
t152 = pkin(3) * t193 - t195 * qJ(4);
t153 = rSges(5,1) * t193 - t195 * rSges(5,3);
t284 = t152 + t153;
t156 = pkin(3) * t195 + qJ(4) * t193;
t157 = rSges(5,1) * t195 + rSges(5,3) * t193;
t374 = t156 + t157;
t370 = -t380 * t188 + (-t142 * t187 + t308 - t396 - t86) * qJD(1);
t180 = t187 * rSges(5,2);
t320 = t297 * t348 + t298 * t376 + t180;
t369 = t363 + t320;
t368 = 0.2e1 * qJD(3);
t274 = qJD(4) * t195;
t182 = t188 * rSges(5,2);
t321 = t187 * t374 - t182;
t25 = -t274 + qJD(2) + (t187 * t321 + t188 * t320) * qJD(3);
t366 = qJD(3) * t25;
t251 = t284 * qJD(3);
t275 = qJD(4) * t193;
t215 = -t251 + t275;
t364 = t215 * t187;
t252 = t188 * rSges(3,1) - rSges(3,2) * t187;
t362 = t191 + t252;
t151 = t188 * t275;
t276 = qJD(3) * t188;
t266 = t195 * t276;
t361 = rSges(5,2) * t278 + t266 * t376 + t151;
t323 = -t149 * t187 - t88;
t327 = -Icges(4,2) * t299 - t164 + t92;
t353 = -t193 * t327 + t195 * t323;
t325 = -t147 * t187 + t82;
t329 = t231 * t187 - t90;
t352 = t193 * t329 + t195 * t325;
t113 = t153 * t187;
t343 = t156 * t278 + (-qJD(3) * t152 + t275) * t187 - qJD(3) * t113 + (t157 * t188 + t180) * qJD(1);
t267 = t193 * t276;
t279 = qJD(1) * t187;
t213 = -t195 * t279 - t267;
t270 = t193 * t279;
t344 = t213 * t348 - t270 * t376 + t361;
t1 = (t275 + t344 * t188 + t343 * t187 + (-t187 * t320 + t188 * t321) * qJD(1)) * qJD(3);
t351 = m(5) * t1;
t194 = sin(qJ(1));
t347 = pkin(1) * t194;
t345 = qJD(1) / 0.2e1;
t340 = rSges(4,1) * t195;
t154 = rSges(4,1) * t193 + rSges(4,2) * t195;
t118 = t154 * t188;
t277 = qJD(3) * t187;
t269 = t154 * t277;
t40 = qJD(1) * t375 - t269;
t338 = t118 * t40;
t185 = t188 * pkin(5);
t125 = pkin(2) * t187 - t185;
t282 = rSges(4,2) * t300 + t188 * rSges(4,3);
t95 = rSges(4,1) * t299 - t282;
t264 = -t95 - t347;
t268 = t154 * t276;
t39 = -t268 + (-t125 + t264) * qJD(1);
t337 = t187 * t39;
t228 = -t188 * t251 + t151;
t250 = -t321 - t347;
t28 = (-t125 + t250) * qJD(1) + t228;
t336 = t188 * t28;
t335 = t188 * t39;
t328 = t231 * t188 - t91;
t326 = -t145 * t188 + t93;
t324 = -Icges(5,1) * t298 + t163 + t83;
t322 = -t149 * t188 - t89;
t294 = -t152 * t187 - t113;
t293 = t284 * t188;
t292 = -qJD(3) * t374 + t274;
t289 = rSges(4,2) * t270 + rSges(4,3) * t278;
t288 = -t231 + t233;
t287 = t140 - t147;
t286 = -t145 + t150;
t285 = t149 + t232;
t281 = qJD(1) * t142;
t280 = qJD(1) * t144;
t273 = t197 * t347;
t272 = t197 * t191;
t265 = -pkin(2) - t340;
t261 = -t277 / 0.2e1;
t258 = t276 / 0.2e1;
t256 = t185 - t347;
t254 = -t84 + t331;
t172 = pkin(5) * t278;
t248 = qJD(1) * (-pkin(2) * t279 + t172) - t273;
t245 = -t195 * t348 - pkin(2);
t124 = rSges(3,1) * t187 + rSges(3,2) * t188;
t242 = -rSges(4,2) * t193 + t340;
t241 = -t187 * t40 - t335;
t240 = t187 * t95 + t188 * t97;
t227 = qJD(3) * (t274 + t292);
t114 = t154 * t187;
t212 = t193 * t328 + t195 * t324;
t211 = -t193 * t326 + t195 * t322;
t210 = (t193 * t287 + t195 * t288) * qJD(1);
t209 = (-t193 * t285 + t195 * t286) * qJD(1);
t66 = rSges(4,1) * t213 - rSges(4,2) * t266 + t289;
t68 = -qJD(3) * t114 + (t188 * t242 + t179) * qJD(1);
t208 = t187 * t68 + t188 * t66 + (-t187 * t97 + t188 * t95) * qJD(1);
t134 = t242 * qJD(3);
t123 = qJD(1) * t125;
t122 = t126 * qJD(1);
t38 = qJD(3) * t240 + qJD(2);
t29 = qJD(1) * t369 + t364;
t27 = -t272 - t134 * t276 + (-t122 - t68 + t269) * qJD(1);
t26 = -t134 * t277 + (t66 - t268) * qJD(1) + t248;
t16 = t208 * qJD(3);
t15 = -t272 + t188 * t227 + (-t122 - t343 - t364) * qJD(1);
t14 = t187 * t227 + (t188 * t215 + t344) * qJD(1) + t248;
t2 = [m(3) * ((-t124 * t197 - t273) * t362 + (-t272 + (-0.2e1 * t252 - t191 + t362) * t197) * (-t124 - t347)) + (((t33 - t71 + (t85 + t332) * t188 + t342) * t188 - t377 + (t30 - t365 + t379) * t187) * qJD(3) + t394) * t258 + (t406 * qJD(3) + t410 * t193 + t411 * t195) * qJD(1) + (t15 * (t182 + t256) + t14 * t369 + t29 * (-t348 * t267 + t172 + t361) + ((-t194 * t29 - t196 * t28) * pkin(1) + (-t193 * t376 + t245) * t336) * qJD(1) + (t15 * t245 + (-t28 * qJD(4) - t15 * t376) * t193 + t28 * (t193 * t348 - t195 * t376) * qJD(3) + (t28 * (-rSges(5,2) - pkin(5)) + t29 * (-pkin(2) - t374)) * qJD(1)) * t187 - (qJD(1) * t250 - t123 + t228 - t28) * t29) * m(5) + (t27 * (t187 * t265 + t256 + t282) + t26 * t375 + t40 * (t172 + t289) + (t154 * t337 - t338) * qJD(3) + ((-t194 * t40 - t196 * t39) * pkin(1) + (-pkin(2) - t242) * t335 + (t39 * (-rSges(4,3) - pkin(5)) + t40 * t265) * t187) * qJD(1) - (qJD(1) * t264 - t123 - t268 - t39) * t40) * m(4) + (((t188 * t254 - t379 + t383) * t188 + (t187 * t254 + t249 + t255 + t384 - t80) * t187) * qJD(3) + t390 - t391) * t261 + (t386 + t387) * t277 / 0.2e1 - (t385 - t388 + t389) * t276 / 0.2e1 + ((t382 + t404) * t187 + (t381 + t403) * t188) * qJD(3) * t345; m(4) * t16 + t351; -((((-t327 + t329) * t188 + (t326 - t328) * t187) * t195 + ((-t323 - t325) * t188 + (t322 + t324) * t187) * t193) * qJD(3) + ((t285 - t287) * t195 + (t286 + t288) * t193) * qJD(1)) * qJD(1) / 0.2e1 + (t388 * t188 + t387 * t187 + (t187 * t382 + t188 * t381) * qJD(1)) * t345 + ((-t277 * t303 + t281) * t187 + (t209 + (-t353 * t188 + (t304 + t211) * t187) * qJD(3)) * t188 + (-t277 * t301 + t280) * t187 + (t210 + (-t352 * t188 + (t302 + t212) * t187) * qJD(3)) * t188) * t261 + ((-t276 * t304 - t281) * t188 + (t209 + (t211 * t187 + (t303 - t353) * t188) * qJD(3)) * t187 + (-t276 * t302 - t280) * t188 + (t210 + (t212 * t187 + (t301 - t352) * t188) * qJD(3)) * t187) * t258 + (-(t193 * t25 + (t187 * t29 + t336) * t195) * qJD(4) - (-t28 * t294 - t29 * t293) * qJD(1) - ((-t25 * t293 - t28 * t374) * t188 + (t25 * t294 - t29 * t374) * t187) * qJD(3) + (-t15 * t284 + t28 * t292 + t1 * t320 + t25 * t344 + (t25 * t321 - t284 * t29) * qJD(1)) * t188 + (-t14 * t284 + t29 * t292 + t1 * t321 + t25 * t343 + (-t25 * t320 + t28 * t284) * qJD(1)) * t187) * m(5) + (t16 * t240 + t38 * t208 + t241 * t134 + (-t26 * t187 - t27 * t188 + (-t188 * t40 + t337) * qJD(1)) * t154 - (t114 * t39 - t338) * qJD(1) - (t38 * (-t114 * t187 - t118 * t188) + t241 * t242) * qJD(3)) * m(4) + (t386 * qJD(1) + (t383 * t278 + (qJD(1) * t384 + t370 * t187 - t378) * t187) * t368) * t187 / 0.2e1 - (t385 * qJD(1) + ((t402 * qJD(1) + t378) * t188 + (t405 * qJD(1) - t370 * t188) * t187) * t368) * t188 / 0.2e1 + (t390 + t392) * t279 / 0.2e1 + (t389 + t393) * t278 / 0.2e1; -t195 * t351 + 0.2e1 * (m(5) * (t14 * t187 + t15 * t188 + t366) / 0.2e1 - m(5) * (t187 ^ 2 + t188 ^ 2) * t366 / 0.2e1) * t193;];
tauc = t2(:);
