% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR8
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:14
% DurationCPUTime: 6.39s
% Computational Cost: add. (5588->478), mult. (8986->653), div. (0->0), fcn. (6976->6), ass. (0->281)
t227 = qJ(3) + qJ(4);
t210 = sin(t227);
t211 = cos(t227);
t229 = sin(qJ(1));
t352 = t211 * t229;
t186 = Icges(5,4) * t352;
t353 = t210 * t229;
t231 = cos(qJ(1));
t356 = Icges(5,5) * t231;
t109 = Icges(5,1) * t353 + t186 + t356;
t226 = qJD(3) + qJD(4);
t171 = t226 * t229;
t172 = t226 * t231;
t359 = Icges(5,4) * t210;
t156 = Icges(5,1) * t211 - t359;
t269 = Icges(5,2) * t211 + t359;
t416 = -t269 + t156;
t358 = Icges(5,4) * t211;
t271 = Icges(5,1) * t210 + t358;
t110 = -Icges(5,5) * t229 + t231 * t271;
t154 = -Icges(5,2) * t210 + t358;
t418 = t154 * t231 + t110;
t395 = qJD(1) * t416 - t171 * t418 + t172 * (-Icges(5,2) * t353 + t109 + t186);
t417 = t154 + t271;
t108 = -Icges(5,6) * t229 + t231 * t269;
t419 = -t156 * t231 + t108;
t107 = Icges(5,6) * t231 + t229 * t269;
t420 = -t156 * t229 + t107;
t396 = qJD(1) * t417 - t171 * t419 + t172 * t420;
t422 = t210 * t396 - t211 * t395;
t230 = cos(qJ(3));
t421 = pkin(3) * t230;
t274 = rSges(5,1) * t210 + rSges(5,2) * t211;
t228 = sin(qJ(3));
t349 = t229 * t230;
t199 = Icges(4,4) * t349;
t351 = t228 * t229;
t357 = Icges(4,5) * t231;
t125 = Icges(4,1) * t351 + t199 + t357;
t360 = Icges(4,4) * t230;
t272 = Icges(4,1) * t228 + t360;
t126 = -Icges(4,5) * t229 + t231 * t272;
t176 = -Icges(4,2) * t228 + t360;
t145 = t176 * t231;
t243 = t229 * (t126 + t145) - t231 * (-Icges(4,2) * t351 + t125 + t199);
t361 = Icges(4,4) * t228;
t270 = Icges(4,2) * t230 + t361;
t123 = Icges(4,6) * t231 + t229 * t270;
t124 = -Icges(4,6) * t229 + t231 * t270;
t178 = Icges(4,1) * t230 - t361;
t146 = t178 * t229;
t147 = t178 * t231;
t244 = t229 * (t124 - t147) - t231 * (t123 - t146);
t415 = -t244 * t228 + t243 * t230;
t332 = t176 + t272;
t333 = -t270 + t178;
t414 = (t228 * t332 - t230 * t333) * qJD(1);
t413 = 2 * qJD(3);
t412 = rSges(4,2) * t230;
t372 = rSges(5,2) * t210;
t374 = rSges(5,1) * t211;
t158 = -t372 + t374;
t136 = t158 * t229;
t184 = t231 * pkin(1) + t229 * qJ(2);
t292 = pkin(5) * t231 + t184;
t316 = qJD(3) * t231;
t194 = t316 * t421;
t209 = qJD(2) * t231;
t330 = t194 + t209;
t223 = t231 * rSges(5,3);
t111 = rSges(5,1) * t353 + rSges(5,2) * t352 + t223;
t202 = pkin(3) * t351;
t232 = -pkin(6) - pkin(5);
t377 = pkin(5) + t232;
t139 = -t377 * t231 + t202;
t342 = -t111 - t139;
t42 = -t158 * t172 + (t292 - t342) * qJD(1) - t330;
t411 = t42 * (qJD(1) * t136 + t172 * t274);
t267 = Icges(5,5) * t210 + Icges(5,6) * t211;
t106 = -Icges(5,3) * t229 + t231 * t267;
t44 = -t231 * t106 - t108 * t352 - t110 * t353;
t261 = t154 * t211 + t156 * t210;
t152 = Icges(5,5) * t211 - Icges(5,6) * t210;
t355 = t152 * t231;
t54 = t229 * t261 + t355;
t410 = t54 * qJD(1) + t171 * t44;
t137 = t158 * t231;
t67 = -t226 * t137 + (t229 * t274 + t223) * qJD(1);
t409 = t136 * t171 + t172 * t137 + t231 * t67;
t265 = t108 * t211 + t110 * t210;
t408 = t231 * t265;
t262 = t124 * t230 + t126 * t228;
t405 = t262 * t231;
t290 = -rSges(3,2) * t231 + t229 * rSges(3,3);
t404 = t184 + t290;
t268 = Icges(4,5) * t228 + Icges(4,6) * t230;
t122 = -Icges(4,3) * t229 + t231 * t268;
t323 = qJD(1) * t122;
t72 = t124 * t228 - t126 * t230;
t77 = qJD(1) * t123 - qJD(3) * t145;
t79 = -qJD(3) * t147 + (t229 * t272 + t357) * qJD(1);
t403 = qJD(3) * t72 + t228 * t79 + t230 * t77 + t323;
t163 = t270 * qJD(3);
t164 = t272 * qJD(3);
t174 = Icges(4,5) * t230 - Icges(4,6) * t228;
t402 = qJD(1) * t174 + qJD(3) * (t176 * t228 - t178 * t230) + t163 * t230 + t164 * t228;
t263 = t123 * t228 - t125 * t230;
t121 = Icges(4,3) * t231 + t229 * t268;
t324 = qJD(1) * t121;
t318 = qJD(3) * t229;
t78 = qJD(1) * t124 + t176 * t318;
t80 = qJD(1) * t126 + qJD(3) * t146;
t401 = qJD(3) * t263 - t228 * t80 - t230 * t78 + t324;
t283 = t417 * t226;
t284 = t416 * t226;
t399 = qJD(1) * t152 + t210 * t283 - t211 * t284;
t286 = -qJD(1) * t107 + t226 * t418;
t288 = (t229 * t271 + t356) * qJD(1) + t419 * t226;
t325 = qJD(1) * t106;
t398 = t210 * t288 - t211 * t286 + t325;
t287 = qJD(1) * t108 + t109 * t226 + t154 * t171;
t289 = -qJD(1) * t110 + t226 * t420;
t105 = Icges(5,3) * t231 + t229 * t267;
t326 = qJD(1) * t105;
t397 = t210 * t289 - t211 * t287 + t326;
t394 = t229 ^ 2;
t321 = qJD(1) * t229;
t160 = t226 * t321;
t393 = -t160 / 0.2e1;
t161 = qJD(1) * t172;
t392 = t161 / 0.2e1;
t391 = -t171 / 0.2e1;
t390 = t171 / 0.2e1;
t389 = -t172 / 0.2e1;
t388 = t172 / 0.2e1;
t387 = t229 / 0.2e1;
t386 = t231 / 0.2e1;
t385 = rSges(3,2) - pkin(1);
t384 = -rSges(5,3) - pkin(1);
t383 = pkin(3) * t228;
t381 = pkin(5) * t229;
t380 = pkin(5) * qJD(1) ^ 2;
t379 = -qJD(1) / 0.2e1;
t378 = qJD(1) / 0.2e1;
t320 = qJD(1) * t231;
t309 = t171 * t374 + t274 * t320;
t311 = t226 * t372;
t68 = (-rSges(5,3) * qJD(1) - t311) * t229 + t309;
t317 = qJD(3) * t230;
t302 = t229 * t317;
t192 = pkin(3) * t302;
t306 = t228 * t320;
t307 = pkin(3) * t306 + t232 * t321 + t192;
t91 = pkin(5) * t321 + t307;
t376 = -t68 - t91;
t370 = rSges(3,3) * t231;
t224 = t231 * rSges(4,3);
t213 = t231 * qJ(2);
t180 = pkin(1) * t229 - t213;
t208 = qJD(2) * t229;
t277 = t171 * t158 + t192 + t208;
t220 = t229 * rSges(5,3);
t112 = t231 * t274 - t220;
t312 = t231 * t383;
t138 = t229 * t377 + t312;
t341 = t112 + t138;
t279 = t341 - t381;
t41 = (-t180 + t279) * qJD(1) + t277;
t368 = t231 * t41;
t183 = rSges(4,1) * t230 - rSges(4,2) * t228;
t150 = t183 * t318;
t221 = t229 * rSges(4,3);
t275 = rSges(4,1) * t228 + t412;
t128 = t275 * t231 - t221;
t293 = t128 - t381;
t58 = t150 + t208 + (-t180 + t293) * qJD(1);
t367 = t231 * t58;
t120 = t274 * t226;
t233 = qJD(3) ^ 2;
t314 = qJD(1) * qJD(2);
t329 = qJ(2) * t320 + t208;
t336 = qJD(1) * (-pkin(1) * t321 + t329) + t229 * t314;
t257 = -t229 * t380 + t336;
t28 = t233 * t312 + t120 * t172 + t158 * t160 + (t192 - t376) * qJD(1) + t257;
t366 = t28 * t231;
t140 = qJD(1) * t184 - t209;
t204 = t231 * t314;
t280 = -t231 * t380 + t204;
t92 = qJD(1) * t139 - t194;
t29 = -t233 * t202 - t120 * t171 + t158 * t161 + (-t140 - t67 - t92 + t194) * qJD(1) + t280;
t365 = t29 * t229;
t166 = t275 * qJD(3);
t308 = t320 * t412 + (t302 + t306) * rSges(4,1);
t319 = qJD(3) * t228;
t84 = (-rSges(4,2) * t319 - rSges(4,3) * qJD(1)) * t229 + t308;
t37 = t166 * t316 + (t84 + t150) * qJD(1) + t257;
t364 = t37 * t231;
t304 = t183 * t316;
t149 = t183 * t231;
t83 = -qJD(3) * t149 + (t229 * t275 + t224) * qJD(1);
t38 = -t166 * t318 + (-t140 - t83 + t304) * qJD(1) + t280;
t363 = t38 * t229;
t362 = t41 * t158;
t354 = t174 * t231;
t350 = t229 * t106;
t130 = t229 * t152;
t142 = t229 * t174;
t55 = t231 * t261 - t130;
t348 = t55 * qJD(1);
t260 = t230 * t176 + t228 * t178;
t74 = t231 * t260 - t142;
t347 = t74 * qJD(1);
t328 = rSges(3,2) * t321 + rSges(3,3) * t320;
t168 = qJD(1) * t180;
t327 = t208 - t168;
t315 = t268 * qJD(1);
t313 = -rSges(4,3) - pkin(1) - pkin(5);
t43 = t231 * t105 + t107 * t352 + t109 * t353;
t47 = t231 * t121 + t123 * t349 + t125 * t351;
t48 = -t231 * t122 - t124 * t349 - t126 * t351;
t127 = rSges(4,1) * t351 + rSges(4,2) * t349 + t224;
t303 = t228 * t318;
t300 = -t321 / 0.2e1;
t299 = t320 / 0.2e1;
t298 = -t318 / 0.2e1;
t296 = -t316 / 0.2e1;
t294 = t158 + t421;
t281 = qJD(1) * t137 - t171 * t274;
t59 = -t304 - t209 + (t127 + t292) * qJD(1);
t273 = t229 * t58 - t231 * t59;
t266 = t107 * t211 + t109 * t210;
t57 = t108 * t210 - t110 * t211;
t264 = t123 * t230 + t125 * t228;
t258 = t43 + t350;
t255 = (t229 * t48 + t231 * t47) * qJD(3);
t113 = t229 * t121;
t49 = -t264 * t231 + t113;
t50 = -t229 * t122 + t405;
t254 = (t229 * t50 + t231 * t49) * qJD(3);
t69 = (-t127 * t229 - t128 * t231) * qJD(3);
t253 = t274 + t383;
t250 = -qJD(1) * t267 + t130 * t172 - t171 * t355;
t96 = t229 * t105;
t45 = -t231 * t266 + t96;
t248 = -qJD(1) * t265 - t226 * t355 + t326;
t247 = qJD(1) * t266 + t130 * t226 + t325;
t246 = -qJD(1) * t262 - qJD(3) * t354 + t324;
t245 = qJD(1) * t264 + qJD(3) * t142 + t323;
t242 = qJD(1) * t261 - t267 * t226;
t241 = t260 * qJD(1) - t268 * qJD(3);
t10 = -t229 * t397 + t247 * t231;
t11 = t229 * t398 + t248 * t231;
t20 = t172 * t43 + t410;
t46 = -t350 + t408;
t21 = t171 * t46 + t172 * t45 - t348;
t26 = t242 * t229 + t231 * t399;
t27 = -t229 * t399 + t242 * t231;
t30 = -t210 * t287 - t211 * t289;
t31 = t210 * t286 + t211 * t288;
t56 = -t107 * t210 + t109 * t211;
t8 = t247 * t229 + t231 * t397;
t9 = t248 * t229 - t231 * t398;
t240 = (qJD(1) * t26 - t160 * t45 + t161 * t46 + t171 * t9 + t172 * t8) * t387 + (-t210 * t395 - t211 * t396) * t379 + t20 * t300 + t21 * t299 + (qJD(1) * t27 + t10 * t172 + t11 * t171 - t160 * t43 + t161 * t44) * t386 + (t229 * t44 + t231 * t43) * t393 + (t229 * t46 + t231 * t45) * t392 + (t229 * t9 + t231 * t8 + (-t229 * t45 + t231 * t46) * qJD(1)) * t390 + (t10 * t231 + t11 * t229 + (-t229 * t43 + t231 * t44) * qJD(1)) * t388 + (t229 * t31 + t231 * t30 + (-t229 * t56 + t231 * t57) * qJD(1)) * t378 + (t250 * t229 + t422 * t231) * t391 + (-t422 * t229 + t250 * t231) * t389;
t181 = rSges(3,2) * t229 + t370;
t148 = t183 * t229;
t99 = t231 * t112;
t95 = qJD(1) * t404 - t209;
t94 = t208 + (-t180 + t181) * qJD(1);
t82 = t204 + (-qJD(1) * t290 - t140) * qJD(1);
t81 = qJD(1) * t328 + t336;
t73 = t229 * t260 + t354;
t70 = t73 * qJD(1);
t36 = -t111 * t171 - t112 * t172 + (-t138 * t231 - t139 * t229) * qJD(3);
t35 = -t229 * t402 + t241 * t231;
t34 = t241 * t229 + t231 * t402;
t33 = t262 * qJD(3) - t228 * t77 + t230 * t79;
t32 = -qJD(3) * t264 - t228 * t78 + t230 * t80;
t23 = t254 - t347;
t22 = t70 + t255;
t13 = -t111 * t161 + t112 * t160 - t171 * t68 + t172 * t67 + (-t229 * t91 + t231 * t92 + (t138 * t229 - t139 * t231) * qJD(1)) * qJD(3);
t1 = [(t70 + ((-t49 + t113 + t48) * t229 + (t50 - t405 + (t122 - t264) * t229 + t47) * t231) * qJD(3)) * t298 + ((t258 + t46 - t408) * t172 + t410) * t391 + t57 * t392 - t161 * t55 / 0.2e1 + (t54 + t56) * t393 + (t348 + (t265 * t229 + t44 - t96) * t172 + (t258 - t43) * t171 + ((t106 + t266) * t172 - t265 * t171) * t231 + t21) * t389 + (t30 + t27) * t388 + (t347 + (t394 * t122 + (-t113 + t48 + (t122 + t264) * t231) * t231) * qJD(3) + t23) * t296 + (-qJD(3) * t260 + t163 * t228 - t164 * t230 - t210 * t284 - t211 * t283) * qJD(1) + (t29 * (t229 * t232 - t180 - t220) + t41 * t330 + t28 * (t202 + t111 + t184) + t42 * (-t229 * t311 + t307 + t309 + t329) + (t226 * t362 - t28 * t232 + t253 * t29) * t231 + ((t232 + t384) * t368 + (t41 * (-qJ(2) - t253) + t42 * t384) * t229) * qJD(1) - (qJD(1) * t279 - t168 + t277 - t41) * t42) * m(5) + (t38 * (-t180 - t221 - t381) + t58 * t209 + t37 * (t127 + t184) + t59 * (-rSges(4,2) * t303 + t308 + t329) + (qJD(3) * t183 * t58 + t37 * pkin(5) + t275 * t38) * t231 + (t313 * t367 + (t58 * (-qJ(2) - t275) + t59 * t313) * t229) * qJD(1) - (qJD(1) * t293 + t150 + t327 - t58) * t59) * m(4) + (t82 * (t229 * t385 + t213 + t370) + t94 * t209 + t81 * t404 + t95 * (t328 + t329) + (t94 * t385 * t231 + (t94 * (-rSges(3,3) - qJ(2)) - t95 * pkin(1)) * t229) * qJD(1) - (qJD(1) * t181 + t327 - t94) * t95) * m(3) + (t31 + t26 + t20) * t390 + (t33 + t34 + t22) * t318 / 0.2e1 + (qJD(1) * t72 + t32 + t35) * t316 / 0.2e1 + (t231 * t74 + (-t263 + t73) * t229) * qJD(3) * t379; 0.2e1 * (t365 / 0.2e1 - t366 / 0.2e1) * m(5) + 0.2e1 * (t363 / 0.2e1 - t364 / 0.2e1) * m(4) + 0.2e1 * (t82 * t387 - t81 * t231 / 0.2e1) * m(3); t240 + (t33 * t229 + t32 * t231 + (t229 * t263 + t72 * t231) * qJD(1)) * t378 + ((-t318 * t354 - t315) * t229 + (t414 + (t229 * t142 + t415) * qJD(3)) * t231) * t298 + ((-t228 * t333 - t230 * t332) * qJD(1) + (t228 * t243 + t230 * t244) * qJD(3)) * t379 + ((t142 * t316 - t315) * t231 + (-t414 + (-t231 * t354 - t415) * qJD(3)) * t229) * t296 + (qJD(1) * t34 + ((t245 * t229 + t231 * t401) * t231 + (t246 * t229 - t231 * t403) * t229 + (-t49 * t229 + t50 * t231) * qJD(1)) * t413) * t387 + (qJD(1) * t35 + ((-t229 * t401 + t245 * t231) * t231 + (t229 * t403 + t246 * t231) * t229 + (-t47 * t229 + t48 * t231) * qJD(1)) * t413) * t386 + (t255 + t22) * t300 + (t254 + t23) * t299 + (-t41 * (-pkin(3) * t303 + t281) - t411 - t13 * t99 + (qJD(1) * t362 + t42 * t120 - t13 * t138 - t28 * t294) * t231 + (t29 * t294 + t41 * (-pkin(3) * t319 - t120) + t13 * t342 + t42 * qJD(1) * t158) * t229 + (-(-t231 ^ 2 - t394) * pkin(3) * t317 + (qJD(1) * t342 + t92) * t231 + (qJD(1) * t341 + t376) * t229 + t409) * t36) * m(5) + (0.2e1 * t69 * (-t229 * t84 + t231 * t83 + (-t127 * t231 + t128 * t229) * qJD(1)) - t273 * t166 + (t363 - t364 + (t229 * t59 + t367) * qJD(1)) * t183 - (t148 * t59 + t149 * t58) * qJD(1) - (t69 * (-t148 * t229 - t149 * t231) - t273 * t275) * qJD(3)) * m(4); t240 + (-t41 * t281 - t411 + t13 * (-t111 * t229 - t99) - (t229 * t41 - t231 * t42) * t120 + (t365 - t366 + (t229 * t42 + t368) * qJD(1)) * t158 + (-t229 * t68 + (-t111 * t231 + t112 * t229) * qJD(1) + t409) * t36) * m(5);];
tauc = t1(:);
