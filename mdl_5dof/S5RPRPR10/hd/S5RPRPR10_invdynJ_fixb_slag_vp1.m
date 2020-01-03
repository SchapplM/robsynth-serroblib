% Calculate vector of inverse dynamics joint torques for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR10_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:54
% EndTime: 2019-12-31 18:26:04
% DurationCPUTime: 7.01s
% Computational Cost: add. (10017->457), mult. (11853->541), div. (0->0), fcn. (11890->8), ass. (0->241)
t222 = qJDD(1) - qJDD(3);
t223 = qJD(1) - qJD(3);
t228 = qJD(1) ^ 2;
t389 = cos(qJ(1));
t216 = t389 * qJ(2);
t388 = sin(qJ(1));
t323 = t388 * pkin(1);
t180 = t323 - t216;
t214 = qJD(2) * t389;
t184 = t389 * pkin(1) + t388 * qJ(2);
t394 = qJD(1) * t184 - t214;
t247 = -qJDD(1) * t180 + qJDD(2) * t388 + (t214 - t394) * qJD(1);
t230 = (-qJDD(1) * t388 - t228 * t389) * pkin(2) + t247;
t225 = sin(qJ(3));
t312 = t388 * t225;
t227 = cos(qJ(3));
t313 = t389 * t227;
t249 = t312 + t313;
t311 = t388 * t227;
t314 = t389 * t225;
t419 = t311 - t314;
t340 = -rSges(4,1) * t419 + rSges(4,2) * t249;
t114 = t223 * t249;
t307 = qJD(1) * t388;
t308 = qJD(1) * t389;
t400 = t419 * qJD(3) + t225 * t308;
t115 = -t227 * t307 + t400;
t59 = t114 * rSges(4,1) - t115 * rSges(4,2);
t34 = t222 * t340 - t223 * t59 + t230;
t420 = t34 - g(1);
t327 = qJ(3) + pkin(8);
t300 = sin(t327);
t301 = cos(t327);
t150 = -t388 * t300 - t389 * t301;
t151 = t389 * t300 - t388 * t301;
t110 = -t151 * pkin(4) - t150 * pkin(7);
t224 = sin(qJ(5));
t347 = t151 * t224;
t319 = rSges(6,2) * t347;
t226 = cos(qJ(5));
t346 = t151 * t226;
t343 = rSges(6,1) * t346 + t150 * rSges(6,3);
t76 = t319 - t343;
t363 = t110 + t76;
t360 = Icges(6,4) * t226;
t275 = -Icges(6,2) * t224 + t360;
t66 = -Icges(6,6) * t151 + t150 * t275;
t367 = t224 * t66;
t361 = Icges(6,4) * t224;
t277 = Icges(6,1) * t226 - t361;
t70 = -Icges(6,5) * t151 + t150 * t277;
t278 = t226 * t70 - t367;
t67 = Icges(6,6) * t150 + t151 * t275;
t368 = t224 * t67;
t71 = Icges(6,5) * t150 + t151 * t277;
t418 = -t226 * t71 + t368;
t273 = Icges(6,5) * t226 - Icges(6,6) * t224;
t63 = Icges(6,3) * t150 + t151 * t273;
t382 = -t150 * t63 - t346 * t71;
t62 = -Icges(6,3) * t151 + t150 * t273;
t417 = -(t62 - t368) * t151 + t382;
t279 = -t224 * t70 - t226 * t66;
t281 = -t224 * t71 - t226 * t67;
t332 = qJD(5) * t151;
t333 = qJD(5) * t150;
t390 = -t223 / 0.2e1;
t416 = ((-t279 * t150 + t281 * t151) * qJD(5) + t279 * t333 - t281 * t332) * t390;
t116 = -rSges(4,1) * t249 - rSges(4,2) * t419;
t220 = t389 * pkin(2);
t213 = qJD(2) * t388;
t336 = qJ(2) * t308 + t213;
t256 = -pkin(1) * t307 + t336;
t253 = qJDD(1) * t184 - qJDD(2) * t389 + (t213 + t256) * qJD(1);
t322 = t388 * pkin(2);
t235 = qJDD(1) * t220 - t228 * t322 + t253;
t60 = t115 * rSges(4,1) + t114 * rSges(4,2);
t35 = -t116 * t222 + t223 * t60 + t235;
t415 = t35 - g(2);
t246 = t213 + (-t322 - t180) * qJD(1);
t209 = pkin(3) * t227 + pkin(2);
t188 = t388 * t209;
t299 = pkin(3) * t314;
t268 = t299 - t188;
t128 = t322 + t268;
t317 = t128 - t363;
t287 = rSges(6,1) * t224 + rSges(6,2) * t226;
t330 = qJD(5) * t287;
t23 = t150 * t330 + t223 * t317 + t246;
t92 = t287 * t151;
t414 = t23 * t92;
t413 = t151 * t62;
t410 = t63 * t151;
t111 = -t150 * pkin(4) + t151 * pkin(7);
t350 = t150 * t224;
t321 = rSges(6,2) * t350;
t349 = t150 * t226;
t342 = -rSges(6,1) * t349 + t151 * rSges(6,3);
t77 = t321 + t342;
t362 = t111 + t77;
t408 = t63 + t367;
t407 = t223 * t340;
t341 = t151 * rSges(5,1) - t150 * rSges(5,2);
t106 = t223 * t150;
t107 = t223 * t151;
t365 = t107 * rSges(5,1) - t106 * rSges(5,2);
t402 = t223 * t341 - t365;
t397 = t220 + t184;
t244 = qJD(1) * t397 - t214;
t272 = Icges(6,5) * t224 + Icges(6,6) * t226;
t87 = t272 * t150;
t86 = t272 * t151;
t334 = t389 * rSges(3,1) + t388 * rSges(3,3);
t398 = t184 + t334;
t201 = pkin(3) * t312;
t401 = pkin(3) * t313 + t201;
t329 = qJD(5) * t224;
t262 = -t106 * t226 + t151 * t329;
t260 = t107 * t226 + t150 * t329;
t328 = qJD(5) * t226;
t261 = -t107 * t224 + t150 * t328;
t399 = t261 * rSges(6,2);
t274 = Icges(6,2) * t226 + t361;
t276 = Icges(6,1) * t224 + t360;
t271 = -t224 * t274 + t226 * t276;
t47 = t151 * t271 + t87;
t43 = t47 * t223;
t48 = t150 * t271 - t86;
t44 = t48 * t223;
t239 = t268 - t180;
t396 = t239 - t110;
t296 = pkin(2) * t307;
t335 = -qJD(1) * t180 + t213;
t257 = -t296 + t335;
t395 = -t223 * t128 + t256 - t257;
t109 = -t150 * rSges(5,1) - t151 * rSges(5,2);
t376 = t274 * t151 - t71;
t378 = -t276 * t151 - t67;
t393 = t224 * t376 + t226 * t378;
t381 = t150 * t62 + t346 * t70;
t380 = t349 * t71 - t410;
t379 = t349 * t70 - t413;
t377 = -t276 * t150 - t66;
t375 = t274 * t150 - t70;
t337 = t389 * t209 + t201;
t129 = -t220 + t337;
t24 = t151 * t330 + (t129 + t362) * t223 + t244;
t366 = t24 * t287;
t364 = t106 * rSges(5,1) + t107 * rSges(5,2);
t345 = t273 * t223;
t344 = t341 + t128;
t339 = -t274 + t277;
t338 = -t275 - t276;
t183 = -rSges(6,1) * t226 + rSges(6,2) * t224;
t162 = t183 * qJD(5);
t331 = qJD(5) * t162;
t326 = -t389 / 0.2e1;
t325 = t388 / 0.2e1;
t294 = -qJD(1) * t201 + t401 * qJD(3) - t209 * t308;
t297 = pkin(2) * t308;
t97 = -t297 - t294;
t324 = -t97 + t364;
t316 = t388 * rSges(3,1);
t306 = t333 / 0.2e1;
t305 = -t332 / 0.2e1;
t304 = t332 / 0.2e1;
t303 = -t106 * pkin(4) + t107 * pkin(7);
t302 = -t107 * pkin(4) - t106 * pkin(7);
t295 = -t337 - t109;
t293 = t400 * pkin(3) - t209 * t307;
t28 = Icges(6,4) * t260 + Icges(6,2) * t261 + Icges(6,6) * t106;
t30 = Icges(6,1) * t260 + Icges(6,4) * t261 + Icges(6,5) * t106;
t237 = qJD(5) * t279 + t224 * t28 - t226 * t30;
t263 = t106 * t224 + t151 * t328;
t27 = Icges(6,4) * t262 + Icges(6,2) * t263 + Icges(6,6) * t107;
t29 = Icges(6,1) * t262 + Icges(6,4) * t263 + Icges(6,5) * t107;
t238 = qJD(5) * t281 + t224 * t27 - t226 * t29;
t25 = Icges(6,5) * t262 + Icges(6,6) * t263 + Icges(6,3) * t107;
t26 = Icges(6,5) * t260 + Icges(6,6) * t261 + Icges(6,3) * t106;
t289 = -(-t106 * t418 - t107 * t63 - t150 * t25 + t151 * t238) * t150 + t151 * (t106 * t278 - t107 * t62 - t150 * t26 + t151 * t237);
t288 = -t150 * (-t106 * t63 + t107 * t418 + t150 * t238 + t151 * t25) + t151 * (-t106 * t62 - t107 * t278 + t150 * t237 + t151 * t26);
t17 = -t347 * t67 - t382;
t18 = -t347 * t66 + t381;
t286 = -t150 * t17 + t151 * t18;
t19 = -t350 * t67 + t380;
t20 = -t350 * t66 + t379;
t285 = -t150 * t19 + t151 * t20;
t284 = -t150 * t23 - t151 * t24;
t258 = t262 * rSges(6,1) + t107 * rSges(6,3);
t31 = rSges(6,2) * t263 + t258;
t259 = t260 * rSges(6,1) + t106 * rSges(6,3);
t32 = t259 + t399;
t283 = -t150 * t32 - t151 * t31;
t282 = t150 * t77 + t151 * t76;
t121 = -t224 * t276 - t226 * t274;
t270 = t401 + t362;
t202 = pkin(3) * t311;
t269 = -t299 + t202;
t255 = -t323 - t322;
t254 = -t316 - t323;
t186 = rSges(2,1) * t389 - rSges(2,2) * t388;
t182 = rSges(2,1) * t388 + rSges(2,2) * t389;
t250 = -t111 - t337 - t342;
t245 = t224 * t375 + t226 * t377;
t242 = t296 + t293;
t241 = t110 - t299 - t343;
t240 = (t224 * t338 + t226 * t339) * t223;
t158 = t275 * qJD(5);
t159 = t277 * qJD(5);
t234 = qJD(5) * t121 - t158 * t224 + t159 * t226;
t233 = -t259 - t293 + t302;
t232 = t258 - t294 + t303;
t231 = t222 * t129 + t223 * t242 + t235;
t11 = -qJD(5) * t418 - t224 * t29 - t226 * t27;
t12 = qJD(5) * t278 - t224 * t30 - t226 * t28;
t157 = t273 * qJD(5);
t13 = t106 * t271 - t107 * t272 + t150 * t157 + t151 * t234;
t14 = -t106 * t272 - t107 * t271 + t150 * t234 - t151 * t157;
t6 = qJD(5) * t286 + t43;
t7 = qJD(5) * t285 + t44;
t79 = qJD(5) * t106 + qJDD(5) * t151;
t80 = qJD(5) * t107 - qJDD(5) * t150;
t229 = t223 * (-t159 * t224 + t274 * t329 + (-qJD(5) * t276 - t158) * t226) + t304 * t6 - (-t279 + t48) * t79 / 0.2e1 - (-t281 + t47) * t80 / 0.2e1 + (t12 + t14) * t305 + (t11 + t13 + t7) * t306 + (-Icges(4,3) - Icges(5,3) + t121) * t222;
t218 = t389 * rSges(3,3);
t211 = rSges(3,3) * t308;
t181 = t316 - t218;
t147 = t223 * t269;
t93 = t287 * t150;
t84 = -t223 * t116 + t244;
t83 = t246 + t407;
t58 = qJDD(1) * t334 + qJD(1) * (-rSges(3,1) * t307 + t211) + t253;
t57 = -qJDD(1) * t181 - t228 * t334 + t247;
t46 = (t109 + t129) * t223 + t244;
t45 = t223 * t344 + t246;
t33 = qJD(5) * t282 - qJD(4);
t16 = t222 * t109 + t223 * t365 + t231;
t15 = t222 * t344 + t223 * t324 + t230;
t10 = (-t150 * t67 - t151 * t66) * t224 + t380 + t381;
t9 = -t151 * t331 + t79 * t287 + (-t302 + t32) * t223 + t362 * t222 + t231;
t8 = -t150 * t331 - t80 * t287 + (-t303 - t31 - t97) * t223 + t317 * t222 + t230;
t5 = qJD(5) * t283 - t76 * t79 + t77 * t80 + qJDD(4);
t1 = [t416 - t229 - m(2) * (-g(1) * t182 + g(2) * t186) + (-t43 + ((t19 + (-t278 + t63) * t151) * t151 + (t408 * t150 + t20 - t379 - t413) * t150) * qJD(5)) * t305 + (t44 + ((t418 * t151 + t17 + t379) * t151 + (-t408 * t151 - t10 + t18) * t150) * qJD(5)) * t306 + (m(2) * (t182 ^ 2 + t186 ^ 2) + Icges(2,3) + Icges(3,2)) * qJDD(1) + (-g(1) * (-t76 + t396) - t366 * t333 + t8 * (t343 + t396) + t23 * (-pkin(1) * t308 - qJ(2) * t307 + t214 - t232) + (-t23 * t263 - t347 * t8) * rSges(6,2) + (-g(2) + t9) * (-t250 + t321 + t184) + (t363 * t223 + t23 - t233 + t395 + t399) * t24) * m(6) + (t45 * (t294 + t364 - t394) + (t16 - g(2)) * (-t295 + t184) + (t293 + t395 - t402 + t45) * t46 + (t15 - g(1)) * (t239 + t341)) * m(5) + (t415 * (-t116 + t397) + (-t244 - t59) * t83 + (qJD(1) * t255 - t257 + t336 - t407 + t60 + t83) * t84 + t420 * (t216 + t255 + t340)) * m(4) + ((t58 - g(2)) * t398 + (t57 - g(1)) * (t216 + t218 + t254) + (t211 - t335 + t336 + (t181 + t254) * qJD(1)) * (qJD(1) * t398 - t214)) * m(3); (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t388 - g(2) * t389) + 0.2e1 * (t325 * t8 + t326 * t9) * m(6) + 0.2e1 * (t15 * t325 + t16 * t326) * m(5) + 0.2e1 * (t325 * t34 + t326 * t35) * m(4) + 0.2e1 * (t325 * t57 + t326 * t58) * m(3); t416 + t229 + (t43 + ((t10 - t19) * t151 + (t278 * t150 - t20 + t417) * t150) * qJD(5)) * t305 + (-t44 + ((-t17 - t417) * t151 + (-t18 + (-t418 + t62) * t150 - t410) * t150) * qJD(5)) * t306 + (-t24 * t147 - (-t150 * t366 + t414) * qJD(5) - (t23 * t270 + t24 * t363) * t223 + t8 * (-t322 + t188 + t241) + t23 * (-t297 + t232) + t9 * (t220 + t250) + t24 * (-t296 + t233) + ((-t150 * t24 + t151 * t23) * t328 + (t106 * t23 + t107 * t24 - t150 * t9 + t151 * t8) * t224) * rSges(6,2) - g(1) * (t202 + t241 + t319) + g(2) * t270) * m(6) + (-t15 * t344 - t45 * t324 + t16 * (t220 + t295) - g(1) * (t269 - t341) + (-t147 - t242 + t402) * t46 + (-t223 * t45 + g(2)) * (t401 + t109)) * m(5) + (t59 * t83 - t60 * t84 - (-t223 * t84 + t420) * t340 + (t223 * t83 + t415) * t116) * m(4); (t5 + g(3)) * m(6) + (qJDD(4) + g(3)) * m(5); t106 * t7 / 0.2e1 + t151 * (qJD(5) * t288 + t14 * t223 + t19 * t80 + t20 * t79 + t222 * t48) / 0.2e1 + t79 * t285 / 0.2e1 + (t106 * t20 + t107 * t19 + t288) * t304 + t107 * t6 / 0.2e1 - t150 * (qJD(5) * t289 + t13 * t223 + t17 * t80 + t18 * t79 + t222 * t47) / 0.2e1 + t80 * t286 / 0.2e1 - (t106 * t18 + t107 * t17 + t289) * t333 / 0.2e1 + t222 * (t150 * t281 - t151 * t279) / 0.2e1 + t223 * (-t106 * t279 - t107 * t281 - t11 * t150 + t12 * t151) / 0.2e1 + ((t332 * t87 - t345) * t151 + (t240 + (-t393 * t150 + (-t86 + t245) * t151) * qJD(5)) * t150) * t305 + ((t333 * t86 + t345) * t150 + (t240 + (t245 * t151 + (-t393 - t87) * t150) * qJD(5)) * t151) * t306 + ((t224 * t339 - t226 * t338) * t223 + ((t150 * t376 - t151 * t375) * t226 + (-t150 * t378 + t151 * t377) * t224) * qJD(5)) * t390 + (-t5 * t282 + t33 * (t106 * t76 - t107 * t77 - t283) + t284 * t162 - (-t24 * t106 + t23 * t107 - t8 * t150 - t151 * t9) * t287 - (t24 * t93 - t414) * t223 - (t33 * (t150 * t93 + t151 * t92) + t284 * t183) * qJD(5) - g(1) * t93 - g(2) * t92 - g(3) * t183) * m(6);];
tau = t1;
