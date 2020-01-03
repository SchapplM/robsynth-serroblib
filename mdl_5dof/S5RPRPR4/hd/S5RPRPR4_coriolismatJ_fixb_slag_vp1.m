% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:17
% EndTime: 2020-01-03 11:38:38
% DurationCPUTime: 6.36s
% Computational Cost: add. (22019->405), mult. (14381->520), div. (0->0), fcn. (13192->10), ass. (0->251)
t270 = qJ(1) + pkin(8);
t261 = sin(t270);
t269 = qJ(3) + pkin(9);
t264 = qJ(5) + t269;
t255 = sin(t264);
t256 = cos(t264);
t212 = rSges(6,1) * t255 + rSges(6,2) * t256;
t260 = sin(t269);
t272 = sin(qJ(3));
t387 = pkin(3) * t272;
t232 = -pkin(4) * t260 - t387;
t285 = t212 - t232;
t122 = t285 * t261;
t412 = m(6) / 0.2e1;
t263 = cos(t270);
t124 = t285 * t263;
t441 = t124 * t263;
t413 = m(5) / 0.2e1;
t262 = cos(t269);
t308 = rSges(5,1) * t260 + rSges(5,2) * t262 + t387;
t432 = t308 * t263;
t433 = t308 * t261;
t444 = (t261 * t433 + t263 * t432) * t413;
t373 = -t444 + (-t122 * t261 - t441) * t412;
t178 = t212 * t261;
t206 = t261 * t232;
t126 = t206 - t178;
t382 = (-t126 * t261 + t441) * t412 + t444;
t23 = t382 - t373;
t445 = t23 * qJD(1);
t367 = Icges(6,4) * t255;
t211 = Icges(6,1) * t256 - t367;
t149 = -Icges(6,5) * t263 + t211 * t261;
t344 = t256 * t261;
t119 = t149 * t344;
t343 = t256 * t263;
t347 = t255 * t263;
t146 = Icges(6,5) * t343 - Icges(6,6) * t347 + Icges(6,3) * t261;
t223 = Icges(6,4) * t347;
t150 = Icges(6,1) * t343 + Icges(6,5) * t261 - t223;
t148 = Icges(6,4) * t343 - Icges(6,2) * t347 + Icges(6,6) * t261;
t362 = t148 * t255;
t289 = t150 * t256 - t362;
t443 = -t146 * t261 - t289 * t263 - t119;
t207 = Icges(6,5) * t256 - Icges(6,6) * t255;
t354 = t207 * t261;
t145 = -Icges(6,3) * t263 + t354;
t442 = t145 * t261 + t149 * t343;
t258 = t261 ^ 2;
t259 = t263 ^ 2;
t311 = t258 + t259;
t440 = t212 * t311;
t249 = Icges(6,4) * t256;
t209 = -Icges(6,2) * t255 + t249;
t426 = Icges(6,1) * t255 + t249;
t439 = t209 + t426;
t274 = cos(qJ(3));
t438 = Icges(4,5) * t272 + Icges(5,5) * t260 + Icges(4,6) * t274 + Icges(5,6) * t262;
t396 = -t261 / 0.2e1;
t435 = t261 / 0.2e1;
t393 = -t263 / 0.2e1;
t154 = t261 * t178;
t179 = t212 * t263;
t305 = t263 * t179 + t154;
t434 = m(6) * t305;
t368 = Icges(5,4) * t260;
t215 = Icges(5,2) * t262 + t368;
t218 = Icges(5,1) * t262 - t368;
t369 = Icges(4,4) * t272;
t235 = Icges(4,2) * t274 + t369;
t238 = Icges(4,1) * t274 - t369;
t431 = -(t218 / 0.2e1 - t215 / 0.2e1) * t260 - (t238 / 0.2e1 - t235 / 0.2e1) * t272;
t379 = rSges(6,1) * t256;
t213 = -rSges(6,2) * t255 + t379;
t351 = t213 * t263;
t352 = t213 * t261;
t291 = Icges(6,5) * t255 + Icges(6,6) * t256;
t166 = t261 * t291;
t167 = t291 * t263;
t325 = Icges(6,2) * t343 - t150 + t223;
t208 = Icges(6,2) * t256 + t367;
t326 = -t208 * t261 + t149;
t327 = t263 * t426 + t148;
t147 = -Icges(6,6) * t263 + t209 * t261;
t328 = -t261 * t426 - t147;
t421 = t255 * (-t325 * t261 - t326 * t263) + t256 * (t327 * t261 + t328 * t263);
t384 = (t258 * t167 + (-t261 * t166 + t421) * t263) * t396 + (-t259 * t166 + (t263 * t167 - t421) * t261) * t393;
t254 = t263 * pkin(6);
t266 = t274 * pkin(3);
t257 = t266 + pkin(2);
t271 = -qJ(4) - pkin(6);
t334 = t263 * t271;
t133 = t261 * (t334 + t254 + (-pkin(2) + t257) * t261);
t348 = t255 * t261;
t314 = -rSges(6,2) * t348 - t263 * rSges(6,3);
t138 = t261 * (rSges(6,1) * t344 + t314);
t233 = t263 * t257;
t144 = pkin(2) * t263 - t233 + (pkin(6) + t271) * t261;
t301 = rSges(6,1) * t343 - rSges(6,2) * t347;
t153 = rSges(6,3) * t261 + t301;
t386 = pkin(4) * t262;
t226 = t257 + t386;
t205 = t263 * t226;
t47 = t133 + t138 + (t226 - t257) * t258 + (-t144 + t205 - t233 + t153) * t263;
t6 = t384 + m(6) * (t122 * t352 + t124 * t351 - t305 * t47);
t429 = t6 * qJD(5);
t428 = t438 * t261;
t427 = t438 * t263;
t250 = Icges(5,4) * t262;
t216 = -Icges(5,2) * t260 + t250;
t425 = Icges(5,1) * t260 + t250;
t265 = Icges(4,4) * t274;
t236 = -Icges(4,2) * t272 + t265;
t424 = Icges(4,1) * t272 + t265;
t420 = t439 * t255 + (t208 - t211) * t256;
t341 = t260 * t263;
t229 = Icges(5,4) * t341;
t335 = t262 * t263;
t160 = Icges(5,1) * t335 + Icges(5,5) * t261 - t229;
t321 = Icges(5,2) * t335 - t160 + t229;
t159 = -Icges(5,5) * t263 + t218 * t261;
t322 = -t215 * t261 + t159;
t158 = Icges(5,4) * t335 - Icges(5,2) * t341 + Icges(5,6) * t261;
t323 = t263 * t425 + t158;
t157 = -Icges(5,6) * t263 + t216 * t261;
t324 = -t261 * t425 - t157;
t419 = -t260 * (-t321 * t261 - t322 * t263) - t262 * (t323 * t261 + t324 * t263);
t333 = t263 * t272;
t245 = Icges(4,4) * t333;
t332 = t263 * t274;
t177 = Icges(4,1) * t332 + Icges(4,5) * t261 - t245;
t317 = -Icges(4,2) * t332 + t177 - t245;
t175 = Icges(4,4) * t332 - Icges(4,2) * t333 + Icges(4,6) * t261;
t319 = t263 * t424 + t175;
t418 = -t272 * t317 - t274 * t319;
t176 = -Icges(4,5) * t263 + t238 * t261;
t318 = -t235 * t261 + t176;
t174 = -Icges(4,6) * t263 + t236 * t261;
t320 = t261 * t424 + t174;
t417 = -t272 * t318 - t274 * t320;
t304 = t439 * t256 / 0.2e1 + (-t208 / 0.2e1 + t211 / 0.2e1) * t255;
t120 = t150 * t344;
t361 = t149 * t256;
t363 = t147 * t255;
t65 = -t145 * t263 - t147 * t348 + t119;
t66 = t146 * t263 + t148 * t348 - t120;
t306 = ((t120 + (t145 - t362) * t261 - t442) * t261 + ((t145 + t289) * t263 + (t361 + t363) * t261 + t443) * t263) * t396 + (-t261 * t66 - t263 * t65) * t435 + ((-t66 + (t146 - t361) * t263 + t442) * t263 + (t65 + (t146 + t363) * t261 + t443) * t261) * t393;
t416 = 0.4e1 * qJD(1);
t415 = 2 * qJD(3);
t414 = m(4) / 0.2e1;
t381 = rSges(4,1) * t274;
t309 = pkin(2) + t381;
t337 = t261 * t272;
t312 = -rSges(4,2) * t337 - t263 * rSges(4,3);
t388 = sin(qJ(1)) * pkin(1);
t128 = t309 * t261 - t254 + t312 + t388;
t247 = rSges(4,2) * t333;
t267 = cos(qJ(1)) * pkin(1);
t129 = -t247 + t267 + t309 * t263 + (rSges(4,3) + pkin(6)) * t261;
t239 = rSges(4,1) * t272 + rSges(4,2) * t274;
t201 = t239 * t261;
t202 = t239 * t263;
t411 = m(4) * (-t128 * t201 - t129 * t202);
t342 = t260 * t261;
t313 = -rSges(5,2) * t342 - t263 * rSges(5,3);
t380 = rSges(5,1) * t262;
t113 = t388 + t334 + (t257 + t380) * t261 + t313;
t302 = rSges(5,1) * t335 - rSges(5,2) * t341;
t114 = t233 + t267 + (rSges(5,3) - t271) * t261 + t302;
t409 = m(5) * (-t113 * t433 - t114 * t432);
t408 = m(5) * (-t113 * t263 + t114 * t261);
t268 = -pkin(7) + t271;
t109 = t388 + t263 * t268 + (t226 + t379) * t261 + t314;
t110 = t205 + t267 + (rSges(6,3) - t268) * t261 + t301;
t300 = t109 * t351 - t110 * t352;
t403 = m(6) * (t122 * t179 - t124 * t178 + t300);
t402 = m(6) * ((t124 * t261 + t126 * t263) * t212 + t300);
t401 = m(6) * (t109 * t126 - t110 * t124);
t400 = m(6) * (-t109 * t178 - t110 * t179);
t399 = m(6) * (-t109 * t263 + t110 * t261);
t392 = t263 / 0.2e1;
t360 = t157 * t260;
t359 = t158 * t260;
t358 = t159 * t262;
t357 = t174 * t272;
t356 = t175 * t272;
t355 = t176 * t274;
t338 = t261 * t262;
t336 = t261 * t274;
t286 = -m(6) * t440 / 0.2e1;
t75 = -t434 / 0.2e1 + t286;
t330 = t75 * qJD(1);
t307 = t213 + t386;
t295 = Icges(4,5) * t274 - Icges(4,6) * t272;
t293 = Icges(5,5) * t262 - Icges(5,6) * t260;
t288 = t160 * t262 - t359;
t287 = t177 * t274 - t356;
t282 = -t306 + (t327 * t255 + t325 * t256 + t420 * t263 - t354) * t396 + (-t207 * t263 + t328 * t255 + t326 * t256 - t420 * t261) * t393;
t280 = -t304 + (t392 + t393) * (t148 * t256 + t150 * t255);
t248 = pkin(3) * t332;
t241 = -rSges(4,2) * t272 + t381;
t220 = -rSges(5,2) * t260 + t380;
t173 = Icges(4,5) * t332 - Icges(4,6) * t333 + Icges(4,3) * t261;
t172 = -Icges(4,3) * t263 + t295 * t261;
t165 = t220 * t263 + t248;
t163 = (-t220 - t266) * t261;
t156 = Icges(5,5) * t335 - Icges(5,6) * t341 + Icges(5,3) * t261;
t155 = -Icges(5,3) * t263 + t293 * t261;
t143 = t174 * t333;
t142 = t177 * t336;
t141 = t176 * t336;
t132 = t157 * t341;
t131 = t160 * t338;
t130 = t159 * t338;
t125 = t307 * t263 + t248;
t123 = (-t307 - t266) * t261;
t116 = -t201 * t261 - t202 * t263;
t103 = qJD(5) * t434;
t98 = t263 * t153 + t138;
t95 = t308 * t311;
t83 = t173 * t261 + t287 * t263;
t82 = -t172 * t261 - t176 * t332 + t143;
t81 = t173 * t263 + t175 * t337 - t142;
t80 = -t172 * t263 - t174 * t337 + t141;
t76 = t434 / 0.2e1 + t286;
t74 = t156 * t261 + t288 * t263;
t73 = -t155 * t261 - t159 * t335 + t132;
t72 = t156 * t263 + t158 * t342 - t131;
t71 = -t155 * t263 - t157 * t342 + t130;
t60 = t261 * t206 - t154 + (-t311 + t258) * t387 + (-t179 + (t232 + t387) * t263) * t263;
t53 = -t261 * t83 - t263 * t82;
t52 = -t261 * t81 - t263 * t80;
t49 = -t261 * t74 - t263 * t73;
t48 = -t261 * t72 - t263 * t71;
t44 = t399 + t408;
t37 = t304 + t400;
t35 = t402 / 0.2e1;
t33 = t403 / 0.2e1;
t25 = t373 + t382;
t22 = (t143 - t81 + (t173 - t355) * t263) * t263 + (-t141 + t80 + (t173 + t357) * t261) * t261;
t21 = (t82 + t142 - t143 + (t172 - t356) * t261) * t261 + (-t141 - t83 + (t172 + t287) * t263 + (t355 + t357) * t261) * t263;
t19 = (t132 - t72 + (t156 - t358) * t263) * t263 + (-t130 + t71 + (t156 + t360) * t261) * t261;
t18 = (t73 + t131 - t132 + (t155 - t359) * t261) * t261 + (-t130 - t74 + (t155 + t288) * t263 + (t358 + t360) * t261) * t263;
t17 = (t424 / 0.2e1 + t236 / 0.2e1) * t274 + (t425 / 0.2e1 + t216 / 0.2e1) * t262 + t411 + t409 + t401 + t304 - t431;
t8 = m(6) * (t213 * t440 - t98 * t305) + t384;
t7 = t8 * qJD(5);
t4 = t33 - t402 / 0.2e1 + t306;
t3 = t35 - t403 / 0.2e1 + t306;
t2 = t33 + t35 + t282;
t1 = (-t22 / 0.2e1 - t53 / 0.2e1 - t19 / 0.2e1 - t49 / 0.2e1) * t263 + (t52 / 0.2e1 - t21 / 0.2e1 + t48 / 0.2e1 - t18 / 0.2e1) * t261 + t306;
t5 = [qJD(3) * t17 + qJD(4) * t44 + qJD(5) * t37, 0, t17 * qJD(1) + t25 * qJD(4) + t2 * qJD(5) + ((t113 * t165 + t114 * t163) * t413 + (t109 * t125 + t110 * t123 + (t122 + t126) * t124) * t412 + ((t128 * t263 - t129 * t261) * t241 + (-t201 * t263 + t202 * t261) * t239) * t414) * t415 + ((t295 + t293) * (t259 / 0.2e1 + t258 / 0.2e1) + t282 + (t21 + t18) * t435 + (t324 * t260 + t322 * t262 - t320 * t272 + t318 * t274) * t393 + (t323 * t260 + t321 * t262 + t319 * t272 - t317 * t274 + t48 + t52) * t396 + (t22 + t53 + t19 + t49) * t392) * qJD(3), qJD(1) * t44 + qJD(3) * t25 + qJD(5) * t76, t37 * qJD(1) + t2 * qJD(3) + t76 * qJD(4) + (t282 + ((-t178 * t263 + t179 * t261) * t212 + t300) * m(6)) * qJD(5); 0, 0, -t103 + (t116 * t414 + t60 * t412 - t95 * t413) * t415, 0, -qJD(3) * t434 - t103; (t280 - (t216 + t425) * t262 / 0.2e1 - (t236 + t424) * t274 / 0.2e1 + t431) * qJD(1) + t1 * qJD(3) - t23 * qJD(4) + t4 * qJD(5) + (-t411 / 0.4e1 - t409 / 0.4e1 - t401 / 0.4e1) * t416, 0, t1 * qJD(1) + (m(6) * (-t122 * t123 + t124 * t125 + t47 * t60) + m(5) * (-t433 * t163 + t432 * t165 - (t133 + t261 * (rSges(5,1) * t338 + t313) + (rSges(5,3) * t261 - t144 + t302) * t263) * t95) + m(4) * (t239 * t241 * t311 + (t263 * (rSges(4,1) * t332 + rSges(4,3) * t261 - t247) + t261 * (rSges(4,1) * t336 + t312)) * t116) + t384 + ((t417 * t263 + (-t418 - t428) * t261 - t419) * t263 + t427 * t258) * t396 + ((t418 * t261 + (-t417 + t427) * t263 + t419) * t261 - t428 * t259) * t393) * qJD(3) + t429, -t445, t4 * qJD(1) + t6 * qJD(3) + t429; (-t408 / 0.4e1 - t399 / 0.4e1) * t416 + t23 * qJD(3) - t75 * qJD(5), 0, t445 + ((-t123 * t263 - t125 * t261) * t412 + (-t163 * t263 - t165 * t261) * t413) * t415, 0, -t330; (t280 - t400) * qJD(1) + t3 * qJD(3) + t75 * qJD(4) + t306 * qJD(5), 0, t3 * qJD(1) + ((t60 * t98 + (-t123 * t261 + t125 * t263) * t212) * m(6) + t384) * qJD(3) + t7, t330, qJD(1) * t306 + qJD(3) * t8 + t7;];
Cq = t5;
