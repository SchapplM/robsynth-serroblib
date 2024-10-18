% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR15_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:02
% EndTime: 2024-09-27 22:24:14
% DurationCPUTime: 7.52s
% Computational Cost: add. (18906->591), mult. (62558->876), div. (0->0), fcn. (50214->12), ass. (0->295)
t290 = sin(qJ(2));
t286 = cos(pkin(5));
t373 = pkin(1) * t286;
t276 = t290 * t373;
t284 = sin(pkin(5));
t294 = cos(qJ(2));
t338 = t284 * t294;
t388 = pkin(8) + pkin(9);
t227 = (t388 * t338 + t276) * qJD(1);
t293 = cos(qJ(3));
t220 = t293 * t227;
t277 = t294 * t373;
t272 = qJD(1) * t277;
t323 = t284 * t388;
t315 = t290 * t323;
t226 = -qJD(1) * t315 + t272;
t289 = sin(qJ(3));
t169 = -t226 * t289 - t220;
t242 = (-t289 * t290 + t293 * t294) * t284;
t237 = qJD(1) * t242;
t370 = pkin(10) * t237;
t155 = t169 - t370;
t217 = t289 * t227;
t170 = t293 * t226 - t217;
t243 = (t289 * t294 + t290 * t293) * t284;
t238 = qJD(1) * t243;
t235 = t238 * pkin(10);
t156 = -t235 + t170;
t280 = pkin(2) * t293 + pkin(3);
t288 = sin(qJ(4));
t292 = cos(qJ(4));
t326 = qJD(4) * t292;
t327 = qJD(4) * t288;
t335 = t288 * t289;
t403 = -t288 * t155 - t292 * t156 + t280 * t326 + (-t289 * t327 + (t292 * t293 - t335) * qJD(3)) * pkin(2);
t184 = t237 * t288 + t238 * t292;
t285 = cos(pkin(6));
t369 = pkin(11) * t285;
t178 = t184 * t369;
t428 = t178 + t403;
t291 = cos(qJ(5));
t316 = t292 * t237 - t238 * t288;
t287 = sin(qJ(5));
t337 = t285 * t287;
t419 = -t184 * t337 + t291 * t316;
t336 = t285 * t291;
t420 = -t184 * t336 - t287 * t316;
t427 = Ifges(6,1) * t419 + Ifges(6,4) * t420;
t426 = Ifges(6,4) * t419 + Ifges(6,2) * t420;
t425 = Ifges(6,5) * t419 + Ifges(6,6) * t420;
t398 = qJD(2) + qJD(3);
t197 = t398 * t242;
t191 = qJD(1) * t197;
t198 = t398 * t243;
t192 = qJD(1) * t198;
t108 = -qJD(4) * t184 - t191 * t288 - t192 * t292;
t105 = Ifges(5,6) * t108;
t107 = qJD(4) * t316 + t191 * t292 - t192 * t288;
t106 = Ifges(5,5) * t107;
t274 = t286 * qJD(1) + qJD(2);
t269 = qJD(3) + t274;
t262 = qJD(4) + t269;
t283 = sin(pkin(6));
t345 = t316 * t285;
t305 = t262 * t283 + t345;
t123 = -t184 * t287 + t291 * t305;
t124 = t184 * t291 + t287 * t305;
t180 = Ifges(5,4) * t316;
t129 = Ifges(5,1) * t184 + Ifges(5,5) * t262 + t180;
t346 = t316 * t283;
t164 = t262 * t285 + qJD(5) - t346;
t260 = (-pkin(2) * t294 - pkin(1)) * t284;
t255 = qJD(1) * t260;
t201 = -t237 * pkin(3) + t255;
t204 = pkin(2) * t274 + t226;
t165 = t293 * t204 - t217;
t146 = t165 - t235;
t141 = pkin(3) * t269 + t146;
t166 = t204 * t289 + t220;
t147 = t166 + t370;
t264 = qJD(2) * t272;
t306 = qJD(2) * t315;
t213 = -qJD(1) * t306 + t264;
t229 = (-t294 * t323 - t276) * qJD(2);
t214 = qJD(1) * t229;
t328 = qJD(3) * t293;
t329 = qJD(3) * t289;
t113 = t204 * t328 + t293 * t213 + t289 * t214 - t227 * t329;
t87 = -pkin(10) * t192 + t113;
t114 = -qJD(3) * t166 - t213 * t289 + t293 * t214;
t88 = -pkin(10) * t191 + t114;
t32 = t141 * t326 - t147 * t327 + t288 * t88 + t292 * t87;
t17 = t108 * t369 + t32;
t421 = t184 * t283;
t111 = -pkin(4) * t316 - pkin(11) * t421 + t201;
t142 = t288 * t147;
t78 = t292 * t141 - t142;
t68 = t78 - t178;
t66 = pkin(4) * t262 + t68;
t311 = t111 * t283 + t285 * t66;
t144 = t292 * t147;
t79 = t141 * t288 + t144;
t65 = pkin(11) * t305 + t79;
t19 = -t287 * t65 + t291 * t311;
t33 = -qJD(4) * t79 - t288 * t87 + t292 * t88;
t18 = -t107 * t369 + t33;
t331 = qJD(2) * t284;
t318 = qJD(1) * t331;
t314 = t290 * t318;
t171 = pkin(2) * t314 + pkin(3) * t192;
t282 = t283 * pkin(11);
t54 = -pkin(4) * t108 - t107 * t282 + t171;
t313 = t18 * t285 + t283 * t54;
t3 = qJD(5) * t19 + t17 * t291 + t287 * t313;
t325 = qJD(5) * t283;
t339 = t283 * t291;
t340 = t283 * t287;
t361 = Ifges(6,4) * t291;
t362 = Ifges(6,4) * t287;
t364 = Ifges(5,4) * t184;
t8 = -t18 * t283 + t285 * t54;
t368 = t283 * t8;
t363 = Ifges(6,4) * t124;
t63 = Ifges(6,2) * t123 + Ifges(6,6) * t164 + t363;
t390 = -t63 / 0.2e1;
t43 = -qJD(5) * t124 - t107 * t287 + t108 * t336;
t391 = t43 / 0.2e1;
t42 = qJD(5) * t123 + t107 * t291 + t108 * t337;
t392 = t42 / 0.2e1;
t349 = t108 * t283;
t416 = -t349 / 0.2e1;
t393 = Ifges(6,1) * t392 + Ifges(6,4) * t391 + Ifges(6,5) * t416;
t394 = Ifges(6,4) * t392 + Ifges(6,2) * t391 + Ifges(6,6) * t416;
t12 = Ifges(6,5) * t42 + Ifges(6,6) * t43 - Ifges(6,3) * t349;
t395 = t12 / 0.2e1;
t20 = t287 * t311 + t291 * t65;
t4 = -qJD(5) * t20 - t17 * t287 + t313 * t291;
t405 = t33 * mrSges(5,1) - t32 * mrSges(5,2);
t47 = t111 * t285 - t283 * t66;
t121 = Ifges(6,4) * t123;
t64 = Ifges(6,1) * t124 + Ifges(6,5) * t164 + t121;
t424 = t3 * (-mrSges(6,2) * t285 + mrSges(6,3) * t339) + t339 * t394 + t105 + t106 + t285 * t395 + t4 * (mrSges(6,1) * t285 - mrSges(6,3) * t340) + (Ifges(6,5) * t285 + (Ifges(6,1) * t287 + t361) * t283) * t392 + (Ifges(6,6) * t285 + (Ifges(6,2) * t291 + t362) * t283) * t391 + (-mrSges(6,1) * t291 + mrSges(6,2) * t287) * t368 + (Ifges(6,3) * t285 + (Ifges(6,5) * t287 + Ifges(6,6) * t291) * t283) * t416 + t340 * t393 + t405 + (t47 * (mrSges(6,1) * t287 + mrSges(6,2) * t291) + t287 * t390) * t325 + (t123 * (-Ifges(6,2) * t287 + t361) + t124 * (Ifges(6,1) * t291 - t362) + t164 * (Ifges(6,5) * t291 - Ifges(6,6) * t287) + t291 * t64) * t325 / 0.2e1 - (Ifges(5,5) * t316 - Ifges(5,6) * t184) * t262 / 0.2e1 - (-Ifges(5,2) * t184 + t129 + t180) * t316 / 0.2e1 - t201 * (mrSges(5,1) * t184 + mrSges(5,2) * t316) - t419 * t64 / 0.2e1 + t420 * t390 - (Ifges(5,1) * t316 - t364) * t184 / 0.2e1 - t47 * (-mrSges(6,1) * t420 + mrSges(6,2) * t419);
t423 = -t294 / 0.2e1;
t422 = pkin(4) * t184;
t356 = t184 * t79;
t128 = Ifges(5,2) * t316 + Ifges(5,6) * t262 + t364;
t417 = t128 / 0.2e1;
t412 = t316 * t78;
t265 = pkin(3) * t288 + t282;
t279 = pkin(3) * t292 + pkin(4);
t231 = t265 * t291 + t279 * t337;
t371 = pkin(3) * t238;
t120 = -t282 * t316 + t371 + t422;
t324 = t316 * t369;
t85 = -t146 * t288 - t144;
t71 = t85 - t324;
t309 = t120 * t283 + t285 * t71;
t358 = pkin(3) * qJD(4);
t86 = t292 * t146 - t142;
t72 = -t178 + t86;
t411 = t287 * t72 - t291 * t309 - t231 * qJD(5) + (-t287 * t292 - t288 * t336) * t358;
t230 = -t265 * t287 + t279 * t336;
t410 = -t287 * t309 - t291 * t72 + t230 * qJD(5) + (-t288 * t337 + t291 * t292) * t358;
t258 = pkin(4) * t337 + pkin(11) * t339;
t136 = -pkin(11) * t346 + t422;
t69 = -pkin(11) * t345 - t79;
t307 = t136 * t283 + t285 * t69;
t409 = -t258 * qJD(5) + t287 * t68 - t291 * t307;
t256 = pkin(4) * t336 - pkin(11) * t340;
t408 = t256 * qJD(5) - t287 * t307 - t291 * t68;
t334 = t289 * t292;
t254 = pkin(2) * t334 + t288 * t280;
t239 = t282 + t254;
t253 = -pkin(2) * t335 + t292 * t280;
t252 = pkin(4) + t253;
t188 = t239 * t291 + t252 * t337;
t203 = -t280 * t327 + (-t289 * t326 + (-t288 * t293 - t334) * qJD(3)) * pkin(2);
t332 = qJD(1) * t284;
t321 = t290 * t332;
t270 = pkin(2) * t321;
t115 = t120 + t270;
t90 = t292 * t155 - t156 * t288;
t73 = t90 - t324;
t310 = t115 * t283 + t285 * t73;
t407 = -qJD(5) * t188 + t203 * t336 - t428 * t287 - t291 * t310;
t187 = -t239 * t287 + t252 * t336;
t406 = qJD(5) * t187 + t203 * t337 - t287 * t310 + t428 * t291;
t404 = -t90 + t203;
t400 = t114 * mrSges(4,1) - t113 * mrSges(4,2);
t225 = pkin(2) * t286 + t277 - t315;
t333 = pkin(8) * t338 + t276;
t236 = pkin(9) * t338 + t333;
t175 = t293 * t225 - t236 * t289;
t157 = pkin(3) * t286 - pkin(10) * t243 + t175;
t176 = t289 * t225 + t293 * t236;
t161 = pkin(10) * t242 + t176;
t100 = t288 * t157 + t292 * t161;
t240 = -pkin(8) * t314 + t264;
t251 = t333 * qJD(2);
t241 = qJD(1) * t251;
t399 = -t241 * mrSges(3,1) - t240 * mrSges(3,2);
t246 = -pkin(8) * t321 + t272;
t247 = t333 * qJD(1);
t320 = t294 * t332;
t268 = Ifges(3,4) * t320;
t359 = Ifges(3,2) * t294;
t360 = Ifges(3,5) * t294;
t365 = Ifges(3,4) * t290;
t374 = -t290 / 0.2e1;
t397 = (t246 * t294 + t247 * t290) * mrSges(3,3) - t274 * (-Ifges(3,6) * t290 + t360) / 0.2e1 - (Ifges(3,6) * t274 + (t359 + t365) * t332) * t374 + (Ifges(3,1) * t321 + t274 * Ifges(3,5) + t268) * t423;
t196 = t242 * t288 + t243 * t292;
t195 = t242 * t292 - t243 * t288;
t304 = t195 * t285 + t283 * t286;
t135 = t196 * t291 + t287 * t304;
t210 = -t242 * pkin(3) + t260;
t122 = -t195 * pkin(4) - t196 * t282 + t210;
t99 = t292 * t157 - t161 * t288;
t77 = pkin(4) * t286 - t196 * t369 + t99;
t308 = t122 * t283 + t285 * t77;
t76 = pkin(11) * t304 + t100;
t35 = t287 * t308 + t291 * t76;
t396 = -0.2e1 * pkin(1);
t387 = -t123 / 0.2e1;
t386 = -t124 / 0.2e1;
t385 = t124 / 0.2e1;
t384 = -t164 / 0.2e1;
t381 = t184 / 0.2e1;
t379 = -t238 / 0.2e1;
t378 = t238 / 0.2e1;
t376 = -t283 / 0.2e1;
t375 = t286 / 0.2e1;
t372 = pkin(2) * t290;
t366 = mrSges(4,3) * t237;
t355 = t238 * mrSges(4,3);
t354 = t238 * Ifges(4,4);
t117 = -qJD(4) * t196 - t197 * t288 - t198 * t292;
t348 = t117 * t283;
t330 = qJD(2) * t290;
t62 = Ifges(6,5) * t124 + Ifges(6,6) * t123 + Ifges(6,3) * t164;
t322 = t62 * t376;
t319 = t284 * t330;
t271 = pkin(2) * t319;
t179 = pkin(3) * t198 + t271;
t116 = qJD(4) * t195 + t197 * t292 - t198 * t288;
t273 = qJD(2) * t277;
t228 = t273 - t306;
t125 = t225 * t328 + t293 * t228 + t289 * t229 - t236 * t329;
t101 = -pkin(10) * t198 + t125;
t126 = -qJD(3) * t176 - t228 * t289 + t293 * t229;
t102 = -pkin(10) * t197 + t126;
t37 = -qJD(4) * t100 - t101 * t288 + t292 * t102;
t22 = -t116 * t369 + t37;
t61 = -pkin(4) * t117 - t116 * t282 + t179;
t312 = t22 * t285 + t283 * t61;
t36 = t292 * t101 + t288 * t102 + t157 * t326 - t161 * t327;
t300 = qJD(5) * (-t19 * t291 - t20 * t287);
t298 = mrSges(6,3) * t300;
t34 = -t287 * t76 + t291 * t308;
t134 = -t196 * t287 + t291 * t304;
t173 = t237 * Ifges(4,2) + t269 * Ifges(4,6) + t354;
t234 = Ifges(4,4) * t237;
t174 = t238 * Ifges(4,1) + t269 * Ifges(4,5) + t234;
t185 = Ifges(4,6) * t192;
t186 = Ifges(4,5) * t191;
t295 = t165 * t366 + mrSges(5,3) * t412 + t173 * t378 + (Ifges(4,1) * t237 - t354) * t379 + t424 + t186 - t185 - t255 * (mrSges(4,1) * t238 + mrSges(4,2) * t237) - t269 * (Ifges(4,5) * t237 - Ifges(4,6) * t238) / 0.2e1 + t400 - (-Ifges(4,2) * t238 + t174 + t234) * t237 / 0.2e1 + t184 * t417 + (Ifges(6,3) * t421 + t425) * t384 + (Ifges(6,5) * t421 + t427) * t386 + (Ifges(6,6) * t421 + t426) * t387 - t20 * (-mrSges(6,2) * t421 + mrSges(6,3) * t420) - t19 * (mrSges(6,1) * t421 - mrSges(6,3) * t419) - t62 * t421 / 0.2e1;
t261 = t318 * t360;
t257 = -pkin(8) * t284 * t290 + t277;
t250 = -pkin(8) * t319 + t273;
t245 = -t274 * mrSges(3,2) + mrSges(3,3) * t320;
t244 = mrSges(3,1) * t274 - mrSges(3,3) * t321;
t211 = t270 + t371;
t206 = mrSges(4,1) * t269 - t355;
t205 = -mrSges(4,2) * t269 + t366;
t194 = -mrSges(4,1) * t237 + mrSges(4,2) * t238;
t177 = -t195 * t283 + t285 * t286;
t168 = mrSges(5,1) * t262 - mrSges(5,3) * t184;
t167 = -mrSges(5,2) * t262 + mrSges(5,3) * t316;
t140 = -mrSges(5,1) * t316 + mrSges(5,2) * t184;
t81 = mrSges(6,1) * t164 - mrSges(6,3) * t124;
t80 = -mrSges(6,2) * t164 + mrSges(6,3) * t123;
t67 = -mrSges(6,1) * t123 + mrSges(6,2) * t124;
t57 = t122 * t285 - t283 * t77;
t56 = t136 * t285 - t283 * t69;
t55 = t115 * t285 - t283 * t73;
t53 = t120 * t285 - t283 * t71;
t51 = -qJD(5) * t135 - t116 * t287 + t117 * t336;
t50 = qJD(5) * t134 + t116 * t291 + t117 * t337;
t28 = mrSges(6,2) * t349 + mrSges(6,3) * t43;
t27 = -mrSges(6,1) * t349 - mrSges(6,3) * t42;
t21 = t117 * t369 + t36;
t16 = -mrSges(6,1) * t43 + mrSges(6,2) * t42;
t15 = -t22 * t283 + t285 * t61;
t6 = -qJD(5) * t35 - t21 * t287 + t312 * t291;
t5 = qJD(5) * t34 + t21 * t291 + t287 * t312;
t1 = [t117 * t417 + t34 * t27 + t35 * t28 + (Ifges(5,1) * t116 + Ifges(5,4) * t117) * t381 + (Ifges(6,1) * t50 + Ifges(6,4) * t51 - Ifges(6,5) * t348) * t385 + (Ifges(6,4) * t135 + Ifges(6,2) * t134 + Ifges(6,6) * t177) * t391 + (Ifges(6,1) * t135 + Ifges(6,4) * t134 + Ifges(6,5) * t177) * t392 + t135 * t393 + t134 * t394 + t177 * t395 + t47 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t57 * t16 + t51 * t63 / 0.2e1 + t50 * t64 / 0.2e1 + m(4) * (t113 * t176 + t114 * t175 + t125 * t166 + t126 * t165 + t255 * t271) + t15 * t67 + t5 * t80 + t6 * t81 + t116 * t129 / 0.2e1 + t8 * (-mrSges(6,1) * t134 + mrSges(6,2) * t135) + t36 * t167 + t37 * t168 + t3 * (-mrSges(6,2) * t177 + mrSges(6,3) * t134) + t4 * (mrSges(6,1) * t177 - mrSges(6,3) * t135) + t179 * t140 + (t240 * t294 + t241 * t290) * t284 * mrSges(3,3) + (-t116 * t78 + t117 * t79 + t195 * t32 - t196 * t33) * mrSges(5,3) + t117 * t322 + m(6) * (t15 * t47 + t19 * t6 + t20 * t5 + t3 * t35 + t34 * t4 + t57 * t8) + m(5) * (t100 * t32 + t171 * t210 + t179 * t201 + t33 * t99 + t36 * t79 + t37 * t78) + t316 * (Ifges(5,4) * t116 + Ifges(5,2) * t117) / 0.2e1 + ((-t257 * mrSges(3,3) + Ifges(3,5) * t375 + (mrSges(3,2) * t396 + 0.3e1 / 0.2e1 * Ifges(3,4) * t294) * t284) * t294 + (-t333 * mrSges(3,3) - Ifges(3,6) * t286 + (mrSges(3,1) * t396 - 0.3e1 / 0.2e1 * t365 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t294) * t284 + (m(4) * t260 - mrSges(4,1) * t242 + mrSges(4,2) * t243) * pkin(2)) * t290) * t318 + (t261 / 0.2e1 + t106 / 0.2e1 + t105 / 0.2e1 + t186 / 0.2e1 - t185 / 0.2e1 + t399 + t400 + t405) * t286 - (-t260 * mrSges(4,1) + t176 * mrSges(4,3) + Ifges(4,4) * t243 + Ifges(4,2) * t242 + Ifges(4,6) * t375) * t192 + (t210 * mrSges(5,2) - t99 * mrSges(5,3) + Ifges(5,1) * t196 + Ifges(5,4) * t195 + Ifges(5,5) * t375) * t107 + (t260 * mrSges(4,2) - t175 * mrSges(4,3) + Ifges(4,1) * t243 + Ifges(4,4) * t242 + Ifges(4,5) * t375) * t191 + ((Ifges(6,5) * t135 + Ifges(6,6) * t134 + Ifges(6,3) * t177) * t376 + t100 * mrSges(5,3) + Ifges(5,4) * t196 + Ifges(5,2) * t195 + Ifges(5,6) * t375 - t210 * mrSges(5,1)) * t108 + t164 * (Ifges(6,5) * t50 + Ifges(6,6) * t51 - Ifges(6,3) * t348) / 0.2e1 + t123 * (Ifges(6,4) * t50 + Ifges(6,2) * t51 - Ifges(6,6) * t348) / 0.2e1 + t19 * (-mrSges(6,1) * t348 - mrSges(6,3) * t50) + t20 * (mrSges(6,2) * t348 + mrSges(6,3) * t51) + t262 * (Ifges(5,5) * t116 + Ifges(5,6) * t117) / 0.2e1 + t250 * t245 - t251 * t244 + t201 * (-mrSges(5,1) * t117 + mrSges(5,2) * t116) + t125 * t205 + t126 * t206 + t171 * (-mrSges(5,1) * t195 + mrSges(5,2) * t196) + t197 * t174 / 0.2e1 - t198 * t173 / 0.2e1 + (t194 * t372 - t397) * t331 + (Ifges(4,1) * t197 - Ifges(4,4) * t198) * t378 + (t113 * t242 - t114 * t243 - t165 * t197 - t166 * t198) * mrSges(4,3) + t269 * (Ifges(4,5) * t197 - Ifges(4,6) * t198) / 0.2e1 + t255 * (mrSges(4,1) * t198 + mrSges(4,2) * t197) + t237 * (Ifges(4,4) * t197 - Ifges(4,2) * t198) / 0.2e1 + m(3) * (t240 * t333 - t241 * t257 - t246 * t251 + t247 * t250); t166 * t355 - m(4) * (t165 * t169 + t166 * t170) - t55 * t67 + t403 * t167 + t404 * t168 + t406 * t80 + (-t252 * t16 - t203 * t67 + t298) * t283 + t295 + t407 * t81 + (-t107 * t253 + t108 * t254 + t356) * mrSges(5,3) + t261 + ((t205 * t293 - t206 * t289) * qJD(3) + (-t191 * t293 - t192 * t289) * mrSges(4,3) + m(4) * (t113 * t289 + t114 * t293 + (-t165 * t289 + t166 * t293) * qJD(3))) * pkin(2) - t246 * t245 + t247 * t244 - t170 * t205 - t169 * t206 - t211 * t140 + t187 * t27 + t188 * t28 + t399 + (t268 * t423 - Ifges(3,6) * t330 + (t290 * t359 / 0.2e1 + (Ifges(3,1) * t294 - t365) * t374 + pkin(1) * (mrSges(3,1) * t290 + mrSges(3,2) * t294)) * t332 + (-m(4) * t255 - t194) * t372 + t397) * t332 + (-t201 * t211 + t253 * t33 + t254 * t32 + t403 * t79 + t404 * t78) * m(5) + (t187 * t4 + t188 * t3 + (-t203 * t47 - t252 * t8) * t283 - t47 * t55 + t406 * t20 + t407 * t19) * m(6); mrSges(5,3) * t356 - m(5) * (t78 * t85 + t79 * t86) - t53 * t67 - t86 * t167 - t85 * t168 + t411 * t81 + t410 * t80 + (-t279 * t16 + t298) * t283 + t295 + (t206 + t355) * t166 + t230 * t27 + t231 * t28 - t165 * t205 + (t19 * t411 + t20 * t410 + t230 * t4 + t231 * t3 - t279 * t368 - t47 * t53) * m(6) + (m(6) * t283 * t47 * t327 - t238 * t140 + (-t107 * t292 + t108 * t288) * mrSges(5,3) + (t167 * t292 + (t283 * t67 - t168) * t288) * qJD(4) + (0.2e1 * t201 * t379 + t288 * t32 + t292 * t33 + t326 * t79 - t327 * t78) * m(5)) * pkin(3); t128 * t381 - t283 * pkin(4) * t16 - t56 * t67 - t78 * t167 + t79 * t168 + t184 * t322 + t424 + t258 * t28 + t256 * t27 + t408 * t80 + t409 * t81 + (-pkin(4) * t368 + t19 * t409 + t20 * t408 + t256 * t4 + t258 * t3 - t47 * t56) * m(6) + (t412 + t356) * mrSges(5,3) + (t19 * t419 - t20 * t420 + t283 * t300) * mrSges(6,3) + t425 * t384 + t426 * t387 + t427 * t386 + (-t19 * mrSges(6,1) + t20 * mrSges(6,2) + Ifges(6,5) * t386 + Ifges(6,6) * t387 + Ifges(6,3) * t384) * t421; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t47 * (mrSges(6,1) * t124 + mrSges(6,2) * t123) + (Ifges(6,1) * t123 - t363) * t386 + t63 * t385 + (Ifges(6,5) * t123 - Ifges(6,6) * t124) * t384 - t19 * t80 + t20 * t81 + (t123 * t19 + t124 * t20) * mrSges(6,3) + t12 + (-Ifges(6,2) * t124 + t121 + t64) * t387;];
tauc = t1(:);
