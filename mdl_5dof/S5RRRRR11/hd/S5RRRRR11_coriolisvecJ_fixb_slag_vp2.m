% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:49
% EndTime: 2019-12-31 22:39:25
% DurationCPUTime: 15.92s
% Computational Cost: add. (12969->728), mult. (33856->1039), div. (0->0), fcn. (26018->10), ass. (0->320)
t267 = sin(qJ(4));
t269 = sin(qJ(2));
t271 = cos(qJ(4));
t264 = sin(pkin(5));
t315 = qJD(1) * t264;
t272 = cos(qJ(3));
t273 = cos(qJ(2));
t324 = t272 * t273;
t193 = (-t267 * t324 + t269 * t271) * t315;
t268 = sin(qJ(3));
t309 = qJD(4) * t271;
t311 = qJD(3) * t272;
t427 = t267 * t311 + t268 * t309 + t193;
t297 = pkin(3) * t268 - pkin(9) * t272;
t244 = t297 * qJD(3);
t249 = -pkin(3) * t272 - pkin(9) * t268 - pkin(2);
t310 = qJD(4) * t267;
t312 = qJD(3) * t268;
t149 = t267 * t244 + t249 * t309 + (-t271 * t312 - t272 * t310) * pkin(8);
t305 = t269 * t315;
t265 = cos(pkin(5));
t365 = pkin(1) * t273;
t308 = t265 * t365;
t223 = -pkin(7) * t305 + qJD(1) * t308;
t286 = (pkin(2) * t269 - pkin(8) * t273) * t264;
t224 = qJD(1) * t286;
t157 = t272 * t223 + t268 * t224;
t143 = pkin(9) * t305 + t157;
t260 = t265 * t269 * pkin(1);
t326 = t264 * t273;
t160 = (t260 + (pkin(7) + t297) * t326) * qJD(1);
t80 = t271 * t143 + t267 * t160;
t426 = -t80 + t149;
t194 = (t267 * t269 + t271 * t324) * t315;
t325 = t271 * t272;
t261 = pkin(8) * t325;
t304 = t273 * t315;
t362 = pkin(8) * t267;
t317 = t271 * t244 + t312 * t362;
t360 = pkin(10) * t268;
t363 = pkin(4) * t268;
t79 = -t267 * t143 + t271 * t160;
t425 = t194 * pkin(10) - t304 * t363 + (-pkin(10) * t325 + t363) * qJD(3) + (-t261 + (-t249 + t360) * t267) * qJD(4) + t317 - t79;
t424 = -t427 * pkin(10) + t426;
t388 = -pkin(10) - pkin(9);
t306 = qJD(4) * t388;
t256 = qJD(1) * t265 + qJD(2);
t205 = t256 * t272 - t268 * t305;
t328 = t205 * t267;
t316 = pkin(7) * t326 + t260;
t226 = t316 * qJD(1);
t186 = t256 * pkin(8) + t226;
t220 = (-pkin(2) * t273 - pkin(8) * t269 - pkin(1)) * t264;
t198 = qJD(1) * t220;
t131 = -t268 * t186 + t198 * t272;
t206 = t256 * t268 + t272 * t305;
t151 = pkin(3) * t206 - pkin(9) * t205;
t69 = t271 * t131 + t267 * t151;
t423 = pkin(10) * t328 + t267 * t306 - t69;
t359 = pkin(10) * t271;
t68 = -t131 * t267 + t271 * t151;
t422 = -pkin(4) * t206 + t205 * t359 + t271 * t306 - t68;
t250 = qJD(3) - t304;
t111 = -pkin(3) * t250 - t131;
t185 = -t256 * pkin(2) - t223;
t109 = -t205 * pkin(3) - t206 * pkin(9) + t185;
t132 = t186 * t272 + t198 * t268;
t112 = pkin(9) * t250 + t132;
t53 = t271 * t109 - t112 * t267;
t54 = t109 * t267 + t112 * t271;
t288 = t267 * t54 + t271 * t53;
t290 = Ifges(5,5) * t271 - Ifges(5,6) * t267;
t350 = Ifges(5,4) * t271;
t292 = -Ifges(5,2) * t267 + t350;
t351 = Ifges(5,4) * t267;
t294 = Ifges(5,1) * t271 - t351;
t295 = mrSges(5,1) * t267 + mrSges(5,2) * t271;
t367 = t271 / 0.2e1;
t368 = -t267 / 0.2e1;
t200 = qJD(4) - t205;
t373 = t200 / 0.2e1;
t162 = t206 * t271 + t250 * t267;
t380 = t162 / 0.2e1;
t161 = -t206 * t267 + t250 * t271;
t382 = t161 / 0.2e1;
t352 = Ifges(5,4) * t162;
t74 = Ifges(5,2) * t161 + Ifges(5,6) * t200 + t352;
t159 = Ifges(5,4) * t161;
t75 = Ifges(5,1) * t162 + Ifges(5,5) * t200 + t159;
t421 = t288 * mrSges(5,3) - t111 * t295 - t290 * t373 - t292 * t382 - t294 * t380 - t367 * t75 - t368 * t74;
t270 = cos(qJ(5));
t266 = sin(qJ(5));
t47 = pkin(10) * t161 + t54;
t331 = t266 * t47;
t46 = -pkin(10) * t162 + t53;
t39 = pkin(4) * t200 + t46;
t13 = t270 * t39 - t331;
t330 = t270 * t47;
t14 = t266 * t39 + t330;
t299 = t270 * t161 - t162 * t266;
t96 = t161 * t266 + t162 * t270;
t366 = Ifges(6,4) * t96;
t189 = qJD(5) + t200;
t377 = -t189 / 0.2e1;
t390 = -t96 / 0.2e1;
t392 = -t299 / 0.2e1;
t92 = Ifges(6,4) * t299;
t45 = Ifges(6,1) * t96 + Ifges(6,5) * t189 + t92;
t70 = -pkin(4) * t161 + t111;
t420 = (Ifges(6,5) * t299 - Ifges(6,6) * t96) * t377 + (t13 * t299 + t14 * t96) * mrSges(6,3) + (-Ifges(6,2) * t96 + t45 + t92) * t392 - t70 * (mrSges(6,1) * t96 + mrSges(6,2) * t299) + (Ifges(6,1) * t299 - t366) * t390;
t313 = qJD(2) * t273;
t303 = t268 * t313;
t168 = t256 * t312 + (t269 * t311 + t303) * t315;
t314 = qJD(2) * t264;
t300 = qJD(1) * t314;
t298 = t269 * t300;
t225 = qJD(2) * t286;
t215 = qJD(1) * t225;
t327 = t264 * t269;
t257 = pkin(7) * t327;
t235 = -t257 + t308;
t227 = t235 * qJD(2);
t216 = qJD(1) * t227;
t71 = -t186 * t312 + t198 * t311 + t268 * t215 + t272 * t216;
t62 = pkin(9) * t298 + t71;
t302 = t272 * t313;
t167 = t256 * t311 + (-t269 * t312 + t302) * t315;
t228 = t316 * qJD(2);
t217 = qJD(1) * t228;
t83 = t168 * pkin(3) - t167 * pkin(9) + t217;
t19 = -qJD(4) * t54 - t267 * t62 + t271 * t83;
t86 = qJD(4) * t161 + t167 * t271 + t267 * t298;
t10 = pkin(4) * t168 - pkin(10) * t86 + t19;
t18 = t109 * t309 - t112 * t310 + t267 * t83 + t271 * t62;
t87 = -qJD(4) * t162 - t167 * t267 + t271 * t298;
t11 = pkin(10) * t87 + t18;
t2 = qJD(5) * t13 + t10 * t266 + t11 * t270;
t3 = -qJD(5) * t14 + t10 * t270 - t11 * t266;
t419 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t239 = t271 * t249;
t171 = -t268 * t359 + t239 + (-pkin(4) - t362) * t272;
t214 = t267 * t249 + t261;
t182 = -t267 * t360 + t214;
t113 = t171 * t270 - t182 * t266;
t418 = qJD(5) * t113 + t425 * t266 + t424 * t270;
t114 = t171 * t266 + t182 * t270;
t417 = -qJD(5) * t114 - t424 * t266 + t425 * t270;
t156 = -t268 * t223 + t224 * t272;
t142 = -pkin(3) * t305 - t156;
t416 = t427 * pkin(4) + pkin(8) * t311 - t142;
t199 = Ifges(4,4) * t205;
t335 = t250 * Ifges(4,5);
t337 = t206 * Ifges(4,1);
t128 = t199 + t335 + t337;
t279 = t131 * mrSges(4,3) - t128 / 0.2e1 - t185 * mrSges(4,2) - t335 / 0.2e1;
t415 = t279 + t421;
t30 = qJD(5) * t299 + t266 * t87 + t270 * t86;
t398 = t30 / 0.2e1;
t31 = -qJD(5) * t96 - t266 * t86 + t270 * t87;
t397 = t31 / 0.2e1;
t44 = Ifges(6,2) * t299 + Ifges(6,6) * t189 + t366;
t413 = t44 / 0.2e1;
t379 = t168 / 0.2e1;
t412 = -Ifges(3,6) * t256 / 0.2e1;
t84 = Ifges(5,6) * t87;
t85 = Ifges(5,5) * t86;
t34 = Ifges(5,3) * t168 + t84 + t85;
t28 = Ifges(6,6) * t31;
t29 = Ifges(6,5) * t30;
t6 = Ifges(6,3) * t168 + t28 + t29;
t411 = t34 + t6;
t410 = -t19 * mrSges(5,1) + t18 * mrSges(5,2) - t419;
t251 = t388 * t267;
t252 = t388 * t271;
t188 = t251 * t266 - t252 * t270;
t405 = -qJD(5) * t188 - t423 * t266 + t422 * t270;
t187 = t251 * t270 + t252 * t266;
t404 = qJD(5) * t187 + t422 * t266 + t423 * t270;
t50 = -mrSges(6,1) * t299 + mrSges(6,2) * t96;
t403 = m(6) * t70 + t50;
t287 = t266 * t267 - t270 * t271;
t230 = t287 * t268;
t218 = t257 + (-pkin(2) - t365) * t265;
t231 = -t265 * t272 + t268 * t327;
t232 = t265 * t268 + t272 * t327;
t139 = t231 * pkin(3) - t232 * pkin(9) + t218;
t219 = pkin(8) * t265 + t316;
t153 = t272 * t219 + t268 * t220;
t141 = -pkin(9) * t326 + t153;
t65 = t267 * t139 + t271 * t141;
t402 = t18 * t271 - t19 * t267;
t401 = qJD(4) + qJD(5);
t400 = Ifges(6,4) * t398 + Ifges(6,2) * t397 + Ifges(6,6) * t379;
t399 = Ifges(6,1) * t398 + Ifges(6,4) * t397 + Ifges(6,5) * t379;
t36 = t86 * Ifges(5,1) + t87 * Ifges(5,4) + t168 * Ifges(5,5);
t396 = t36 / 0.2e1;
t395 = -t74 / 0.2e1;
t394 = t86 / 0.2e1;
t393 = t87 / 0.2e1;
t391 = t299 / 0.2e1;
t389 = t96 / 0.2e1;
t386 = pkin(1) * mrSges(3,1);
t385 = pkin(1) * mrSges(3,2);
t383 = -t161 / 0.2e1;
t381 = -t162 / 0.2e1;
t376 = t189 / 0.2e1;
t375 = -t199 / 0.2e1;
t374 = -t200 / 0.2e1;
t372 = -t231 / 0.2e1;
t370 = t232 / 0.2e1;
t369 = t265 / 0.2e1;
t364 = pkin(4) * t267;
t361 = pkin(8) * t272;
t358 = t71 * mrSges(4,2);
t72 = -t186 * t311 - t198 * t312 + t215 * t272 - t268 * t216;
t357 = t72 * mrSges(4,1);
t356 = t299 * Ifges(6,6);
t355 = t96 * Ifges(6,5);
t354 = Ifges(3,4) * t269;
t353 = Ifges(4,4) * t206;
t348 = t161 * Ifges(5,6);
t347 = t162 * Ifges(5,5);
t346 = t167 * Ifges(4,1);
t345 = t167 * Ifges(4,4);
t344 = t168 * Ifges(4,4);
t342 = t189 * Ifges(6,3);
t340 = t200 * Ifges(5,3);
t339 = t205 * Ifges(4,2);
t334 = t250 * Ifges(4,6);
t101 = -mrSges(5,1) * t161 + mrSges(5,2) * t162;
t170 = mrSges(4,1) * t250 - mrSges(4,3) * t206;
t323 = t101 - t170;
t241 = t266 * t271 + t267 * t270;
t174 = t401 * t241;
t121 = -t174 * t268 - t287 * t311;
t137 = t193 * t266 + t194 * t270;
t322 = t121 - t137;
t122 = t230 * t401 - t241 * t311;
t136 = t193 * t270 - t194 * t266;
t321 = t122 - t136;
t134 = t241 * t205;
t320 = -t134 + t174;
t135 = t287 * t205;
t173 = t401 * t287;
t319 = -t135 + t173;
t318 = -mrSges(3,1) * t256 - mrSges(4,1) * t205 + mrSges(4,2) * t206 + mrSges(3,3) * t305;
t307 = Ifges(4,5) * t167 - Ifges(4,6) * t168 + Ifges(4,3) * t298;
t301 = t269 * t314;
t64 = t271 * t139 - t141 * t267;
t152 = -t268 * t219 + t220 * t272;
t140 = pkin(3) * t326 - t152;
t296 = mrSges(5,1) * t271 - mrSges(5,2) * t267;
t293 = Ifges(5,1) * t267 + t350;
t291 = Ifges(5,2) * t271 + t351;
t289 = Ifges(5,5) * t267 + Ifges(5,6) * t271;
t285 = -t232 * t271 + t267 * t326;
t51 = pkin(4) * t231 + pkin(10) * t285 + t64;
t178 = -t232 * t267 - t271 * t326;
t55 = pkin(10) * t178 + t65;
t22 = -t266 * t55 + t270 * t51;
t23 = t266 * t51 + t270 * t55;
t117 = t178 * t270 + t266 * t285;
t118 = t178 * t266 - t270 * t285;
t89 = -t219 * t311 - t220 * t312 + t225 * t272 - t268 * t227;
t176 = qJD(3) * t232 + t264 * t303;
t177 = -qJD(3) * t231 + t264 * t302;
t102 = t176 * pkin(3) - t177 * pkin(9) + t228;
t88 = -t219 * t312 + t220 * t311 + t268 * t225 + t272 * t227;
t77 = pkin(9) * t301 + t88;
t26 = t267 * t102 + t139 * t309 - t141 * t310 + t271 * t77;
t253 = Ifges(3,4) * t304;
t281 = -t223 * mrSges(3,3) + t256 * Ifges(3,5) + Ifges(3,1) * t305 / 0.2e1 + t253 / 0.2e1;
t78 = -pkin(3) * t301 - t89;
t27 = -qJD(4) * t65 + t271 * t102 - t267 * t77;
t63 = -pkin(3) * t298 - t72;
t278 = t131 * mrSges(4,1) + t250 * Ifges(4,3) + t206 * Ifges(4,5) + t205 * Ifges(4,6) + t412 - (Ifges(3,2) * t273 + t354) * t315 / 0.2e1 - t132 * mrSges(4,2) - t226 * mrSges(3,3);
t127 = t334 + t339 + t353;
t43 = t342 + t355 + t356;
t73 = t340 + t347 + t348;
t276 = -t185 * mrSges(4,1) + t334 / 0.2e1 - t348 / 0.2e1 - t347 / 0.2e1 - t340 / 0.2e1 + t14 * mrSges(6,2) - t13 * mrSges(6,1) + t54 * mrSges(5,2) - t53 * mrSges(5,1) - t43 / 0.2e1 - t73 / 0.2e1 - t356 / 0.2e1 - t355 / 0.2e1 - t342 / 0.2e1 + t127 / 0.2e1 + t132 * mrSges(4,3) + t353 / 0.2e1;
t274 = t276 + t339 / 0.2e1;
t263 = -pkin(4) * t271 - pkin(3);
t248 = Ifges(3,5) * t273 * t300;
t245 = (pkin(8) + t364) * t268;
t229 = t241 * t268;
t222 = -t256 * mrSges(3,2) + mrSges(3,3) * t304;
t213 = -t267 * t361 + t239;
t169 = -mrSges(4,2) * t250 + mrSges(4,3) * t205;
t150 = -qJD(4) * t214 + t317;
t147 = -mrSges(4,2) * t298 - mrSges(4,3) * t168;
t146 = mrSges(4,1) * t298 - mrSges(4,3) * t167;
t116 = mrSges(5,1) * t200 - mrSges(5,3) * t162;
t115 = -mrSges(5,2) * t200 + mrSges(5,3) * t161;
t108 = qJD(4) * t178 + t177 * t271 + t267 * t301;
t107 = qJD(4) * t285 - t177 * t267 + t271 * t301;
t105 = mrSges(4,1) * t168 + mrSges(4,2) * t167;
t98 = pkin(4) * t328 + t132;
t97 = -pkin(4) * t178 + t140;
t91 = Ifges(4,5) * t298 - t344 + t346;
t90 = -t168 * Ifges(4,2) + Ifges(4,6) * t298 + t345;
t67 = mrSges(6,1) * t189 - mrSges(6,3) * t96;
t66 = -mrSges(6,2) * t189 + mrSges(6,3) * t299;
t59 = -mrSges(5,2) * t168 + mrSges(5,3) * t87;
t58 = mrSges(5,1) * t168 - mrSges(5,3) * t86;
t49 = -pkin(4) * t107 + t78;
t48 = -mrSges(5,1) * t87 + mrSges(5,2) * t86;
t40 = -pkin(4) * t87 + t63;
t38 = -qJD(5) * t118 + t107 * t270 - t108 * t266;
t37 = qJD(5) * t117 + t107 * t266 + t108 * t270;
t35 = t86 * Ifges(5,4) + t87 * Ifges(5,2) + t168 * Ifges(5,6);
t21 = -mrSges(6,2) * t168 + mrSges(6,3) * t31;
t20 = mrSges(6,1) * t168 - mrSges(6,3) * t30;
t17 = pkin(10) * t107 + t26;
t16 = t270 * t46 - t331;
t15 = -t266 * t46 - t330;
t12 = pkin(4) * t176 - pkin(10) * t108 + t27;
t9 = -mrSges(6,1) * t31 + mrSges(6,2) * t30;
t5 = -qJD(5) * t23 + t12 * t270 - t17 * t266;
t4 = qJD(5) * t22 + t12 * t266 + t17 * t270;
t1 = [t64 * t58 + t65 * t59 + t4 * t66 + t5 * t67 + t49 * t50 + t63 * (-mrSges(5,1) * t178 - mrSges(5,2) * t285) + t19 * (mrSges(5,1) * t231 + mrSges(5,3) * t285) + (-Ifges(5,4) * t285 + Ifges(5,2) * t178 + Ifges(5,6) * t231) * t393 + (-Ifges(5,1) * t285 + Ifges(5,4) * t178 + Ifges(5,5) * t231) * t394 + (-Ifges(5,5) * t285 + Ifges(6,5) * t118 + Ifges(5,6) * t178 + Ifges(6,6) * t117 + (Ifges(5,3) + Ifges(6,3)) * t231) * t379 - t285 * t396 + t205 * (Ifges(4,4) * t177 - Ifges(4,2) * t176) / 0.2e1 + t250 * (Ifges(4,5) * t177 - Ifges(4,6) * t176) / 0.2e1 - t326 * t357 - t307 * t326 / 0.2e1 - t168 * (Ifges(4,4) * t232 - Ifges(4,2) * t231 - Ifges(4,6) * t326) / 0.2e1 + t167 * (Ifges(4,1) * t232 - Ifges(4,4) * t231 - Ifges(4,5) * t326) / 0.2e1 + t216 * (-t265 * mrSges(3,2) + mrSges(3,3) * t326) + (-mrSges(3,1) * t265 + mrSges(4,1) * t231 + mrSges(4,2) * t232 + mrSges(3,3) * t327) * t217 + (t73 + t43) * t176 / 0.2e1 + m(6) * (t13 * t5 + t14 * t4 + t2 * t23 + t22 * t3 + t40 * t97 + t49 * t70) + m(5) * (t111 * t78 + t140 * t63 + t18 * t65 + t19 * t64 + t26 * t54 + t27 * t53) + m(4) * (t131 * t89 + t132 * t88 + t152 * t72 + t153 * t71 + t185 * t228 + t217 * t218) + m(3) * (t216 * t316 - t217 * t235 - t223 * t228 + t226 * t227) + ((Ifges(3,5) * t369 - t235 * mrSges(3,3) + (-0.2e1 * t385 + 0.3e1 / 0.2e1 * Ifges(3,4) * t273) * t264) * t273 + (Ifges(4,5) * t370 + Ifges(4,6) * t372 - Ifges(3,6) * t265 - t316 * mrSges(3,3) + (-0.2e1 * t386 - 0.3e1 / 0.2e1 * t354 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1) * t273) * t264) * t269) * t300 + t206 * (Ifges(4,1) * t177 - Ifges(4,4) * t176) / 0.2e1 + t411 * t231 / 0.2e1 + (t281 * t273 + (t412 + t278) * t269) * t314 + t37 * t45 / 0.2e1 + t318 * t228 + t22 * t20 + t23 * t21 + t70 * (-mrSges(6,1) * t38 + mrSges(6,2) * t37) + t97 * t9 + t78 * t101 + t107 * t74 / 0.2e1 + t108 * t75 / 0.2e1 + t111 * (-mrSges(5,1) * t107 + mrSges(5,2) * t108) + t26 * t115 + t27 * t116 + t40 * (-mrSges(6,1) * t117 + mrSges(6,2) * t118) + t140 * t48 + t152 * t146 + t153 * t147 + t88 * t169 + t89 * t170 + t13 * (mrSges(6,1) * t176 - mrSges(6,3) * t37) + t14 * (-mrSges(6,2) * t176 + mrSges(6,3) * t38) + t54 * (-mrSges(5,2) * t176 + mrSges(5,3) * t107) + t53 * (mrSges(5,1) * t176 - mrSges(5,3) * t108) - t176 * t127 / 0.2e1 + t177 * t128 / 0.2e1 + t178 * t35 / 0.2e1 + t326 * t358 + t248 * t369 + t91 * t370 + t90 * t372 + (Ifges(5,5) * t108 + Ifges(5,6) * t107 + Ifges(5,3) * t176) * t373 + (Ifges(6,5) * t37 + Ifges(6,6) * t38 + Ifges(6,3) * t176) * t376 + t185 * (mrSges(4,1) * t176 + mrSges(4,2) * t177) + t218 * t105 + t227 * t222 + t2 * (-mrSges(6,2) * t231 + mrSges(6,3) * t117) + t3 * (mrSges(6,1) * t231 - mrSges(6,3) * t118) + t18 * (-mrSges(5,2) * t231 + mrSges(5,3) * t178) + (Ifges(5,1) * t108 + Ifges(5,4) * t107 + Ifges(5,5) * t176) * t380 + (Ifges(5,4) * t108 + Ifges(5,2) * t107 + Ifges(5,6) * t176) * t382 + (Ifges(6,1) * t37 + Ifges(6,4) * t38 + Ifges(6,5) * t176) * t389 + (Ifges(6,4) * t37 + Ifges(6,2) * t38 + Ifges(6,6) * t176) * t391 + (Ifges(6,4) * t118 + Ifges(6,2) * t117 + Ifges(6,6) * t231) * t397 + (Ifges(6,1) * t118 + Ifges(6,4) * t117 + Ifges(6,5) * t231) * t398 + t118 * t399 + t117 * t400 + t38 * t413 + (-t131 * t177 - t132 * t176 - t231 * t71 - t232 * t72) * mrSges(4,3); (-t13 * t322 + t14 * t321 - t2 * t229 + t230 * t3) * mrSges(6,3) + t40 * (mrSges(6,1) * t229 - mrSges(6,2) * t230) + (-Ifges(6,5) * t230 - Ifges(6,6) * t229) * t379 + (-Ifges(6,4) * t230 - Ifges(6,2) * t229) * t397 + (-Ifges(6,1) * t230 - Ifges(6,4) * t229) * t398 + t417 * t67 + (t113 * t3 + t114 * t2 + t13 * t417 + t14 * t418 + t245 * t40 + t416 * t70) * m(6) + t418 * t66 + (-mrSges(6,1) * t321 + mrSges(6,2) * t322) * t70 + (t121 / 0.2e1 - t137 / 0.2e1) * t45 + (t122 / 0.2e1 - t136 / 0.2e1) * t44 + m(4) * (-pkin(2) * t217 + t361 * t71) + t416 * t50 + (t36 * t367 + t35 * t368 + t63 * t295 + t294 * t394 + t292 * t393 + t290 * t379 + t217 * mrSges(4,2) - t72 * mrSges(4,3) + t346 / 0.2e1 - t344 / 0.2e1 + t91 / 0.2e1 + (-t18 * t267 - t19 * t271) * mrSges(5,3) + (-m(4) * t72 + m(5) * t63 - t146 + t48) * pkin(8) + (t291 * t383 + t293 * t381 + t289 * t374 + t111 * t296 + t75 * t368 + t271 * t395 + (t267 * t53 - t271 * t54) * mrSges(5,3)) * qJD(4)) * t268 + m(5) * (t149 * t54 + t150 * t53 + t18 * t214 + t19 * t213) + ((qJD(2) * (Ifges(4,5) * t268 + Ifges(4,6) * t272) / 0.2e1 + (t386 + t354 / 0.2e1) * t315 + (t256 / 0.2e1 - qJD(2)) * Ifges(3,6) - t278) * t269 + ((t385 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t269) * t315 - t253 / 0.2e1 + (t375 - t337 / 0.2e1 + t279) * t272 + t274 * t268 - t281) * t273) * t315 + ((-t274 + (-m(4) * t132 - t169) * pkin(8)) * t268 + (t199 / 0.2e1 + t337 / 0.2e1 + (-m(4) * t131 + m(5) * t111 + t323) * pkin(8) - t415) * t272) * qJD(3) + t248 + (-t193 * t54 + t194 * t53) * mrSges(5,3) - t318 * t226 - m(5) * (t111 * t142 + t53 * t79 + t54 * t80) - m(4) * (t131 * t156 + t132 * t157 + t185 * t226) + (t71 * mrSges(4,3) + pkin(8) * t147 - t217 * mrSges(4,1) - t85 / 0.2e1 - t84 / 0.2e1 + t345 / 0.2e1 - t34 / 0.2e1 - t6 / 0.2e1 + t90 / 0.2e1 - t29 / 0.2e1 - t28 / 0.2e1 + (-Ifges(4,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t168 + t410) * t272 - pkin(2) * t105 + t113 * t20 + t114 * t21 - t142 * t101 - t157 * t169 - t156 * t170 + (Ifges(5,5) * t194 + Ifges(5,6) * t193) * t374 + (Ifges(6,5) * t121 + Ifges(6,6) * t122) * t376 + (Ifges(6,5) * t137 + Ifges(6,6) * t136) * t377 - t194 * t75 / 0.2e1 - t111 * (-mrSges(5,1) * t193 + mrSges(5,2) * t194) + t213 * t58 + t214 * t59 - t216 * mrSges(3,2) - t217 * mrSges(3,1) - t223 * t222 + (-t79 + t150) * t116 + t245 * t9 + (Ifges(5,1) * t194 + Ifges(5,4) * t193) * t381 + (Ifges(5,4) * t194 + Ifges(5,2) * t193) * t383 + (Ifges(6,1) * t121 + Ifges(6,4) * t122) * t389 + (Ifges(6,1) * t137 + Ifges(6,4) * t136) * t390 + (Ifges(6,4) * t121 + Ifges(6,2) * t122) * t391 + (Ifges(6,4) * t137 + Ifges(6,2) * t136) * t392 + t193 * t395 - t230 * t399 - t229 * t400 + t426 * t115; t404 * t66 + (t13 * t405 + t14 * t404 + t187 * t3 + t188 * t2 + t263 * t40 - t70 * t98) * m(6) + t405 * t67 + (Ifges(6,5) * t241 - Ifges(6,6) * t287 + t289) * t379 + (t13 * t319 - t14 * t320 - t2 * t287 - t241 * t3) * mrSges(6,3) + t40 * (mrSges(6,1) * t287 + mrSges(6,2) * t241) + (Ifges(6,4) * t241 - Ifges(6,2) * t287) * t397 + (Ifges(6,1) * t241 - Ifges(6,4) * t287) * t398 + t357 - t323 * t132 - t63 * t296 + t402 * mrSges(5,3) + (m(5) * t402 - t267 * t58 + t271 * t59 + (-m(5) * t288 - t267 * t115 - t271 * t116) * qJD(4)) * pkin(9) + t276 * t206 + (-pkin(3) * t63 - t111 * t132 - t53 * t68 - t54 * t69) * m(5) + (-Ifges(6,5) * t135 - Ifges(6,6) * t134) * t377 + (-Ifges(6,1) * t135 - Ifges(6,4) * t134) * t390 + (-Ifges(6,4) * t135 - Ifges(6,2) * t134) * t392 + (t375 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t206 + t415) * t205 - pkin(3) * t48 + t307 - t287 * t400 + (t135 / 0.2e1 - t173 / 0.2e1) * t45 + (t134 / 0.2e1 - t174 / 0.2e1) * t44 + (-Ifges(6,5) * t173 - Ifges(6,6) * t174) * t376 + (-Ifges(6,1) * t173 - Ifges(6,4) * t174) * t389 + (-Ifges(6,4) * t173 - Ifges(6,2) * t174) * t391 + (mrSges(6,1) * t320 - mrSges(6,2) * t319) * t70 - t358 + (t403 * t364 - t421) * qJD(4) - t98 * t50 - t69 * t115 - t68 * t116 - t131 * t169 + t35 * t367 + t187 * t20 + t188 * t21 + t263 * t9 + t291 * t393 + t293 * t394 + t267 * t396 + t241 * t399; -t16 * t66 - t15 * t67 + (t270 * t20 + t266 * t21 + m(6) * (t2 * t266 + t270 * t3) - t403 * t162 + (-t266 * t67 + t270 * t66 + m(6) * (-t13 * t266 + t14 * t270)) * qJD(5)) * pkin(4) - m(6) * (t13 * t15 + t14 * t16) + t96 * t413 - t410 + t411 + (-Ifges(5,2) * t162 + t159 + t75) * t383 + (t161 * t53 + t162 * t54) * mrSges(5,3) - t53 * t115 + t54 * t116 - t111 * (mrSges(5,1) * t162 + mrSges(5,2) * t161) + (Ifges(5,5) * t161 - Ifges(5,6) * t162) * t374 + t74 * t380 + (Ifges(5,1) * t161 - t352) * t381 + t420; -t13 * t66 + t14 * t67 + t44 * t389 + t419 + t420 + t6;];
tauc = t1(:);
