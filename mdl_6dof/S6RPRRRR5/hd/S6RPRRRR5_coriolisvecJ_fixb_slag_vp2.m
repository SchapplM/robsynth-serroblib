% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:34:56
% EndTime: 2018-11-23 16:35:07
% DurationCPUTime: 10.78s
% Computational Cost: add. (21357->596), mult. (55965->821), div. (0->0), fcn. (44600->10), ass. (0->278)
t268 = sin(qJ(6));
t269 = sin(qJ(5));
t272 = cos(qJ(6));
t273 = cos(qJ(5));
t247 = t268 * t273 + t269 * t272;
t400 = qJD(5) + qJD(6);
t209 = t400 * t247;
t266 = sin(pkin(11));
t271 = sin(qJ(3));
t267 = cos(pkin(11));
t275 = cos(qJ(3));
t327 = t267 * t275;
t242 = -t266 * t271 + t327;
t232 = t242 * qJD(1);
t243 = t266 * t275 + t267 * t271;
t233 = t243 * qJD(1);
t270 = sin(qJ(4));
t274 = cos(qJ(4));
t403 = t274 * t232 - t270 * t233;
t426 = t247 * t403;
t406 = t426 - t209;
t293 = t268 * t269 - t272 * t273;
t208 = t400 * t293;
t425 = t293 * t403;
t405 = t425 - t208;
t259 = pkin(3) * t270 + pkin(9);
t359 = -pkin(10) - t259;
t312 = qJD(5) * t359;
t352 = pkin(3) * qJD(4);
t319 = t274 * t352;
t424 = t403 * t269;
t436 = pkin(10) * t424;
t360 = pkin(7) + qJ(2);
t252 = t360 * t266;
t244 = qJD(1) * t252;
t253 = t360 * t267;
t245 = qJD(1) * t253;
t207 = -t244 * t271 + t245 * t275;
t183 = pkin(8) * t232 + t207;
t177 = t270 * t183;
t206 = -t275 * t244 - t245 * t271;
t182 = -pkin(8) * t233 + t206;
t122 = t182 * t274 - t177;
t295 = t232 * t270 + t274 * t233;
t160 = pkin(4) * t295 - pkin(9) * t403;
t131 = pkin(3) * t233 + t160;
t79 = t273 * t122 + t269 * t131;
t440 = t269 * t312 + t273 * t319 + t436 - t79;
t262 = t273 * pkin(10);
t418 = pkin(5) * t295 - t262 * t403;
t78 = -t122 * t269 + t273 * t131;
t439 = -t269 * t319 + t273 * t312 - t418 - t78;
t234 = t242 * qJD(3);
t226 = qJD(1) * t234;
t235 = t243 * qJD(3);
t227 = qJD(1) * t235;
t155 = qJD(4) * t403 + t226 * t274 - t227 * t270;
t156 = qJD(4) * t295 + t226 * t270 + t274 * t227;
t265 = qJD(3) + qJD(4);
t184 = t265 * t273 - t269 * t295;
t178 = t274 * t183;
t179 = qJD(3) * pkin(3) + t182;
t118 = t179 * t270 + t178;
t112 = pkin(9) * t265 + t118;
t317 = -pkin(2) * t267 - pkin(1);
t251 = qJD(1) * t317 + qJD(2);
t210 = -pkin(3) * t232 + t251;
t119 = -pkin(4) * t403 - pkin(9) * t295 + t210;
t69 = t112 * t273 + t119 * t269;
t55 = pkin(10) * t184 + t69;
t337 = t268 * t55;
t195 = qJD(5) - t403;
t185 = t265 * t269 + t273 * t295;
t68 = -t112 * t269 + t273 * t119;
t54 = -pkin(10) * t185 + t68;
t47 = pkin(5) * t195 + t54;
t20 = t272 * t47 - t337;
t335 = t272 * t55;
t21 = t268 * t47 + t335;
t107 = -qJD(5) * t185 - t155 * t269;
t321 = qJD(5) * t273;
t322 = qJD(5) * t269;
t117 = t179 * t274 - t177;
t258 = qJD(2) * t327;
t320 = qJD(1) * qJD(2);
t323 = qJD(3) * t275;
t175 = -t244 * t323 + qJD(1) * t258 + (-qJD(3) * t245 - t266 * t320) * t271;
t166 = -pkin(8) * t227 + t175;
t287 = t243 * qJD(2);
t176 = -qJD(1) * t287 - qJD(3) * t207;
t282 = -pkin(8) * t226 + t176;
t60 = qJD(4) * t117 + t274 * t166 + t270 * t282;
t364 = pkin(3) * t227;
t84 = pkin(4) * t156 - pkin(9) * t155 + t364;
t17 = -t112 * t322 + t119 * t321 + t269 * t84 + t273 * t60;
t12 = pkin(10) * t107 + t17;
t106 = qJD(5) * t184 + t155 * t273;
t18 = -qJD(5) * t69 - t269 * t60 + t273 * t84;
t7 = pkin(5) * t156 - pkin(10) * t106 + t18;
t3 = qJD(6) * t20 + t12 * t272 + t268 * t7;
t301 = Ifges(6,5) * t269 + Ifges(6,6) * t273;
t355 = Ifges(6,4) * t269;
t303 = Ifges(6,2) * t273 + t355;
t354 = Ifges(6,4) * t273;
t305 = Ifges(6,1) * t269 + t354;
t308 = mrSges(6,1) * t273 - mrSges(6,2) * t269;
t365 = t273 / 0.2e1;
t61 = qJD(4) * t118 + t166 * t270 - t274 * t282;
t37 = -pkin(5) * t107 + t61;
t193 = qJD(6) + t195;
t373 = t193 / 0.2e1;
t374 = -t193 / 0.2e1;
t378 = t156 / 0.2e1;
t127 = t184 * t268 + t185 * t272;
t379 = t127 / 0.2e1;
t380 = -t127 / 0.2e1;
t310 = t272 * t184 - t185 * t268;
t381 = t310 / 0.2e1;
t382 = -t310 / 0.2e1;
t383 = t107 / 0.2e1;
t384 = t106 / 0.2e1;
t120 = Ifges(7,4) * t310;
t65 = t127 * Ifges(7,1) + t193 * Ifges(7,5) + t120;
t388 = t65 / 0.2e1;
t389 = -t65 / 0.2e1;
t353 = Ifges(7,4) * t127;
t64 = Ifges(7,2) * t310 + t193 * Ifges(7,6) + t353;
t390 = t64 / 0.2e1;
t391 = -t64 / 0.2e1;
t41 = -qJD(6) * t127 - t106 * t268 + t107 * t272;
t392 = t41 / 0.2e1;
t40 = qJD(6) * t310 + t106 * t272 + t107 * t268;
t393 = t40 / 0.2e1;
t394 = Ifges(7,1) * t393 + Ifges(7,4) * t392 + Ifges(7,5) * t378;
t395 = Ifges(7,4) * t393 + Ifges(7,2) * t392 + Ifges(7,6) * t378;
t4 = -qJD(6) * t21 - t12 * t268 + t272 * t7;
t356 = Ifges(6,4) * t185;
t102 = t184 * Ifges(6,2) + t195 * Ifges(6,6) + t356;
t181 = Ifges(6,4) * t184;
t103 = t185 * Ifges(6,1) + t195 * Ifges(6,5) + t181;
t111 = -pkin(4) * t265 - t117;
t302 = Ifges(6,5) * t273 - Ifges(6,6) * t269;
t304 = -Ifges(6,2) * t269 + t354;
t306 = Ifges(6,1) * t273 - t355;
t307 = mrSges(6,1) * t269 + mrSges(6,2) * t273;
t366 = -t269 / 0.2e1;
t375 = t185 / 0.2e1;
t417 = t111 * t307 + t195 * t302 / 0.2e1 + t306 * t375 + t184 * t304 / 0.2e1 + t103 * t365 + t102 * t366;
t421 = t17 * t273 - t18 * t269;
t44 = t106 * Ifges(6,4) + t107 * Ifges(6,2) + t156 * Ifges(6,6);
t45 = t106 * Ifges(6,1) + t107 * Ifges(6,4) + t156 * Ifges(6,5);
t92 = -pkin(5) * t184 + t111;
t438 = (Ifges(7,4) * t247 - Ifges(7,2) * t293) * t392 + (Ifges(7,1) * t247 - Ifges(7,4) * t293) * t393 + t37 * (mrSges(7,1) * t293 + mrSges(7,2) * t247) - t293 * t395 + (Ifges(7,5) * t247 - Ifges(7,6) * t293 + t301) * t378 - t426 * t391 + (-Ifges(7,4) * t425 - Ifges(7,2) * t426) * t382 + (-Ifges(7,5) * t425 - Ifges(7,6) * t426) * t374 + (-Ifges(7,1) * t425 - Ifges(7,4) * t426) * t380 - t425 * t389 + t303 * t383 + t305 * t384 + t417 * qJD(5) + (-t308 - mrSges(5,1)) * t61 + t44 * t365 + (-t20 * t405 + t21 * t406 - t4 * t247 - t293 * t3) * mrSges(7,3) - t208 * t388 - t209 * t390 + t247 * t394 + t421 * mrSges(6,3) + (-Ifges(7,5) * t208 - Ifges(7,6) * t209) * t373 + (-Ifges(7,1) * t208 - Ifges(7,4) * t209) * t379 + (-Ifges(7,4) * t208 - Ifges(7,2) * t209) * t381 + Ifges(5,5) * t155 - Ifges(5,6) * t156 - t60 * mrSges(5,2) + (-t406 * mrSges(7,1) + mrSges(7,2) * t405) * t92 + t269 * t45 / 0.2e1;
t437 = pkin(5) * t424;
t387 = -pkin(10) - pkin(9);
t316 = qJD(5) * t387;
t81 = t273 * t117 + t269 * t160;
t434 = t269 * t316 + t436 - t81;
t80 = -t117 * t269 + t273 * t160;
t433 = t273 * t316 - t418 - t80;
t151 = Ifges(7,3) * t156;
t431 = t151 + (Ifges(7,5) * t310 - Ifges(7,6) * t127) * t374 + (t127 * t21 + t20 * t310) * mrSges(7,3) - t92 * (mrSges(7,1) * t127 + mrSges(7,2) * t310) + (Ifges(7,1) * t310 - t353) * t380;
t237 = t359 * t269;
t238 = t259 * t273 + t262;
t203 = t237 * t268 + t238 * t272;
t429 = -qJD(6) * t203 - t268 * t440 + t439 * t272;
t202 = t237 * t272 - t238 * t268;
t428 = qJD(6) * t202 + t439 * t268 + t272 * t440;
t53 = -mrSges(6,1) * t107 + mrSges(6,2) * t106;
t427 = -m(6) * t61 - t53;
t121 = t182 * t270 + t178;
t318 = pkin(5) * t322;
t423 = t270 * t352 - t121 + t318 - t437;
t299 = t269 * t69 + t273 * t68;
t291 = t299 * mrSges(6,3);
t194 = Ifges(5,4) * t403;
t340 = t295 * Ifges(5,1);
t398 = t194 / 0.2e1 + t340 / 0.2e1;
t422 = t210 * mrSges(5,2) + t265 * Ifges(5,5) - t291 + t398 + t417;
t420 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,5) * t40 + Ifges(7,6) * t41;
t419 = -Ifges(7,2) * t127 + t120;
t415 = -t403 * Ifges(5,2) / 0.2e1;
t255 = t387 * t269;
t256 = pkin(9) * t273 + t262;
t214 = t255 * t268 + t256 * t272;
t408 = -qJD(6) * t214 - t268 * t434 + t272 * t433;
t213 = t255 * t272 - t256 * t268;
t407 = qJD(6) * t213 + t268 * t433 + t272 * t434;
t205 = t242 * t270 + t243 * t274;
t162 = t293 * t205;
t211 = -t275 * t252 - t253 * t271;
t191 = -pkin(8) * t243 + t211;
t212 = -t271 * t252 + t275 * t253;
t192 = pkin(8) * t242 + t212;
t142 = t191 * t270 + t192 * t274;
t134 = t273 * t142;
t220 = -pkin(3) * t242 + t317;
t294 = t274 * t242 - t243 * t270;
t143 = -pkin(4) * t294 - pkin(9) * t205 + t220;
t86 = t269 * t143 + t134;
t404 = t274 * t191 - t192 * t270;
t402 = t117 * mrSges(5,3) - t194 / 0.2e1 - t422;
t401 = -t269 * t68 + t273 * t69;
t397 = (m(3) * qJ(2) + mrSges(3,3)) * (t266 ^ 2 + t267 ^ 2);
t396 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t106 + Ifges(6,6) * t107 + t420;
t377 = -t184 / 0.2e1;
t376 = -t185 / 0.2e1;
t371 = -t195 / 0.2e1;
t369 = -t233 / 0.2e1;
t368 = t234 / 0.2e1;
t367 = -t235 / 0.2e1;
t363 = pkin(3) * t235;
t362 = pkin(3) * t274;
t358 = Ifges(4,4) * t233;
t351 = t118 * mrSges(5,3);
t348 = t404 * t61;
t328 = t205 * t269;
t326 = -mrSges(5,1) * t265 - mrSges(6,1) * t184 + mrSges(6,2) * t185 + mrSges(5,3) * t295;
t261 = -pkin(5) * t273 - pkin(4);
t188 = -t252 * t323 + t258 + (-qJD(2) * t266 - qJD(3) * t253) * t271;
t170 = -pkin(8) * t235 + t188;
t189 = -qJD(3) * t212 - t287;
t171 = -pkin(8) * t234 + t189;
t75 = qJD(4) * t404 + t170 * t274 + t171 * t270;
t164 = qJD(4) * t294 + t234 * t274 - t235 * t270;
t165 = qJD(4) * t205 + t234 * t270 + t274 * t235;
t91 = pkin(4) * t165 - pkin(9) * t164 + t363;
t313 = -t269 * t75 + t273 * t91;
t311 = t156 * mrSges(5,1) + t155 * mrSges(5,2);
t85 = -t142 * t269 + t273 * t143;
t300 = -t17 * t269 - t18 * t273;
t62 = -pkin(5) * t294 - t205 * t262 + t85;
t74 = -pkin(10) * t328 + t86;
t30 = -t268 * t74 + t272 * t62;
t31 = t268 * t62 + t272 * t74;
t70 = mrSges(6,1) * t156 - mrSges(6,3) * t106;
t71 = -mrSges(6,2) * t156 + mrSges(6,3) * t107;
t297 = -t269 * t70 + t273 * t71;
t135 = -mrSges(6,2) * t195 + mrSges(6,3) * t184;
t136 = mrSges(6,1) * t195 - mrSges(6,3) * t185;
t296 = t273 * t135 - t269 * t136;
t186 = -mrSges(5,2) * t265 + mrSges(5,3) * t403;
t292 = t186 + t296;
t290 = t164 * t269 + t205 * t321;
t22 = -t142 * t322 + t143 * t321 + t269 * t91 + t273 * t75;
t76 = qJD(4) * t142 + t170 * t270 - t274 * t171;
t283 = m(6) * (-qJD(5) * t299 + t421);
t280 = t210 * mrSges(5,1) + t68 * mrSges(6,1) + t20 * mrSges(7,1) - t69 * mrSges(6,2) - t21 * mrSges(7,2) - Ifges(5,4) * t295 + t185 * Ifges(6,5) + t127 * Ifges(7,5) - t265 * Ifges(5,6) + t184 * Ifges(6,6) + t310 * Ifges(7,6) + t195 * Ifges(6,3) + t193 * Ifges(7,3) + t415;
t278 = t280 + t415;
t254 = t261 - t362;
t228 = Ifges(4,4) * t232;
t221 = t226 * mrSges(4,2);
t219 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t233;
t218 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t232;
t197 = t233 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t228;
t196 = t232 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t358;
t161 = t247 * t205;
t159 = -mrSges(5,1) * t403 + mrSges(5,2) * t295;
t152 = Ifges(6,3) * t156;
t108 = pkin(5) * t328 - t404;
t97 = mrSges(7,1) * t193 - mrSges(7,3) * t127;
t96 = -mrSges(7,2) * t193 + mrSges(7,3) * t310;
t93 = t118 + t437;
t77 = -mrSges(7,1) * t310 + mrSges(7,2) * t127;
t49 = t162 * t400 - t247 * t164;
t48 = -t164 * t293 - t205 * t209;
t46 = pkin(5) * t290 + t76;
t33 = -mrSges(7,2) * t156 + mrSges(7,3) * t41;
t32 = mrSges(7,1) * t156 - mrSges(7,3) * t40;
t25 = t272 * t54 - t337;
t24 = -t268 * t54 - t335;
t23 = -qJD(5) * t86 + t313;
t15 = -pkin(10) * t290 + t22;
t14 = -t164 * t262 + pkin(5) * t165 + (-t134 + (pkin(10) * t205 - t143) * t269) * qJD(5) + t313;
t13 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t6 = -qJD(6) * t31 + t14 * t272 - t15 * t268;
t5 = qJD(6) * t30 + t14 * t268 + t15 * t272;
t1 = [(-t117 * t164 - t118 * t165 - t142 * t156 - t155 * t404) * mrSges(5,3) - t404 * t53 + qJD(3) * (Ifges(4,5) * t234 - Ifges(4,6) * t235) / 0.2e1 + t251 * (mrSges(4,1) * t235 + mrSges(4,2) * t234) + (-t161 * t3 + t162 * t4 - t20 * t48 + t21 * t49) * mrSges(7,3) + (-Ifges(7,5) * t162 - Ifges(7,6) * t161) * t378 + (-Ifges(7,4) * t162 - Ifges(7,2) * t161) * t392 + (-Ifges(7,1) * t162 - Ifges(7,4) * t161) * t393 + t37 * (mrSges(7,1) * t161 - mrSges(7,2) * t162) + (t45 * t365 + t44 * t366 + mrSges(5,2) * t364 + Ifges(5,1) * t155 - Ifges(5,4) * t156 + t306 * t384 + t304 * t383 + t302 * t378 + (mrSges(5,3) + t307) * t61 + t300 * mrSges(6,3) + (t103 * t366 - t273 * t102 / 0.2e1 + t111 * t308 + t303 * t377 + t305 * t376 + t301 * t371 - t401 * mrSges(6,3)) * qJD(5)) * t205 + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t373 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t379 + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t381 - (mrSges(5,1) * t364 - t60 * mrSges(5,3) + t151 / 0.2e1 + t152 / 0.2e1 - Ifges(5,4) * t155 + (Ifges(7,3) / 0.2e1 + Ifges(5,2) + Ifges(6,3) / 0.2e1) * t156 + t396) * t294 + 0.2e1 * t397 * t320 + t196 * t367 + t197 * t368 + t278 * t165 + t48 * t388 + t49 * t390 - t162 * t394 - t161 * t395 + (t398 + t422) * t164 + (t175 * t242 - t176 * t243 - t206 * t234 - t207 * t235 - t211 * t226 - t212 * t227) * mrSges(4,3) + t317 * (t227 * mrSges(4,1) + t221) + (-t242 * t227 + t232 * t367) * Ifges(4,2) + (t242 * t226 - t243 * t227 + t232 * t368 + t233 * t367) * Ifges(4,4) + m(5) * (-t117 * t76 + t118 * t75 - t348 + t142 * t60 + (t210 * t235 + t220 * t227) * pkin(3)) + m(6) * (t111 * t76 + t17 * t86 + t18 * t85 + t22 * t69 + t23 * t68 - t348) + t30 * t32 + t31 * t33 + m(4) * (t175 * t212 + t176 * t211 + t188 * t207 + t189 * t206) + m(7) * (t108 * t37 + t20 * t6 + t21 * t5 + t3 * t31 + t30 * t4 + t46 * t92) + t220 * t311 + t46 * t77 + t85 * t70 + t86 * t71 + t92 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t5 * t96 + t6 * t97 + t326 * t76 + t108 * t13 + t22 * t135 + t23 * t136 + t159 * t363 + t75 * t186 + (t243 * t226 + t233 * t368) * Ifges(4,1) + t188 * t218 + t189 * t219; (-t77 - t326) * t295 + t296 * qJD(5) + t221 - t292 * t403 - m(4) * (-t206 * t233 + t207 * t232) - m(5) * (-t117 * t295 + t118 * t403) + t406 * t97 + t405 * t96 - (-m(5) * pkin(3) - mrSges(4,1)) * t227 + t311 - t232 * t218 + t233 * t219 - t293 * t32 + t247 * t33 + t269 * t71 + t273 * t70 - t397 * qJD(1) ^ 2 + (t20 * t406 + t21 * t405 + t247 * t3 - t293 * t4 - t295 * t92) * m(7) + (-t111 * t295 + t195 * t401 - t300) * m(6); -(-Ifges(4,2) * t233 + t197 + t228) * t232 / 0.2e1 + (-t340 / 0.2e1 + t402) * t403 - t291 * qJD(5) + (t297 + t283 + (-t135 * t269 - t136 * t273) * qJD(5)) * t259 + (-t278 + t351) * t295 - m(5) * (-t117 * t121 + t118 * t122) + (Ifges(4,1) * t232 - t358) * t369 + t438 + t423 * t77 - t427 * (-pkin(4) - t362) + t428 * t96 + t429 * t97 + (t20 * t429 + t202 * t4 + t203 * t3 + t21 * t428 + t254 * t37 + t423 * t92) * m(7) - m(6) * (t111 * t121 + t68 * t78 + t69 * t79) + (-t233 * t159 + (-t155 * t274 - t156 * t270) * mrSges(5,3) + (t326 * t270 + t292 * t274 + m(6) * (t111 * t270 + t274 * t401)) * qJD(4) + (t270 * t60 - t274 * t61 + 0.2e1 * t210 * t369 + (-t117 * t270 + t118 * t274) * qJD(4)) * m(5)) * pkin(3) + (t206 * t232 + t207 * t233) * mrSges(4,3) - t326 * t121 - t79 * t135 - t78 * t136 - t175 * mrSges(4,2) + t176 * mrSges(4,1) - t122 * t186 + t202 * t32 + t203 * t33 - t206 * t218 + t207 * t219 + Ifges(4,5) * t226 - Ifges(4,6) * t227 + t233 * t196 / 0.2e1 - qJD(3) * (Ifges(4,5) * t232 - Ifges(4,6) * t233) / 0.2e1 - t251 * (mrSges(4,1) * t233 + mrSges(4,2) * t232) + t254 * t13; ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t295 + t402) * t403 + t297 * pkin(9) - m(6) * (t111 * t118 + t68 * t80 + t69 * t81) + (-t280 + t351) * t295 + pkin(9) * t283 + t427 * pkin(4) + ((-t68 * mrSges(6,3) - pkin(9) * t136) * t273 + (-t69 * mrSges(6,3) + pkin(5) * t77 - pkin(9) * t135) * t269) * qJD(5) + t407 * t96 + t408 * t97 + (t213 * t4 + t214 * t3 + t261 * t37 + (t318 - t93) * t92 + t407 * t21 + t408 * t20) * m(7) - t93 * t77 - t326 * t118 - t81 * t135 - t80 * t136 - t117 * t186 + t213 * t32 + t214 * t33 + t261 * t13 + t438; (t184 * t68 + t185 * t69) * mrSges(6,3) - t127 * t391 + t310 * t389 + (Ifges(6,5) * t184 - Ifges(6,6) * t185) * t371 + t102 * t375 + (Ifges(6,1) * t184 - t356) * t376 + t419 * t382 + t396 + t152 + (-Ifges(6,2) * t185 + t103 + t181) * t377 - m(7) * (t20 * t24 + t21 * t25) + (-t185 * t77 + t268 * t33 + t272 * t32 + (-t268 * t97 + t272 * t96) * qJD(6) + (-t185 * t92 + t268 * t3 + t272 * t4 + (-t20 * t268 + t21 * t272) * qJD(6)) * m(7)) * pkin(5) - t25 * t96 - t24 * t97 - t68 * t135 + t69 * t136 - t111 * (mrSges(6,1) * t185 + mrSges(6,2) * t184) + t431; t64 * t379 - t20 * t96 + t21 * t97 + (t419 + t65) * t382 + t420 + t431;];
tauc  = t1(:);
