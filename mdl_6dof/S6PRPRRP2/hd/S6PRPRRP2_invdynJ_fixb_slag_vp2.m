% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:59
% EndTime: 2019-03-08 20:00:25
% DurationCPUTime: 15.77s
% Computational Cost: add. (5052->601), mult. (11725->799), div. (0->0), fcn. (9160->12), ass. (0->275)
t434 = Ifges(6,1) + Ifges(7,1);
t433 = Ifges(7,4) + Ifges(6,5);
t432 = Ifges(6,6) - Ifges(7,6);
t199 = sin(pkin(6));
t197 = sin(pkin(11));
t200 = cos(pkin(11));
t205 = sin(qJ(2));
t208 = cos(qJ(2));
t234 = t197 * t208 + t200 * t205;
t143 = t234 * t199;
t134 = qJD(1) * t143;
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t260 = pkin(4) * t204 - pkin(9) * t207;
t437 = t260 * qJD(4) - t134;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t307 = t206 * qJD(4);
t315 = qJD(2) * t204;
t166 = t203 * t315 - t307;
t374 = t166 / 0.2e1;
t314 = qJD(4) * t203;
t167 = t206 * t315 + t314;
t436 = -t167 / 0.2e1;
t306 = t207 * qJD(2);
t189 = qJD(5) - t306;
t370 = t189 / 0.2e1;
t435 = qJD(4) / 0.2e1;
t431 = -Ifges(6,3) - Ifges(7,2);
t395 = -t203 * t432 + t206 * t433;
t347 = Ifges(7,5) * t203;
t350 = Ifges(6,4) * t203;
t393 = t206 * t434 + t347 - t350;
t305 = qJD(2) * qJD(4);
t175 = qJDD(2) * t204 + t207 * t305;
t311 = qJD(5) * t166;
t102 = qJDD(4) * t203 + t175 * t206 - t311;
t379 = t102 / 0.2e1;
t103 = qJD(5) * t167 - t206 * qJDD(4) + t175 * t203;
t378 = -t103 / 0.2e1;
t174 = qJDD(2) * t207 - t204 * t305;
t163 = qJDD(5) - t174;
t376 = t163 / 0.2e1;
t330 = t199 * t205;
t288 = qJD(1) * t330;
t178 = t197 * t288;
t317 = qJD(1) * t208;
t287 = t199 * t317;
t137 = t200 * t287 - t178;
t366 = pkin(4) * t207;
t233 = -pkin(9) * t204 - pkin(3) - t366;
t367 = pkin(2) * t200;
t160 = t233 - t367;
t308 = qJD(5) * t206;
t322 = t206 * t207;
t430 = -t137 * t322 + t160 * t308 + t203 * t437;
t324 = t203 * t207;
t192 = pkin(2) * t197 + pkin(8);
t402 = t203 * t160 + t192 * t322;
t429 = -qJD(5) * t402 + t137 * t324 + t206 * t437;
t57 = -mrSges(6,2) * t163 - mrSges(6,3) * t103;
t58 = -mrSges(7,2) * t103 + mrSges(7,3) * t163;
t356 = t57 + t58;
t428 = t356 * t206;
t55 = mrSges(6,1) * t163 - mrSges(6,3) * t102;
t56 = -t163 * mrSges(7,1) + t102 * mrSges(7,2);
t357 = -t55 + t56;
t427 = t357 * t203;
t176 = qJD(2) * pkin(2) + t287;
t120 = t197 * t176 + t200 * t288;
t117 = qJD(2) * pkin(8) + t120;
t202 = cos(pkin(6));
t187 = qJD(1) * t202 + qJD(3);
t334 = t187 * t204;
t77 = t117 * t207 + t334;
t68 = qJD(4) * pkin(9) + t77;
t119 = t176 * t200 - t178;
t84 = qJD(2) * t233 - t119;
t19 = -t203 * t68 + t206 * t84;
t20 = t203 * t84 + t206 * t68;
t266 = qJD(2) * t288;
t328 = t199 * t208;
t148 = qJDD(1) * t328 - t266;
t339 = qJDD(2) * pkin(2);
t140 = t148 + t339;
t282 = qJD(2) * t317;
t149 = (qJDD(1) * t205 + t282) * t199;
t69 = t140 * t200 - t197 * t149;
t62 = -qJDD(2) * pkin(3) - t69;
t30 = -pkin(4) * t174 - pkin(9) * t175 + t62;
t310 = qJD(5) * t203;
t186 = qJDD(1) * t202 + qJDD(3);
t312 = qJD(4) * t207;
t313 = qJD(4) * t204;
t70 = t197 * t140 + t200 * t149;
t63 = qJDD(2) * pkin(8) + t70;
t11 = -t117 * t313 + t204 * t186 + t187 * t312 + t207 * t63;
t9 = qJDD(4) * pkin(9) + t11;
t3 = t203 * t30 + t206 * t9 + t84 * t308 - t310 * t68;
t4 = -qJD(5) * t20 - t203 * t9 + t206 * t30;
t258 = -t203 * t4 + t206 * t3;
t426 = -t19 * t308 - t20 * t310 + t258;
t404 = qJD(6) - t19;
t13 = -pkin(5) * t189 + t404;
t14 = qJ(6) * t189 + t20;
t1 = qJ(6) * t163 + qJD(6) * t189 + t3;
t2 = -pkin(5) * t163 + qJDD(6) - t4;
t259 = t1 * t206 + t2 * t203;
t425 = t13 * t308 - t14 * t310 + t259;
t198 = sin(pkin(10));
t201 = cos(pkin(10));
t326 = t202 * t208;
t424 = -t198 * t205 + t201 * t326;
t194 = Ifges(5,4) * t306;
t162 = Ifges(6,4) * t166;
t348 = Ifges(7,5) * t166;
t409 = t167 * t434 + t189 * t433 - t162 + t348;
t161 = Ifges(7,5) * t167;
t91 = t189 * Ifges(7,6) + t166 * Ifges(7,3) + t161;
t423 = Ifges(5,1) * t315 + Ifges(5,5) * qJD(4) + t203 * t91 + t206 * t409 + t194;
t353 = Ifges(5,4) * t204;
t408 = Ifges(5,2) * t207;
t246 = t353 + t408;
t422 = Ifges(5,6) * t435 + qJD(2) * t246 / 0.2e1 + t431 * t370 + t433 * t436 + t432 * t374;
t421 = -m(6) - m(7);
t420 = t174 / 0.2e1;
t419 = t175 / 0.2e1;
t418 = mrSges(5,3) - mrSges(4,2);
t417 = mrSges(6,3) + mrSges(7,2);
t415 = t433 * t163 + (-Ifges(6,4) + Ifges(7,5)) * t103 + t434 * t102;
t414 = (-t192 * t310 - qJD(6)) * t207 + (-t192 * t206 + qJ(6)) * t313 + t430;
t271 = t192 * t203 + pkin(5);
t413 = -t271 * t313 - t429;
t412 = (-t204 * t307 - t207 * t310) * t192 + t430;
t286 = t192 * t313;
t411 = t203 * t286 + t429;
t235 = pkin(5) * t203 - qJ(6) * t206;
t407 = qJD(5) * t235 - qJD(6) * t203 - t334 - (qJD(2) * t235 + t117) * t207;
t255 = mrSges(5,1) * t207 - mrSges(5,2) * t204;
t406 = -t255 - mrSges(4,1);
t33 = mrSges(6,1) * t103 + mrSges(6,2) * t102;
t405 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t175 + t33;
t298 = mrSges(5,3) * t315;
t401 = -qJD(4) * mrSges(5,1) + mrSges(6,1) * t166 + mrSges(6,2) * t167 + t298;
t400 = -t204 * t431 + t207 * t395;
t399 = t204 * t433 + t207 * t393;
t321 = t208 * t200;
t142 = t197 * t330 - t199 * t321;
t398 = (mrSges(3,1) * t208 - mrSges(3,2) * t205) * t199 - mrSges(4,1) * t142 - mrSges(4,2) * t143;
t251 = mrSges(7,1) * t203 - mrSges(7,3) * t206;
t253 = mrSges(6,1) * t203 + mrSges(6,2) * t206;
t76 = -t204 * t117 + t187 * t207;
t67 = -qJD(4) * pkin(4) - t76;
t31 = pkin(5) * t166 - qJ(6) * t167 + t67;
t397 = t31 * t251 + t67 * t253;
t396 = t203 * t433 + t206 * t432;
t346 = Ifges(7,5) * t206;
t349 = Ifges(6,4) * t206;
t394 = t203 * t434 - t346 + t349;
t390 = t102 * t433 - t103 * t432 - t163 * t431;
t12 = -t117 * t312 + t186 * t207 - t187 * t313 - t204 * t63;
t389 = t11 * t207 - t12 * t204;
t388 = -t102 * Ifges(7,5) / 0.2e1 - t163 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t379 + Ifges(6,6) * t376 + (Ifges(7,3) + Ifges(6,2)) * t378;
t236 = pkin(5) * t206 + qJ(6) * t203;
t252 = -t206 * mrSges(7,1) - t203 * mrSges(7,3);
t254 = mrSges(6,1) * t206 - mrSges(6,2) * t203;
t387 = -m(7) * t236 - mrSges(5,1) + t252 - t254;
t268 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t267 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t384 = -m(6) * t67 - t401;
t115 = t143 * t207 + t202 * t204;
t331 = t199 * t204;
t145 = t234 * t202;
t164 = t197 * t205 - t321;
t87 = t145 * t201 - t164 * t198;
t52 = -t201 * t331 + t207 * t87;
t88 = t198 * t145 + t164 * t201;
t54 = t198 * t331 - t207 * t88;
t383 = -g(1) * t54 - g(2) * t52 - g(3) * t115;
t382 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t381 = qJD(2) ^ 2;
t351 = Ifges(6,4) * t167;
t94 = -t166 * Ifges(6,2) + t189 * Ifges(6,6) + t351;
t380 = -t94 / 0.2e1;
t377 = t103 / 0.2e1;
t375 = -t166 / 0.2e1;
t372 = t167 / 0.2e1;
t369 = t203 / 0.2e1;
t359 = -qJD(2) / 0.2e1;
t358 = qJD(5) / 0.2e1;
t355 = mrSges(6,3) * t166;
t354 = mrSges(6,3) * t167;
t352 = Ifges(5,4) * t207;
t10 = -qJDD(4) * pkin(4) - t12;
t345 = t10 * t204;
t144 = t164 * t202;
t86 = -t144 * t201 - t198 * t234;
t341 = t204 * t86;
t89 = t144 * t198 - t201 * t234;
t340 = t204 * t89;
t170 = t260 * qJD(2);
t39 = t203 * t170 + t206 * t76;
t337 = t142 * t204;
t336 = t160 * t206;
t329 = t199 * t207;
t327 = t202 * t205;
t323 = t204 * t206;
t125 = -mrSges(6,2) * t189 - t355;
t128 = -mrSges(7,2) * t166 + mrSges(7,3) * t189;
t320 = t125 + t128;
t126 = mrSges(6,1) * t189 - t354;
t127 = -mrSges(7,1) * t189 + mrSges(7,2) * t167;
t319 = -t126 + t127;
t309 = qJD(5) * t204;
t188 = pkin(2) * t328;
t32 = mrSges(7,1) * t103 - mrSges(7,3) * t102;
t304 = t32 + t405;
t303 = pkin(9) * t310;
t302 = pkin(9) * t308;
t299 = m(4) + m(5) - t421;
t297 = mrSges(5,3) * t306;
t112 = mrSges(7,1) * t166 - mrSges(7,3) * t167;
t291 = t112 + t401;
t284 = t203 * t312;
t273 = t308 / 0.2e1;
t269 = t305 / 0.2e1;
t262 = t424 * pkin(2);
t261 = -t142 * pkin(3) + pkin(8) * t143 + t188;
t245 = -Ifges(6,2) * t203 + t349;
t244 = Ifges(6,2) * t206 + t350;
t241 = Ifges(5,5) * t207 - Ifges(5,6) * t204;
t238 = Ifges(7,3) * t203 + t346;
t237 = -Ifges(7,3) * t206 + t347;
t38 = t170 * t206 - t203 * t76;
t46 = t115 * t206 + t142 * t203;
t45 = t115 * t203 - t142 * t206;
t114 = t143 * t204 - t202 * t207;
t232 = t192 + t235;
t230 = -t198 * t326 - t201 * t205;
t116 = -qJD(2) * pkin(3) - t119;
t226 = t116 * (mrSges(5,1) * t204 + mrSges(5,2) * t207);
t225 = t204 * (Ifges(5,1) * t207 - t353);
t223 = t86 * pkin(3) + pkin(8) * t87 + t262;
t222 = t230 * pkin(2);
t221 = -t203 * t309 + t207 * t307;
t220 = t204 * t308 + t284;
t215 = Ifges(6,6) * t204 + t207 * t245;
t214 = Ifges(7,6) * t204 + t207 * t238;
t212 = t89 * pkin(3) - pkin(8) * t88 + t222;
t193 = -pkin(3) - t367;
t182 = -qJD(4) * mrSges(5,2) + t297;
t179 = -pkin(4) - t236;
t168 = t255 * qJD(2);
t150 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t174;
t136 = t164 * t199 * qJD(2);
t135 = qJD(2) * t143;
t124 = t232 * t204;
t118 = -mrSges(5,1) * t174 + mrSges(5,2) * t175;
t111 = pkin(5) * t167 + qJ(6) * t166;
t109 = -t192 * t324 + t336;
t105 = t207 * t271 - t336;
t104 = -qJ(6) * t207 + t402;
t61 = (qJD(5) * t236 - qJD(6) * t206) * t204 + t232 * t312;
t53 = t198 * t329 + t204 * t88;
t51 = -t201 * t329 - t204 * t87;
t44 = -qJD(4) * t114 - t136 * t207;
t43 = qJD(4) * t115 - t136 * t204;
t36 = -pkin(5) * t315 - t38;
t35 = qJ(6) * t315 + t39;
t17 = t203 * t54 + t89 * t206;
t15 = t203 * t52 + t86 * t206;
t7 = -qJD(5) * t45 + t135 * t203 + t206 * t44;
t6 = qJD(5) * t46 - t135 * t206 + t203 * t44;
t5 = pkin(5) * t103 - qJ(6) * t102 - qJD(6) * t167 + t10;
t8 = [m(2) * qJDD(1) + t115 * t150 + t142 * t118 - t135 * t168 + t44 * t182 + t320 * t7 + t319 * t6 + t356 * t46 + t357 * t45 + (-mrSges(3,1) * t205 - mrSges(3,2) * t208) * t381 * t199 + (-mrSges(4,1) * t135 + mrSges(4,2) * t136) * qJD(2) + t291 * t43 + t304 * t114 + t398 * qJDD(2) + (-m(2) - m(3) - t299) * g(3) + m(7) * (t1 * t46 + t114 * t5 + t13 * t6 + t14 * t7 + t2 * t45 + t31 * t43) + m(6) * (t10 * t114 - t19 * t6 + t20 * t7 + t3 * t46 - t4 * t45 + t43 * t67) + m(5) * (t11 * t115 - t114 * t12 + t116 * t135 + t142 * t62 - t43 * t76 + t44 * t77) + m(4) * (-t119 * t135 - t120 * t136 - t142 * t69 + t143 * t70 + t186 * t202) + m(3) * (qJDD(1) * t202 ^ 2 + (t148 * t208 + t149 * t205) * t199); t423 * t312 / 0.2e1 + (-t424 * mrSges(3,1) - (-t198 * t208 - t201 * t327) * mrSges(3,2) - m(5) * t223 - m(4) * t262 + t406 * t86 - t418 * t87 - t417 * t341 + t421 * (pkin(9) * t341 + t86 * t366 + t223) - t268 * (t203 * t87 + t322 * t86) + t267 * (-t206 * t87 + t324 * t86)) * g(2) + (t19 * mrSges(6,1) - t13 * mrSges(7,1) - t20 * mrSges(6,2) - t77 * mrSges(5,3) + t14 * mrSges(7,3) - t422) * t313 + (t109 * t4 + t411 * t19 + t412 * t20 + t3 * t402) * m(6) + t402 * t57 + (t266 + t148) * mrSges(3,1) + (qJD(2) * t134 + t200 * t339 + t69) * mrSges(4,1) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t13 * t221 - t14 * t220 + t2 * t323) * mrSges(7,2) + t253 * t345 + t352 * t419 + t246 * t420 + (t199 * t282 - t149) * mrSges(3,2) + (t401 * t312 + t405 * t204 + ((-t204 * t77 - t207 * t76) * qJD(4) + t389) * m(5) + (t312 * t67 + t345) * m(6) + t150 * t207) * t192 + t67 * (mrSges(6,1) * t220 + mrSges(6,2) * t221) + t31 * (mrSges(7,1) * t220 - mrSges(7,3) * t221) + (-t116 * t134 + t193 * t62) * m(5) + t284 * t380 + (-t19 * t221 - t20 * t220 - t323 * t4) * mrSges(6,3) + t225 * t269 - (t203 * t409 + t206 * t94) * t309 / 0.2e1 + (t119 * t134 - t120 * t137 + (t197 * t70 + t200 * t69) * pkin(2)) * m(4) - t182 * t286 + t193 * t118 + (-t312 * t76 + t389) * mrSges(5,3) + t134 * t168 + (Ifges(5,4) * t419 + Ifges(5,2) * t420 - Ifges(7,6) * t377 - Ifges(6,6) * t378 + t352 * t269 - t433 * t379 + t431 * t376 + t382 + Ifges(5,6) * qJDD(4) - t390 / 0.2e1 + (-m(5) * t77 - t182) * t137) * t207 + (-t237 * t374 - t244 * t375 - t370 * t396 - t372 * t394) * t309 + (t214 * t374 + t215 * t375 + t241 * t435 + t400 * t370 + t399 * t372 + t226) * qJD(4) + t411 * t126 + t412 * t125 + t413 * t127 + t414 * t128 + (t1 * t104 + t105 * t2 + t124 * t5 + t13 * t413 + t14 * t414 + t31 * t61) * m(7) + t415 * t323 / 0.2e1 - t62 * t255 + (-m(4) * t188 - m(5) * t261 - t143 * mrSges(5,3) + t142 * t255 + t421 * (-pkin(9) * t337 - t142 * t366 + t261) - t268 * (-t142 * t322 + t143 * t203) + t267 * (-t142 * t324 - t143 * t206) + t417 * t337 - t398) * g(3) + (-m(4) * t222 - t230 * mrSges(3,1) - (t198 * t327 - t201 * t208) * mrSges(3,2) - m(5) * t212 + t406 * t89 + t418 * t88 - t417 * t340 + t421 * (pkin(9) * t340 + t89 * t366 + t212) - t268 * (-t203 * t88 + t322 * t89) + t267 * (t88 * t206 + t324 * t89)) * g(1) + (Ifges(5,5) * qJDD(4) + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) - t388) * t203 + Ifges(5,1) * t175 + Ifges(5,4) * t420 + t238 * t377 + t245 * t378 + t5 * t251 - t269 * t408 + t91 * t273 + t376 * t395 + t379 * t393 + (m(5) * t76 - m(7) * t31 - t112 + t384) * t137) * t204 + t104 * t58 + t105 * t56 + t109 * t55 + t61 * t112 + t124 * t32 + (qJD(2) * t137 - t197 * t339 - t70) * mrSges(4,2); m(4) * t186 + ((t203 * t319 + t206 * t320 + t182) * qJD(4) + m(5) * (qJD(4) * t77 + t12) + m(7) * (t13 * t314 + t14 * t307 - t5) + m(6) * (-t19 * t314 + t20 * t307 - t10) - t304) * t207 + (t150 + t428 + t427 + t291 * qJD(4) + (-t203 * t320 + t206 * t319) * qJD(5) + m(5) * (-qJD(4) * t76 + t11) + m(7) * (qJD(4) * t31 + t425) + m(6) * (qJD(4) * t67 + t426)) * t204 + (-t202 * g(3) + (-g(1) * t198 + g(2) * t201) * t199) * t299; -(-Ifges(5,2) * t315 + t194 + t423) * t306 / 0.2e1 + (t383 + t425) * mrSges(7,2) + t422 * t315 + (t298 + t384) * t77 + (t358 * t393 + t359 * t399) * t167 + (-t245 / 0.2e1 + t238 / 0.2e1) * t311 + (-t182 + t297) * t76 + pkin(9) * t427 + pkin(9) * t428 + (-t19 * t38 - t20 * t39 - pkin(4) * t10 + ((-t19 * t206 - t20 * t203) * qJD(5) + t258) * pkin(9)) * m(6) + (t358 * t395 + t359 * t400) * t189 + t237 * t377 + t244 * t378 + (t380 + t91 / 0.2e1) * t310 + (-t19 * (mrSges(6,1) * t204 - mrSges(6,3) * t322) - t13 * (-mrSges(7,1) * t204 + mrSges(7,2) * t322) - t20 * (-mrSges(6,2) * t204 - mrSges(6,3) * t324) - t14 * (-mrSges(7,2) * t324 + mrSges(7,3) * t204) - t226 + (t215 / 0.2e1 - t214 / 0.2e1) * t166) * qJD(2) - t381 * t225 / 0.2e1 + Ifges(5,3) * qJDD(4) - t241 * t305 / 0.2e1 + (-t303 - t35) * t128 + Ifges(5,6) * t174 + Ifges(5,5) * t175 + t179 * t32 + t388 * t206 + t394 * t379 + t396 * t376 + (t369 * t94 - t397) * t306 + t397 * qJD(5) + (-t303 - t39) * t125 + t407 * t112 + (-t13 * t36 - t14 * t35 + t179 * t5 + ((t13 * t206 - t14 * t203) * qJD(5) + t259) * pkin(9) + t407 * t31) * m(7) - t11 * mrSges(5,2) + t12 * mrSges(5,1) + t409 * t273 + t415 * t369 + t5 * t252 - t10 * t254 + (mrSges(5,2) * t115 + t421 * (-t114 * pkin(4) + pkin(9) * t115) - t387 * t114) * g(3) - pkin(4) * t33 + (mrSges(5,2) * t54 + t421 * (t53 * pkin(4) + pkin(9) * t54) + t387 * t53) * g(1) + (mrSges(5,2) * t52 + t421 * (t51 * pkin(4) + pkin(9) * t52) + t387 * t51) * g(2) + (-t302 - t38) * t126 + (t383 + t426) * mrSges(6,3) + (t302 - t36) * t127; (-Ifges(6,2) * t167 - t162 + t409) * t374 + t94 * t372 - (-t166 * t433 - t167 * t432) * t189 / 0.2e1 + (-t166 * t434 + t161 - t351 + t91) * t436 - t382 + (Ifges(7,3) * t167 - t348) * t375 + (-t319 + t354) * t20 + (-t320 - t355) * t19 + (t13 * t166 + t14 * t167) * mrSges(7,2) + (-pkin(5) * t2 + qJ(6) * t1 - t111 * t31 - t13 * t20 + t14 * t404) * m(7) + (t267 * (-t203 * t86 + t206 * t52) + t268 * t15) * g(2) + (t267 * (-t203 * t89 + t206 * t54) + t268 * t17) * g(1) + (t267 * t46 + t268 * t45) * g(3) - t31 * (mrSges(7,1) * t167 + mrSges(7,3) * t166) - t67 * (mrSges(6,1) * t167 - mrSges(6,2) * t166) - pkin(5) * t56 + qJ(6) * t58 - t111 * t112 + qJD(6) * t128 + t390; t167 * t112 - t189 * t128 + (-g(1) * t17 - g(2) * t15 - g(3) * t45 - t14 * t189 + t31 * t167 + t2) * m(7) + t56;];
tau  = t8;
