% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:24
% EndTime: 2019-03-09 02:20:45
% DurationCPUTime: 14.58s
% Computational Cost: add. (10586->647), mult. (23503->841), div. (0->0), fcn. (17228->18), ass. (0->293)
t415 = mrSges(6,3) + mrSges(7,3);
t222 = sin(qJ(5));
t229 = -pkin(9) - pkin(8);
t269 = qJD(5) * t229;
t227 = cos(qJ(4));
t218 = cos(pkin(11));
t285 = qJD(1) * t218;
t216 = sin(pkin(11));
t223 = sin(qJ(4));
t292 = t216 * t223;
t170 = -qJD(1) * t292 + t227 * t285;
t305 = t170 * t222;
t179 = t216 * t227 + t218 * t223;
t171 = t179 * qJD(1);
t122 = pkin(4) * t171 - pkin(8) * t170;
t226 = cos(qJ(5));
t217 = sin(pkin(10));
t194 = pkin(1) * t217 + qJ(3);
t187 = t194 * qJD(1);
t203 = t218 * qJD(2);
t149 = t203 + (-pkin(7) * qJD(1) - t187) * t216;
t162 = t216 * qJD(2) + t218 * t187;
t150 = pkin(7) * t285 + t162;
t83 = t149 * t227 - t223 * t150;
t58 = t222 * t122 + t226 * t83;
t414 = pkin(9) * t305 + t222 * t269 - t58;
t304 = t170 * t226;
t57 = t226 * t122 - t222 * t83;
t413 = -pkin(5) * t171 + pkin(9) * t304 + t226 * t269 - t57;
t199 = pkin(5) * t226 + pkin(4);
t215 = qJ(5) + qJ(6);
t208 = sin(t215);
t209 = cos(t215);
t254 = -mrSges(6,1) * t226 + mrSges(6,2) * t222;
t412 = -m(6) * pkin(4) - m(7) * t199 - mrSges(7,1) * t209 + mrSges(7,2) * t208 + t254;
t411 = -m(6) * pkin(8) + m(7) * t229 - t415;
t176 = qJD(1) * qJD(3) + qJDD(1) * t194;
t282 = qJD(5) * t222;
t410 = t282 - t305;
t197 = pkin(3) * t218 + pkin(2);
t219 = cos(pkin(10));
t336 = pkin(1) * t219;
t186 = -t197 - t336;
t169 = qJD(1) * t186 + qJD(3);
t387 = Ifges(5,5) * qJD(4);
t409 = -t169 * mrSges(5,2) - t387 / 0.2e1;
t173 = t179 * qJD(4);
t225 = cos(qJ(6));
t221 = sin(qJ(6));
t151 = qJD(4) * t226 - t171 * t222;
t84 = t149 * t223 + t150 * t227;
t78 = qJD(4) * pkin(8) + t84;
t95 = -pkin(4) * t170 - pkin(8) * t171 + t169;
t51 = t222 * t95 + t226 * t78;
t36 = pkin(9) * t151 + t51;
t313 = t221 * t36;
t167 = qJD(5) - t170;
t152 = qJD(4) * t222 + t171 * t226;
t50 = -t222 * t78 + t226 * t95;
t35 = -pkin(9) * t152 + t50;
t34 = pkin(5) * t167 + t35;
t13 = t225 * t34 - t313;
t311 = t225 * t36;
t14 = t221 * t34 + t311;
t386 = Ifges(5,6) * qJD(4);
t408 = -t169 * mrSges(5,1) - t50 * mrSges(6,1) - t13 * mrSges(7,1) + t51 * mrSges(6,2) + t14 * mrSges(7,2) + t386 / 0.2e1;
t165 = Ifges(5,4) * t170;
t406 = t151 * Ifges(6,6);
t405 = t167 * Ifges(6,3);
t403 = t170 * Ifges(5,2);
t402 = t216 ^ 2 + t218 ^ 2;
t366 = m(7) * pkin(5);
t261 = t225 * t151 - t152 * t221;
t178 = -t227 * t218 + t292;
t172 = t178 * qJD(4);
t124 = -qJD(1) * t172 + qJDD(1) * t179;
t74 = qJD(5) * t151 + qJDD(4) * t222 + t124 * t226;
t75 = -qJD(5) * t152 + qJDD(4) * t226 - t124 * t222;
t27 = qJD(6) * t261 + t221 * t75 + t225 * t74;
t365 = t27 / 0.2e1;
t88 = t151 * t221 + t152 * t225;
t28 = -qJD(6) * t88 - t221 * t74 + t225 * t75;
t364 = t28 / 0.2e1;
t357 = t74 / 0.2e1;
t356 = t75 / 0.2e1;
t400 = -m(4) - m(3);
t277 = qJDD(1) * t218;
t278 = qJDD(1) * t216;
t125 = -qJD(1) * t173 - t223 * t278 + t227 * t277;
t121 = qJDD(5) - t125;
t119 = qJDD(6) + t121;
t351 = t119 / 0.2e1;
t350 = t121 / 0.2e1;
t188 = t229 * t222;
t189 = t229 * t226;
t148 = t188 * t221 - t189 * t225;
t395 = -qJD(6) * t148 - t221 * t414 + t413 * t225;
t147 = t188 * t225 + t189 * t221;
t394 = qJD(6) * t147 + t413 * t221 + t225 * t414;
t163 = qJD(6) + t167;
t393 = t152 * Ifges(6,5) + t88 * Ifges(7,5) + Ifges(7,6) * t261 + t163 * Ifges(7,3) + t405 + t406;
t392 = mrSges(6,1) + t366;
t391 = pkin(5) * t410 - t84;
t37 = -mrSges(6,1) * t75 + mrSges(6,2) * t74;
t390 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t124 + t37;
t325 = mrSges(5,3) * t171;
t389 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t151 + mrSges(6,2) * t152 + t325;
t181 = t221 * t226 + t222 * t225;
t377 = qJD(5) + qJD(6);
t130 = t377 * t181;
t99 = t181 * t170;
t388 = t99 - t130;
t243 = t221 * t222 - t225 * t226;
t115 = t243 * t179;
t100 = t243 * t170;
t129 = t377 * t243;
t385 = t100 - t129;
t327 = pkin(7) + t194;
t174 = t327 * t216;
t175 = t327 * t218;
t117 = -t174 * t223 + t175 * t227;
t103 = t226 * t117;
t104 = pkin(4) * t178 - pkin(8) * t179 + t186;
t65 = t222 * t104 + t103;
t384 = -t227 * t174 - t175 * t223;
t281 = qJD(5) * t226;
t239 = -t222 * t172 + t179 * t281;
t201 = t218 * qJDD(2);
t153 = -t176 * t216 + t201;
t154 = t216 * qJDD(2) + t218 * t176;
t382 = -t153 * t216 + t154 * t218;
t54 = mrSges(6,1) * t121 - mrSges(6,3) * t74;
t55 = -mrSges(6,2) * t121 + mrSges(6,3) * t75;
t381 = -t222 * t54 + t226 * t55;
t134 = t201 + (-pkin(7) * qJDD(1) - t176) * t216;
t135 = pkin(7) * t277 + t154;
t283 = qJD(4) * t227;
t284 = qJD(4) * t223;
t47 = t223 * t134 + t227 * t135 + t149 * t283 - t150 * t284;
t43 = qJDD(4) * pkin(8) + t47;
t166 = qJDD(1) * t186 + qJDD(3);
t62 = -pkin(4) * t125 - pkin(8) * t124 + t166;
t11 = t222 * t62 + t226 * t43 + t95 * t281 - t282 * t78;
t12 = -qJD(5) * t51 - t222 * t43 + t226 * t62;
t380 = t11 * t226 - t12 * t222;
t379 = -m(5) - m(7) - m(6);
t77 = -qJD(4) * pkin(4) - t83;
t63 = -pkin(5) * t151 + t77;
t375 = -mrSges(7,1) * t63 + t14 * mrSges(7,3);
t374 = mrSges(7,2) * t63 - t13 * mrSges(7,3);
t373 = -m(6) * t77 - t389;
t6 = pkin(5) * t121 - pkin(9) * t74 + t12;
t9 = pkin(9) * t75 + t11;
t2 = qJD(6) * t13 + t221 * t6 + t225 * t9;
t3 = -qJD(6) * t14 - t221 * t9 + t225 * t6;
t372 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t371 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t213 = pkin(11) + qJ(4);
t204 = sin(t213);
t206 = cos(t213);
t256 = mrSges(5,1) * t206 - mrSges(5,2) * t204;
t257 = -mrSges(4,1) * t218 + mrSges(4,2) * t216;
t370 = m(4) * pkin(2) + t204 * t415 + mrSges(3,1) + t256 - t257;
t369 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t368 = Ifges(7,4) * t365 + Ifges(7,2) * t364 + Ifges(7,6) * t351;
t367 = Ifges(7,1) * t365 + Ifges(7,4) * t364 + Ifges(7,5) * t351;
t363 = Ifges(6,1) * t357 + Ifges(6,4) * t356 + Ifges(6,5) * t350;
t337 = Ifges(7,4) * t88;
t40 = Ifges(7,2) * t261 + Ifges(7,6) * t163 + t337;
t362 = -t40 / 0.2e1;
t361 = t40 / 0.2e1;
t82 = Ifges(7,4) * t261;
t41 = Ifges(7,1) * t88 + Ifges(7,5) * t163 + t82;
t360 = -t41 / 0.2e1;
t359 = t41 / 0.2e1;
t355 = -t261 / 0.2e1;
t354 = t261 / 0.2e1;
t353 = -t88 / 0.2e1;
t352 = t88 / 0.2e1;
t349 = -t151 / 0.2e1;
t348 = -t152 / 0.2e1;
t347 = t152 / 0.2e1;
t346 = -t163 / 0.2e1;
t345 = t163 / 0.2e1;
t344 = -t167 / 0.2e1;
t343 = -t170 / 0.2e1;
t341 = t171 / 0.2e1;
t338 = t226 / 0.2e1;
t224 = sin(qJ(1));
t335 = pkin(1) * t224;
t334 = pkin(5) * t152;
t331 = g(3) * t204;
t228 = cos(qJ(1));
t210 = t228 * pkin(1);
t326 = mrSges(5,3) * t170;
t324 = mrSges(6,3) * t151;
t323 = mrSges(6,3) * t152;
t322 = Ifges(5,4) * t171;
t321 = Ifges(6,4) * t152;
t320 = Ifges(6,4) * t222;
t319 = Ifges(6,4) * t226;
t302 = t179 * t222;
t301 = t179 * t226;
t214 = qJ(1) + pkin(10);
t205 = sin(t214);
t300 = t205 * t208;
t299 = t205 * t209;
t298 = t205 * t222;
t297 = t205 * t226;
t207 = cos(t214);
t296 = t207 * t208;
t295 = t207 * t209;
t294 = t207 * t222;
t293 = t207 * t226;
t288 = t226 * t172;
t143 = t206 * t300 + t295;
t144 = -t206 * t299 + t296;
t287 = -t143 * mrSges(7,1) + t144 * mrSges(7,2);
t145 = -t206 * t296 + t299;
t146 = t206 * t295 + t300;
t286 = t145 * mrSges(7,1) - t146 * mrSges(7,2);
t276 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t119;
t53 = -mrSges(7,1) * t261 + mrSges(7,2) * t88;
t275 = t53 + t389;
t274 = Ifges(6,5) * t74 + Ifges(6,6) * t75 + Ifges(6,3) * t121;
t272 = m(4) - t379;
t142 = Ifges(6,4) * t151;
t71 = t152 * Ifges(6,1) + t167 * Ifges(6,5) + t142;
t270 = t71 * t338;
t198 = -pkin(2) - t336;
t264 = -t282 / 0.2e1;
t263 = -t125 * mrSges(5,1) + t124 * mrSges(5,2);
t123 = pkin(4) * t173 + pkin(8) * t172;
t90 = -t178 * qJD(3) + qJD(4) * t384;
t262 = t226 * t123 - t222 * t90;
t64 = t226 * t104 - t117 * t222;
t260 = pkin(4) * t206 + pkin(8) * t204;
t258 = -mrSges(4,1) * t277 + mrSges(4,2) * t278;
t253 = mrSges(6,1) * t222 + mrSges(6,2) * t226;
t252 = -mrSges(7,1) * t208 - mrSges(7,2) * t209;
t251 = Ifges(6,1) * t226 - t320;
t250 = -Ifges(6,2) * t222 + t319;
t249 = Ifges(6,5) * t226 - Ifges(6,6) * t222;
t52 = pkin(5) * t178 - pkin(9) * t301 + t64;
t56 = -pkin(9) * t302 + t65;
t22 = -t221 * t56 + t225 * t52;
t23 = t221 * t52 + t225 * t56;
t247 = -t50 * t222 + t51 * t226;
t96 = -mrSges(6,2) * t167 + t324;
t97 = mrSges(6,1) * t167 - t323;
t246 = -t222 * t97 + t226 * t96;
t245 = -(-t187 * t216 + t203) * t216 + t162 * t218;
t244 = t199 * t206 - t204 * t229;
t48 = t134 * t227 - t223 * t135 - t149 * t284 - t150 * t283;
t242 = t276 + t372;
t159 = -qJD(4) * mrSges(5,2) + t326;
t241 = -t159 - t246;
t157 = -t206 * t294 + t297;
t155 = t206 * t298 + t293;
t240 = t77 * t253;
t238 = t179 * t282 + t288;
t30 = t104 * t281 - t117 * t282 + t222 * t123 + t226 * t90;
t44 = -qJDD(4) * pkin(4) - t48;
t234 = (-t222 * t51 - t226 * t50) * qJD(5) + t380;
t91 = qJD(3) * t179 + qJD(4) * t117;
t220 = -pkin(7) - qJ(3);
t185 = qJDD(1) * t198 + qJDD(3);
t158 = t206 * t293 + t298;
t156 = -t206 * t297 + t294;
t114 = t181 * t179;
t108 = t171 * Ifges(5,1) + t165 + t387;
t107 = t322 + t386 + t403;
t106 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t125;
t92 = pkin(5) * t302 - t384;
t70 = t151 * Ifges(6,2) + t167 * Ifges(6,6) + t321;
t67 = mrSges(7,1) * t163 - mrSges(7,3) * t88;
t66 = -mrSges(7,2) * t163 + mrSges(7,3) * t261;
t59 = pkin(5) * t239 + t91;
t46 = t115 * t377 + t181 * t172;
t45 = -t130 * t179 + t172 * t243;
t32 = t74 * Ifges(6,4) + t75 * Ifges(6,2) + t121 * Ifges(6,6);
t31 = -qJD(5) * t65 + t262;
t29 = -pkin(5) * t75 + t44;
t24 = -pkin(9) * t239 + t30;
t21 = pkin(9) * t288 + pkin(5) * t173 + (-t103 + (pkin(9) * t179 - t104) * t222) * qJD(5) + t262;
t20 = -mrSges(7,2) * t119 + mrSges(7,3) * t28;
t19 = mrSges(7,1) * t119 - mrSges(7,3) * t27;
t16 = t225 * t35 - t313;
t15 = -t221 * t35 - t311;
t10 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t5 = -qJD(6) * t23 + t21 * t225 - t221 * t24;
t4 = qJD(6) * t22 + t21 * t221 + t225 * t24;
t1 = [(-t270 - Ifges(5,1) * t341 + t83 * mrSges(5,3) - t165 / 0.2e1 - t108 / 0.2e1 + t409) * t172 + m(4) * (t245 * qJD(3) + t185 * t198 + t194 * t382) - t239 * t70 / 0.2e1 - (-m(5) * t48 + m(6) * t44 + t390) * t384 + (-t298 * t366 - mrSges(2,1) * t228 - t158 * mrSges(6,1) - t146 * mrSges(7,1) + mrSges(2,2) * t224 - t157 * mrSges(6,2) - t145 * mrSges(7,2) + t379 * (t207 * t197 - t205 * t220 + t210) + t400 * t210 + t369 * t205 + (-m(6) * t260 - m(7) * t244 - t370) * t207) * g(2) + (-t294 * t366 + mrSges(2,1) * t224 - t156 * mrSges(6,1) - t144 * mrSges(7,1) + mrSges(2,2) * t228 - t155 * mrSges(6,2) - t143 * mrSges(7,2) - t400 * t335 + t379 * (-t207 * t220 - t335) + t369 * t207 + (m(5) * t197 - m(7) * (-t197 - t244) - m(6) * (-t197 - t260) + t370) * t205) * g(1) + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t354 + (t176 * t402 + t382) * mrSges(4,3) + t151 * (-Ifges(6,4) * t238 - Ifges(6,2) * t239) / 0.2e1 + (Ifges(7,5) * t45 + Ifges(7,6) * t46) * t345 + t167 * (-Ifges(6,5) * t238 - Ifges(6,6) * t239) / 0.2e1 + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t352 + m(6) * (t11 * t65 + t12 * t64 + t30 * t51 + t31 * t50) + m(5) * (t117 * t47 + t166 * t186 + t84 * t90) + (-Ifges(6,1) * t238 - Ifges(6,4) * t239) * t347 + (-t11 * t302 - t12 * t301 + t238 * t50 - t239 * t51) * mrSges(6,3) - t115 * t367 - t114 * t368 + t45 * t359 + t46 * t361 + t301 * t363 + m(7) * (t13 * t5 + t14 * t4 + t2 * t23 + t22 * t3 + t29 * t92 + t59 * t63) + t77 * (mrSges(6,1) * t239 - mrSges(6,2) * t238) + t185 * t257 + (mrSges(5,1) * t166 - mrSges(5,3) * t47 - Ifges(5,4) * t124 + Ifges(6,5) * t357 + Ifges(7,5) * t365 - Ifges(5,2) * t125 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t356 + Ifges(7,6) * t364 + Ifges(6,3) * t350 + Ifges(7,3) * t351 + t371 + t372) * t178 + (-m(5) * t83 - t373) * t91 + (Ifges(7,5) * t352 + Ifges(7,6) * t354 + Ifges(7,3) * t345 + Ifges(6,5) * t347 + t393 / 0.2e1 + t405 / 0.2e1 + t406 / 0.2e1 - Ifges(5,4) * t341 - t84 * mrSges(5,3) - t403 / 0.2e1 - t107 / 0.2e1 - t408) * t173 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t219 - 0.2e1 * mrSges(3,2) * t217 + m(3) * (t217 ^ 2 + t219 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-Ifges(7,4) * t115 - Ifges(7,2) * t114) * t364 + (-Ifges(7,5) * t115 - Ifges(7,6) * t114) * t351 + (-Ifges(7,1) * t115 - Ifges(7,4) * t114) * t365 + (-t114 * t2 + t115 * t3 - t13 * t45 + t14 * t46) * mrSges(7,3) + t29 * (mrSges(7,1) * t114 - mrSges(7,2) * t115) + (t276 + t274) * t178 / 0.2e1 + t198 * t258 + (t166 * mrSges(5,2) - t48 * mrSges(5,3) + Ifges(5,1) * t124 + Ifges(5,4) * t125 + Ifges(5,5) * qJDD(4) + t249 * t350 + t250 * t356 + t251 * t357 + t253 * t44 + t264 * t71) * t179 + (Ifges(4,4) * t216 + Ifges(4,2) * t218) * t277 + (Ifges(4,1) * t216 + Ifges(4,4) * t218) * t278 + t186 * t263 + t22 * t19 + t23 * t20 + t59 * t53 + t63 * (-mrSges(7,1) * t46 + mrSges(7,2) * t45) + t64 * t54 + t65 * t55 + t4 * t66 + t5 * t67 + t92 * t10 + t30 * t96 + t31 * t97 + t117 * t106 + t90 * t159 - t32 * t302 / 0.2e1; m(3) * qJDD(2) - t114 * t19 - t115 * t20 + t45 * t66 + t46 * t67 + (t10 + t390) * t178 + t275 * t173 + t241 * t172 + (t106 + (-t222 * t96 - t226 * t97) * qJD(5) + t381) * t179 + (-m(3) - t272) * g(3) + m(7) * (-t114 * t3 - t115 * t2 + t13 * t46 + t14 * t45 + t173 * t63 + t178 * t29) + m(5) * (-t172 * t84 - t173 * t83 - t178 * t48 + t179 * t47) + m(4) * (t153 * t218 + t154 * t216) + m(6) * (-t172 * t247 + t173 * t77 + t178 * t44 + t179 * t234); t263 + t258 - t402 * qJD(1) ^ 2 * mrSges(4,3) + t241 * t170 + t246 * qJD(5) + t388 * t67 + t385 * t66 + t226 * t54 + t222 * t55 - t275 * t171 - t243 * t19 + t181 * t20 + (-g(1) * t205 + g(2) * t207) * t272 + (t13 * t388 + t14 * t385 - t171 * t63 + t181 * t2 - t243 * t3) * m(7) + (t11 * t222 + t12 * t226 + t167 * t247 - t171 * t77) * m(6) + (-t170 * t84 + t171 * t83 + t166) * m(5) + (-qJD(1) * t245 + t185) * m(4); (t249 * t344 + t250 * t349 + t251 * t348 - t240 + t409) * t170 + (-t410 * t51 + (-t281 + t304) * t50 + t380) * mrSges(6,3) + (m(6) * t234 - t281 * t97 - t282 * t96 + t381) * pkin(8) - t243 * t368 + (Ifges(7,4) * t181 - Ifges(7,2) * t243) * t364 + (Ifges(7,1) * t181 - Ifges(7,4) * t243) * t365 + (Ifges(7,5) * t181 - Ifges(7,6) * t243) * t351 + (-t100 * t13 + t14 * t99 - t181 * t3 - t2 * t243) * mrSges(7,3) + t29 * (mrSges(7,1) * t243 + mrSges(7,2) * t181) + (t207 * g(1) + t205 * g(2)) * ((mrSges(5,2) + t411) * t206 + (mrSges(5,1) - t412) * t204) + (t204 * t411 + t206 * t412 - t256) * g(3) + (-pkin(4) * t44 - t50 * t57 - t51 * t58) * m(6) + t222 * t363 + t181 * t367 + (Ifges(6,5) * t222 + Ifges(6,6) * t226) * t350 + (Ifges(6,2) * t226 + t320) * t356 + (Ifges(6,1) * t222 + t319) * t357 - t100 * t360 - t99 * t362 + t391 * t53 + (t264 + t305 / 0.2e1) * t70 + t32 * t338 + t107 * t341 - (Ifges(5,1) * t170 - t322 + t393) * t171 / 0.2e1 + t394 * t66 + t395 * t67 + (t13 * t395 + t14 * t394 + t147 * t3 + t148 * t2 - t199 * t29 + t391 * t63) * m(7) + (t270 + t240) * qJD(5) + t44 * t254 + (t325 + t373) * t84 - (Ifges(7,1) * t352 + Ifges(7,4) * t354 + Ifges(7,5) * t345 + t359 + t374) * t129 - (Ifges(7,4) * t352 + Ifges(7,2) * t354 + Ifges(7,6) * t345 + t361 + t375) * t130 + (Ifges(6,5) * t348 + Ifges(7,5) * t353 - Ifges(5,2) * t343 + Ifges(6,6) * t349 + Ifges(7,6) * t355 + Ifges(6,3) * t344 + Ifges(7,3) * t346 + t408) * t171 + (t165 + t108) * t343 + (t326 - t159) * t83 + (-Ifges(7,5) * t100 - Ifges(7,6) * t99) * t346 + (-Ifges(7,4) * t100 - Ifges(7,2) * t99) * t355 + (-Ifges(7,1) * t100 - Ifges(7,4) * t99) * t353 - t63 * (mrSges(7,1) * t99 - mrSges(7,2) * t100) + (t151 * t250 + t152 * t251 + t167 * t249) * qJD(5) / 0.2e1 + Ifges(5,3) * qJDD(4) - t199 * t10 - pkin(4) * t37 - t47 * mrSges(5,2) + t48 * mrSges(5,1) - t58 * t96 - t57 * t97 + Ifges(5,5) * t124 + Ifges(5,6) * t125 + t147 * t19 + t148 * t20 - t71 * t304 / 0.2e1; t371 + (t323 + t97) * t51 + (t2 * t221 + t225 * t3 + (-t13 * t221 + t14 * t225) * qJD(6)) * t366 + (t222 * t366 - t252 + t253) * t331 + (t324 - t96) * t50 + (Ifges(6,5) * t151 - Ifges(6,6) * t152) * t344 + t70 * t347 + (Ifges(6,1) * t151 - t321) * t348 + (-Ifges(6,2) * t152 + t142 + t71) * t349 + (mrSges(6,2) * t158 - t157 * t392 - t286) * g(1) + t274 + (-mrSges(6,2) * t156 + t155 * t392 - t287) * g(2) - t53 * t334 - m(7) * (t13 * t15 + t14 * t16 + t334 * t63) - (Ifges(7,4) * t353 + Ifges(7,2) * t355 + Ifges(7,6) * t346 + t362 - t375) * t88 + (Ifges(7,1) * t353 + Ifges(7,4) * t355 + Ifges(7,5) * t346 + t360 - t374) * t261 + t242 - t16 * t66 - t15 * t67 - t77 * (mrSges(6,1) * t152 + mrSges(6,2) * t151) + ((-t221 * t67 + t225 * t66) * qJD(6) + t19 * t225 + t20 * t221) * pkin(5); -t63 * (mrSges(7,1) * t88 + mrSges(7,2) * t261) + (Ifges(7,1) * t261 - t337) * t353 + t40 * t352 + (Ifges(7,5) * t261 - Ifges(7,6) * t88) * t346 - t13 * t66 + t14 * t67 - g(1) * t286 - g(2) * t287 - t252 * t331 + (t13 * t261 + t14 * t88) * mrSges(7,3) + t242 + (-Ifges(7,2) * t88 + t41 + t82) * t355;];
tau  = t1;
