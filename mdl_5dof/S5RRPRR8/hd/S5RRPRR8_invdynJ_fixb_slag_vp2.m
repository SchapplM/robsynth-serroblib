% Calculate vector of inverse dynamics joint torques for
% S5RRPRR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:54
% EndTime: 2019-12-31 20:17:17
% DurationCPUTime: 13.20s
% Computational Cost: add. (8471->570), mult. (19860->773), div. (0->0), fcn. (14710->14), ass. (0->278)
t219 = qJ(2) + pkin(9);
t214 = qJ(4) + t219;
t203 = sin(t214);
t204 = cos(t214);
t227 = cos(qJ(5));
t323 = t227 * mrSges(6,1);
t376 = m(6) * pkin(4);
t232 = mrSges(5,2) * t204 + (mrSges(5,1) + t323 + t376) * t203;
t405 = -m(6) - m(5);
t223 = sin(qJ(5));
t342 = mrSges(6,2) * t223;
t416 = -t203 * t342 + t204 * (-m(6) * pkin(8) - mrSges(6,3));
t388 = t204 * pkin(4) + t203 * pkin(8);
t415 = m(6) * t388;
t221 = cos(pkin(9));
t359 = pkin(2) * t221;
t205 = pkin(3) + t359;
t224 = sin(qJ(4));
t228 = cos(qJ(4));
t220 = sin(pkin(9));
t360 = pkin(2) * t220;
t152 = t224 * t205 + t228 * t360;
t222 = -qJ(3) - pkin(6);
t225 = sin(qJ(2));
t191 = t222 * t225;
t178 = qJD(1) * t191;
t229 = cos(qJ(2));
t193 = t222 * t229;
t179 = qJD(1) * t193;
t309 = t221 * t179;
t121 = -t178 * t220 + t309;
t170 = -t220 * t225 + t221 * t229;
t156 = t170 * qJD(1);
t356 = pkin(7) * t156;
t247 = t121 - t356;
t160 = t220 * t179;
t122 = t221 * t178 + t160;
t301 = qJD(1) * t229;
t302 = qJD(1) * t225;
t157 = -t220 * t301 - t221 * t302;
t355 = pkin(7) * t157;
t92 = t122 + t355;
t390 = -t152 * qJD(4) + t224 * t92 - t228 * t247;
t192 = -t229 * mrSges(3,1) + t225 * mrSges(3,2);
t212 = sin(t219);
t213 = cos(t219);
t414 = -t213 * mrSges(4,1) + t212 * mrSges(4,2) + t192;
t413 = -t204 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t203;
t265 = mrSges(3,1) * t225 + mrSges(3,2) * t229;
t358 = pkin(2) * t225;
t412 = m(4) * t358 + mrSges(4,1) * t212 + mrSges(4,2) * t213 + t265 + t405 * (-pkin(3) * t212 - t358) + t232;
t167 = qJD(2) * pkin(2) + t178;
t115 = t221 * t167 + t160;
t84 = qJD(2) * pkin(3) + t115 + t355;
t116 = t220 * t167 - t309;
t87 = t116 + t356;
t49 = -t224 * t87 + t228 * t84;
t411 = t323 - t342;
t215 = t229 * pkin(2);
t206 = t215 + pkin(1);
t183 = -qJD(1) * t206 + qJD(3);
t125 = -pkin(3) * t156 + t183;
t218 = qJD(2) + qJD(4);
t328 = t218 * Ifges(5,6);
t251 = t156 * t224 - t228 * t157;
t337 = Ifges(5,4) * t251;
t50 = t224 * t84 + t228 * t87;
t349 = t50 * mrSges(5,3);
t270 = t228 * t156 + t157 * t224;
t101 = qJD(5) - t270;
t334 = t101 * Ifges(6,3);
t91 = t218 * t223 + t227 * t251;
t347 = t91 * Ifges(6,5);
t90 = t218 * t227 - t223 * t251;
t348 = t90 * Ifges(6,6);
t39 = t334 + t347 + t348;
t45 = pkin(8) * t218 + t50;
t61 = -pkin(4) * t270 - pkin(8) * t251 + t125;
t21 = t223 * t61 + t227 * t45;
t397 = t21 * mrSges(6,2);
t20 = -t223 * t45 + t227 * t61;
t398 = t20 * mrSges(6,1);
t392 = t270 * Ifges(5,2);
t71 = t328 + t337 + t392;
t410 = -t125 * mrSges(5,1) + t349 - t39 / 0.2e1 + t397 - t398 + t71 / 0.2e1 + t337 / 0.2e1 + t328 / 0.2e1;
t100 = Ifges(5,4) * t270;
t262 = mrSges(6,1) * t223 + mrSges(6,2) * t227;
t44 = -pkin(4) * t218 - t49;
t243 = t44 * t262;
t88 = Ifges(6,4) * t90;
t41 = t91 * Ifges(6,1) + t101 * Ifges(6,5) + t88;
t321 = t227 * t41;
t329 = t218 * Ifges(5,5);
t350 = t49 * mrSges(5,3);
t362 = t223 / 0.2e1;
t361 = Ifges(6,4) * t91;
t40 = Ifges(6,2) * t90 + Ifges(6,6) * t101 + t361;
t393 = t251 * Ifges(5,1);
t72 = t100 + t329 + t393;
t409 = -t125 * mrSges(5,2) - t243 + t40 * t362 - t321 / 0.2e1 + t350 - t72 / 0.2e1 - t329 / 0.2e1 - t100 / 0.2e1;
t294 = qJD(1) * qJD(2);
t277 = t225 * t294;
t293 = qJDD(1) * t229;
t180 = -t277 + t293;
t181 = qJDD(1) * t225 + t229 * t294;
t124 = t180 * t220 + t181 * t221;
t169 = t181 * pkin(6);
t299 = qJD(3) * t225;
t113 = qJDD(2) * pkin(2) - qJ(3) * t181 - qJD(1) * t299 - t169;
t207 = pkin(6) * t293;
t300 = qJD(2) * t225;
t289 = pkin(6) * t300;
t298 = qJD(3) * t229;
t119 = qJ(3) * t180 + t207 + (-t289 + t298) * qJD(1);
t75 = t221 * t113 - t119 * t220;
t53 = qJDD(2) * pkin(3) - pkin(7) * t124 + t75;
t123 = t180 * t221 - t181 * t220;
t76 = t220 * t113 + t221 * t119;
t60 = pkin(7) * t123 + t76;
t16 = -qJD(4) * t50 - t224 * t60 + t228 * t53;
t216 = qJDD(2) + qJDD(4);
t58 = qJD(4) * t270 + t123 * t224 + t124 * t228;
t408 = m(5) * t16 + mrSges(5,1) * t216 - mrSges(5,3) * t58;
t403 = t180 / 0.2e1;
t402 = t181 / 0.2e1;
t346 = qJD(2) / 0.2e1;
t396 = -mrSges(5,1) * t218 - mrSges(6,1) * t90 + mrSges(6,2) * t91 + mrSges(5,3) * t251;
t171 = t220 * t229 + t221 * t225;
t395 = Ifges(4,5) * t171;
t394 = Ifges(4,6) * t170;
t318 = qJDD(2) / 0.2e1;
t66 = -mrSges(6,2) * t101 + mrSges(6,3) * t90;
t67 = mrSges(6,1) * t101 - mrSges(6,3) * t91;
t254 = -t223 * t67 + t227 * t66;
t93 = -mrSges(5,2) * t218 + mrSges(5,3) * t270;
t391 = -t254 - t93;
t387 = -t204 * t411 + t413;
t118 = t170 * t224 + t171 * t228;
t296 = qJD(5) * t227;
t158 = t171 * qJD(2);
t159 = t170 * qJD(2);
t250 = t228 * t170 - t171 * t224;
t77 = qJD(4) * t250 - t158 * t224 + t159 * t228;
t245 = t118 * t296 + t223 * t77;
t297 = qJD(5) * t223;
t386 = -t20 * t296 - t21 * t297;
t168 = -pkin(6) * t277 + t207;
t385 = t168 * t229 + t169 * t225;
t383 = -m(3) * pkin(6) + m(4) * t222 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t59 = -qJD(4) * t251 + t123 * t228 - t124 * t224;
t35 = qJD(5) * t90 + t216 * t223 + t227 * t58;
t57 = qJDD(5) - t59;
t17 = mrSges(6,1) * t57 - mrSges(6,3) * t35;
t36 = -qJD(5) * t91 + t216 * t227 - t223 * t58;
t18 = -mrSges(6,2) * t57 + mrSges(6,3) * t36;
t382 = -t223 * t17 + t227 * t18 - t67 * t296 - t66 * t297;
t15 = qJD(4) * t49 + t224 * t53 + t228 * t60;
t11 = pkin(8) * t216 + t15;
t317 = qJDD(1) * pkin(1);
t143 = -pkin(2) * t180 + qJDD(3) - t317;
t89 = -pkin(3) * t123 + t143;
t19 = -pkin(4) * t59 - pkin(8) * t58 + t89;
t2 = qJD(5) * t20 + t11 * t227 + t19 * t223;
t3 = -qJD(5) * t21 - t11 * t223 + t19 * t227;
t381 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t74 = pkin(4) * t251 - pkin(8) * t270;
t380 = m(3) * pkin(1) + m(4) * t206 + mrSges(2,1) - t413 - t414;
t375 = t35 / 0.2e1;
t374 = t36 / 0.2e1;
t371 = t57 / 0.2e1;
t368 = -t90 / 0.2e1;
t367 = -t91 / 0.2e1;
t366 = t91 / 0.2e1;
t365 = -t101 / 0.2e1;
t363 = -t157 / 0.2e1;
t357 = pkin(6) * t229;
t352 = t2 * t227;
t351 = t3 * t223;
t341 = mrSges(6,3) * t223;
t340 = mrSges(6,3) * t227;
t339 = Ifges(3,4) * t225;
t338 = Ifges(3,4) * t229;
t336 = Ifges(6,4) * t223;
t335 = Ifges(6,4) * t227;
t333 = t156 * mrSges(4,3);
t332 = t157 * mrSges(4,3);
t331 = t157 * Ifges(4,4);
t316 = t118 * t223;
t315 = t118 * t227;
t230 = cos(qJ(1));
t308 = t223 * t230;
t226 = sin(qJ(1));
t307 = t226 * t223;
t306 = t226 * t227;
t305 = t227 * t230;
t274 = qJD(2) * t222;
t153 = t225 * t274 + t298;
t154 = t229 * t274 - t299;
t97 = t221 * t153 + t220 * t154;
t127 = t220 * t191 - t221 * t193;
t303 = pkin(3) * t213 + t215;
t292 = Ifges(6,5) * t35 + Ifges(6,6) * t36 + Ifges(6,3) * t57;
t291 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t302) * t357;
t210 = pkin(2) * t300;
t209 = pkin(2) * t302;
t281 = t321 / 0.2e1;
t278 = -t59 * mrSges(5,1) + t58 * mrSges(5,2);
t275 = -t297 / 0.2e1;
t128 = -pkin(3) * t157 + t209;
t129 = pkin(3) * t158 + t210;
t272 = -t123 * mrSges(4,1) + t124 * mrSges(4,2);
t96 = -t153 * t220 + t221 * t154;
t126 = t221 * t191 + t193 * t220;
t132 = -pkin(3) * t170 - t206;
t268 = t416 * t226;
t267 = t416 * t230;
t261 = Ifges(6,1) * t227 - t336;
t260 = Ifges(3,2) * t229 + t339;
t259 = -Ifges(6,2) * t223 + t335;
t258 = Ifges(3,5) * t229 - Ifges(3,6) * t225;
t257 = Ifges(6,5) * t227 - Ifges(6,6) * t223;
t256 = t20 * t227 + t21 * t223;
t255 = -t20 * t223 + t21 * t227;
t68 = -pkin(4) * t250 - pkin(8) * t118 + t132;
t98 = -pkin(7) * t171 + t126;
t99 = pkin(7) * t170 + t127;
t70 = t224 * t98 + t228 * t99;
t28 = -t223 * t70 + t227 * t68;
t29 = t223 * t68 + t227 * t70;
t69 = t224 * t99 - t228 * t98;
t151 = t205 * t228 - t224 * t360;
t248 = -pkin(7) * t159 + t96;
t246 = pkin(1) * t265;
t244 = t118 * t297 - t227 * t77;
t242 = t90 * t259;
t241 = t91 * t261;
t240 = t101 * t257;
t239 = t225 * (Ifges(3,1) * t229 - t339);
t233 = -qJD(5) * t256 - t351;
t12 = -pkin(4) * t216 - t16;
t8 = t35 * Ifges(6,4) + t36 * Ifges(6,2) + t57 * Ifges(6,6);
t9 = Ifges(6,1) * t35 + Ifges(6,4) * t36 + Ifges(6,5) * t57;
t231 = -t15 * mrSges(5,2) + t2 * t340 + t227 * t8 / 0.2e1 + t9 * t362 - t12 * t411 + t16 * mrSges(5,1) + Ifges(5,3) * t216 + (Ifges(6,1) * t223 + t335) * t375 + (Ifges(6,2) * t227 + t336) * t374 + t40 * t275 + (Ifges(6,5) * t223 + Ifges(6,6) * t227) * t371 + Ifges(5,6) * t59 + Ifges(5,5) * t58 + (t281 + t243) * qJD(5) + (t242 + t241 + t240) * qJD(5) / 0.2e1;
t217 = -pkin(7) + t222;
t208 = Ifges(3,4) * t301;
t190 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t301;
t177 = pkin(1) + t303;
t166 = Ifges(3,1) * t302 + Ifges(3,5) * qJD(2) + t208;
t165 = Ifges(3,6) * qJD(2) + qJD(1) * t260;
t150 = Ifges(4,4) * t156;
t148 = -pkin(4) - t151;
t147 = t204 * t305 + t307;
t146 = -t204 * t308 + t306;
t145 = -t204 * t306 + t308;
t144 = t204 * t307 + t305;
t131 = qJD(2) * mrSges(4,1) + t332;
t130 = -qJD(2) * mrSges(4,2) + t333;
t114 = -mrSges(4,1) * t156 - mrSges(4,2) * t157;
t111 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t124;
t110 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t123;
t105 = -t157 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t150;
t104 = t156 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t331;
t83 = -pkin(7) * t158 + t97;
t78 = qJD(4) * t118 + t228 * t158 + t159 * t224;
t73 = -mrSges(5,1) * t270 + mrSges(5,2) * t251;
t64 = t128 + t74;
t63 = t224 * t247 + t228 * t92;
t47 = -mrSges(5,2) * t216 + mrSges(5,3) * t59;
t32 = pkin(4) * t78 - pkin(8) * t77 + t129;
t27 = t223 * t74 + t227 * t49;
t26 = -t223 * t49 + t227 * t74;
t24 = -qJD(4) * t69 + t224 * t248 + t228 * t83;
t23 = t223 * t64 + t227 * t63;
t22 = -t223 * t63 + t227 * t64;
t13 = -mrSges(6,1) * t36 + mrSges(6,2) * t35;
t5 = -qJD(5) * t29 - t223 * t24 + t227 * t32;
t4 = qJD(5) * t28 + t223 * t32 + t227 * t24;
t1 = [t251 * (Ifges(5,1) * t77 - Ifges(5,4) * t78) / 0.2e1 + (-t115 * t159 - t116 * t158 + t170 * t76 - t171 * t75) * mrSges(4,3) + (Ifges(4,1) * t159 - Ifges(4,4) * t158) * t363 + (Ifges(4,5) * t159 - Ifges(4,6) * t158) * t346 + t156 * (Ifges(4,4) * t159 - Ifges(4,2) * t158) / 0.2e1 + t183 * (mrSges(4,1) * t158 + mrSges(4,2) * t159) + (t229 * (-Ifges(3,2) * t225 + t338) + t239) * t294 / 0.2e1 + t270 * (Ifges(5,4) * t77 - Ifges(5,2) * t78) / 0.2e1 + (-t147 * mrSges(6,1) - t146 * mrSges(6,2) + t405 * (t230 * t177 - t226 * t217) + t383 * t226 + (-t380 - t415) * t230) * g(2) + (m(6) * t12 + t13 - t408) * t69 + (t89 * mrSges(5,2) - t16 * mrSges(5,3) + Ifges(5,1) * t58 + Ifges(5,4) * t59 + Ifges(5,5) * t216 + t12 * t262 + t257 * t371 + t259 * t374 + t261 * t375 + t275 * t41) * t118 + (-Ifges(6,1) * t244 - Ifges(6,4) * t245 + Ifges(6,5) * t78) * t366 + t44 * (mrSges(6,1) * t245 - mrSges(6,2) * t244) + t101 * (-Ifges(6,5) * t244 - Ifges(6,6) * t245 + Ifges(6,3) * t78) / 0.2e1 + t90 * (-Ifges(6,4) * t244 - Ifges(6,2) * t245 + Ifges(6,6) * t78) / 0.2e1 + (Ifges(3,4) * t402 + Ifges(3,2) * t403 + Ifges(3,6) * t318 + t166 * t346) * t229 + (Ifges(3,1) * t181 - pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t181) + Ifges(3,4) * t403 + 0.2e1 * Ifges(3,5) * t318) * t225 - t78 * t397 + (-m(5) * t49 + m(6) * t44 + t396) * (qJD(4) * t70 + t224 * t83 - t228 * t248) + (t258 * t346 - t291) * qJD(2) + (t394 / 0.2e1 - mrSges(3,2) * t357 + t395 / 0.2e1 + Ifges(3,6) * t229 / 0.2e1) * qJDD(2) + (t394 + t395) * t318 + t78 * t398 + t338 * t402 + t260 * t403 - t245 * t40 / 0.2e1 + t9 * t315 / 0.2e1 - t165 * t300 / 0.2e1 - t8 * t316 / 0.2e1 - t192 * t317 + m(6) * (t2 * t29 + t20 * t5 + t21 * t4 + t28 * t3) + m(5) * (t125 * t129 + t132 * t89 + t15 * t70 + t24 * t50) + (t180 * t357 + t385) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t385) + t114 * t210 + t78 * t39 / 0.2e1 - t78 * t71 / 0.2e1 + t77 * t72 / 0.2e1 + t4 * t66 + t5 * t67 + t70 * t47 + t28 * t17 + t29 * t18 - t206 * t272 + (-t2 * t316 + t20 * t244 - t21 * t245 - t3 * t315) * mrSges(6,3) + t132 * t278 + Ifges(2,3) * qJDD(1) + t24 * t93 - t78 * t349 - t77 * t350 + m(4) * (t115 * t96 + t116 * t97 + t126 * t75 + t127 * t76 - t143 * t206 + t183 * t210) - t190 * t289 - t246 * t294 + t125 * (mrSges(5,1) * t78 + mrSges(5,2) * t77) + t126 * t111 + t127 * t110 + t129 * t73 + t97 * t130 + t96 * t131 - t158 * t104 / 0.2e1 + t159 * t105 / 0.2e1 + t77 * t281 + t123 * (Ifges(4,4) * t171 + Ifges(4,2) * t170) + t124 * (Ifges(4,1) * t171 + Ifges(4,4) * t170) + t143 * (-mrSges(4,1) * t170 + mrSges(4,2) * t171) - pkin(1) * (-mrSges(3,1) * t180 + mrSges(3,2) * t181) + t218 * (Ifges(5,5) * t77 - Ifges(5,6) * t78) / 0.2e1 - (Ifges(6,5) * t375 + Ifges(6,3) * t371 + Ifges(6,6) * t374 + t292 / 0.2e1 + t89 * mrSges(5,1) - Ifges(5,4) * t58 - Ifges(5,2) * t59 - Ifges(5,6) * t216 - t15 * mrSges(5,3) + t381) * t250 + (-t145 * mrSges(6,1) - t144 * mrSges(6,2) + (-t217 * t405 + t383) * t230 + (m(5) * t177 - m(6) * (-t177 - t388) + t380) * t226) * g(1); -(-Ifges(3,2) * t302 + t166 + t208) * t301 / 0.2e1 - (Ifges(4,2) * t157 + t105 + t150) * t156 / 0.2e1 + ((t220 * t76 + t221 * t75) * pkin(2) - t115 * t121 - t116 * t122 - t183 * t209) * m(4) - t390 * t396 + (t226 * t412 + t268) * g(2) + (t230 * t412 + t267) * g(1) + (-m(5) * t303 - m(6) * (t303 + t388) - m(4) * t215 + t387 + t414) * g(3) + ((m(5) * t50 + m(6) * t255 - t391) * qJD(4) + t408) * t151 + (t257 * t365 + t261 * t367 + t259 * t368 + t20 * t340 + t21 * t341 - t393 / 0.2e1 + t409) * t270 + (Ifges(6,3) * t365 + Ifges(6,5) * t367 + Ifges(6,6) * t368 + t392 / 0.2e1 + t410) * t251 + (m(6) * (t233 + t352) + t382) * (pkin(8) + t152) + (-t125 * t128 + t15 * t152 + t390 * t49 - t50 * t63) * m(5) + t231 + (t165 / 0.2e1 + pkin(6) * t190) * t302 + t104 * t363 + t111 * t359 + t110 * t360 - t3 * t341 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t157 * (Ifges(4,1) * t156 + t331) / 0.2e1 - t116 * t332 + t386 * mrSges(6,3) + t75 * mrSges(4,1) - t76 * mrSges(4,2) - t23 * t66 - t22 * t67 + (t12 * t148 - t20 * t22 - t21 * t23 - t390 * t44) * m(6) + (t291 + (-t239 / 0.2e1 + t246) * qJD(1)) * qJD(1) - t63 * t93 - t114 * t209 - t258 * t294 / 0.2e1 + Ifges(4,6) * t123 + Ifges(4,5) * t124 - t128 * t73 - t122 * t130 - t121 * t131 + t148 * t13 + t152 * t47 - qJD(2) * (Ifges(4,5) * t156 + Ifges(4,6) * t157) / 0.2e1 - t168 * mrSges(3,2) - t169 * mrSges(3,1) + Ifges(3,6) * t180 + Ifges(3,5) * t181 + t115 * t333 - t183 * (-mrSges(4,1) * t157 + mrSges(4,2) * t156); -t156 * t130 - t157 * t131 + t227 * t17 + t223 * t18 - t396 * t251 + t254 * qJD(5) + t391 * t270 + t272 + t278 + (-g(1) * t226 + g(2) * t230) * (m(4) - t405) + (t101 * t255 + t2 * t223 + t227 * t3 - t251 * t44) * m(6) + (t251 * t49 - t270 * t50 + t89) * m(5) + (-t115 * t157 - t116 * t156 + t143) * m(4); t231 + (-t242 / 0.2e1 - t241 / 0.2e1 - t240 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t251 + t256 * mrSges(6,3) + t409) * t270 + t233 * mrSges(6,3) + (t230 * t232 + t267) * g(1) + (t226 * t232 + t268) * g(2) - t396 * t50 + (m(6) * (-t351 + t352 + t386) + t382) * pkin(8) + (-t348 / 0.2e1 - t347 / 0.2e1 - t334 / 0.2e1 + t410) * t251 - t12 * t376 - m(6) * (t20 * t26 + t21 * t27 + t44 * t50) - t27 * t66 - t26 * t67 - pkin(4) * t13 + (t387 - t415) * g(3) - t49 * t93; -t44 * (mrSges(6,1) * t91 + mrSges(6,2) * t90) + (Ifges(6,1) * t90 - t361) * t367 + t40 * t366 + (Ifges(6,5) * t90 - Ifges(6,6) * t91) * t365 - t20 * t66 + t21 * t67 - g(1) * (mrSges(6,1) * t146 - mrSges(6,2) * t147) - g(2) * (-mrSges(6,1) * t144 + mrSges(6,2) * t145) + g(3) * t262 * t203 + (t20 * t90 + t21 * t91) * mrSges(6,3) + t292 + (-Ifges(6,2) * t91 + t41 + t88) * t368 + t381;];
tau = t1;
