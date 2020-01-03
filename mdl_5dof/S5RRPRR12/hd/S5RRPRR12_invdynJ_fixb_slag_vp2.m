% Calculate vector of inverse dynamics joint torques for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:03
% EndTime: 2019-12-31 20:29:31
% DurationCPUTime: 15.75s
% Computational Cost: add. (4876->560), mult. (10297->745), div. (0->0), fcn. (6445->8), ass. (0->253)
t380 = Ifges(4,4) + Ifges(3,5);
t379 = Ifges(4,6) - Ifges(3,6);
t187 = sin(qJ(5));
t191 = cos(qJ(5));
t135 = -mrSges(6,1) * t191 + mrSges(6,2) * t187;
t201 = m(6) * pkin(4) - t135;
t385 = -mrSges(5,1) - t201;
t190 = sin(qJ(1));
t193 = cos(qJ(2));
t151 = t193 * t190 * qJ(3);
t189 = sin(qJ(2));
t195 = -pkin(2) - pkin(3);
t259 = t195 * t189;
t384 = t190 * t259 + t151;
t278 = qJD(1) * t189;
t167 = pkin(6) * t278;
t122 = pkin(7) * t278 - t167;
t383 = -qJD(3) + t122;
t180 = -qJDD(2) + qJDD(4);
t268 = qJD(1) * qJD(2);
t124 = -t193 * qJDD(1) + t189 * t268;
t125 = qJDD(1) * t189 + t193 * t268;
t188 = sin(qJ(4));
t192 = cos(qJ(4));
t111 = t188 * t189 + t192 * t193;
t202 = t111 * qJD(4);
t41 = -qJD(1) * t202 + t124 * t188 + t125 * t192;
t34 = mrSges(5,1) * t180 - mrSges(5,3) * t41;
t181 = -qJD(2) + qJD(4);
t277 = qJD(1) * t193;
t99 = -t188 * t277 + t192 * t278;
t67 = t181 * t191 - t187 * t99;
t21 = qJD(5) * t67 + t180 * t187 + t191 * t41;
t68 = t181 * t187 + t191 * t99;
t22 = -qJD(5) * t68 + t180 * t191 - t187 * t41;
t7 = -mrSges(6,1) * t22 + mrSges(6,2) * t21;
t382 = t34 - t7;
t97 = -t188 * t278 - t192 * t277;
t88 = Ifges(5,4) * t97;
t381 = t97 * Ifges(5,2);
t309 = t99 * mrSges(5,3);
t361 = -mrSges(5,1) * t181 - mrSges(6,1) * t67 + mrSges(6,2) * t68 + t309;
t165 = Ifges(3,4) * t277;
t302 = Ifges(4,5) * t193;
t227 = Ifges(4,1) * t189 - t302;
t378 = Ifges(3,1) * t278 + qJD(1) * t227 + qJD(2) * t380 + t165;
t377 = Ifges(5,5) * t181;
t376 = t181 * Ifges(5,6);
t228 = mrSges(6,1) * t187 + mrSges(6,2) * t191;
t184 = qJD(2) * qJ(3);
t168 = pkin(6) * t277;
t236 = -pkin(7) * t277 + t168;
t100 = t184 + t236;
t252 = t195 * qJD(2);
t77 = t252 - t383;
t48 = -t100 * t188 + t192 * t77;
t45 = -pkin(4) * t181 - t48;
t375 = t228 * t45;
t256 = mrSges(4,2) * t278;
t373 = mrSges(3,3) * t278 + t256 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t372 = t189 * t379 + t193 * t380;
t229 = t193 * mrSges(4,1) + t189 * mrSges(4,3);
t231 = mrSges(3,1) * t193 - mrSges(3,2) * t189;
t371 = t229 + t231;
t212 = t188 * t193 - t189 * t192;
t239 = t111 * mrSges(5,1) - mrSges(5,2) * t212;
t370 = t111 * t201 + t239;
t270 = qJD(5) * t191;
t63 = qJD(2) * t111 - t202;
t207 = t187 * t63 - t212 * t270;
t107 = t124 * pkin(6);
t108 = t125 * pkin(6);
t369 = -t107 * t193 + t108 * t189;
t76 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t107;
t244 = qJDD(3) + t108;
t80 = -qJDD(2) * pkin(2) + t244;
t368 = t189 * t80 + t193 * t76;
t367 = -mrSges(2,1) - t371;
t260 = m(6) * pkin(8) + mrSges(6,3);
t366 = mrSges(5,3) - mrSges(3,3) + mrSges(2,2) - mrSges(4,2);
t365 = -mrSges(5,2) + t260;
t272 = qJD(4) * t192;
t273 = qJD(4) * t188;
t56 = -pkin(7) * t125 + qJDD(2) * t195 + t244;
t59 = pkin(7) * t124 + t76;
t13 = -t100 * t273 + t188 * t56 + t192 * t59 + t272 * t77;
t11 = pkin(8) * t180 + t13;
t101 = -qJD(1) * pkin(1) - pkin(2) * t277 - qJ(3) * t278;
t75 = pkin(3) * t277 - t101;
t36 = -pkin(4) * t97 - pkin(8) * t99 + t75;
t49 = t100 * t192 + t188 * t77;
t46 = pkin(8) * t181 + t49;
t15 = -t187 * t46 + t191 * t36;
t172 = t189 * qJD(3);
t185 = qJDD(1) * pkin(1);
t54 = pkin(2) * t124 - qJ(3) * t125 - qJD(1) * t172 - t185;
t39 = -pkin(3) * t124 - t54;
t42 = qJD(1) * qJD(4) * t212 + t124 * t192 - t125 * t188;
t8 = -pkin(4) * t42 - pkin(8) * t41 + t39;
t1 = qJD(5) * t15 + t11 * t191 + t187 * t8;
t16 = t187 * t36 + t191 * t46;
t2 = -qJD(5) * t16 - t11 * t187 + t191 * t8;
t364 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t363 = m(5) * t48 - m(6) * t45 - t361;
t14 = -qJD(4) * t49 - t188 * t59 + t192 * t56;
t12 = -pkin(4) * t180 - t14;
t221 = Ifges(6,5) * t191 - Ifges(6,6) * t187;
t304 = Ifges(6,4) * t191;
t223 = -Ifges(6,2) * t187 + t304;
t305 = Ifges(6,4) * t187;
t226 = Ifges(6,1) * t191 - t305;
t271 = qJD(5) * t187;
t245 = -t271 / 0.2e1;
t64 = Ifges(6,4) * t67;
t91 = qJD(5) - t97;
t30 = Ifges(6,1) * t68 + Ifges(6,5) * t91 + t64;
t325 = t191 / 0.2e1;
t254 = t30 * t325;
t28 = Ifges(6,5) * t68 + Ifges(6,6) * t67 + Ifges(6,3) * t91;
t323 = Ifges(6,4) * t68;
t29 = Ifges(6,2) * t67 + Ifges(6,6) * t91 + t323;
t294 = t191 * t97;
t298 = t187 * t97;
t40 = qJDD(5) - t42;
t3 = Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t40;
t310 = t97 * mrSges(5,3);
t329 = t99 / 0.2e1;
t332 = t91 / 0.2e1;
t333 = -t91 / 0.2e1;
t334 = t68 / 0.2e1;
t335 = -t68 / 0.2e1;
t336 = t67 / 0.2e1;
t337 = -t67 / 0.2e1;
t338 = t40 / 0.2e1;
t340 = t22 / 0.2e1;
t341 = t21 / 0.2e1;
t342 = Ifges(6,1) * t341 + Ifges(6,4) * t340 + Ifges(6,5) * t338;
t354 = t1 * t191 - t187 * t2;
t89 = Ifges(5,4) * t99;
t50 = t376 + t89 + t381;
t51 = Ifges(5,1) * t99 + t377 + t88;
t362 = (Ifges(6,5) * t99 + t226 * t97) * t335 + (Ifges(6,6) * t99 + t223 * t97) * t337 + (Ifges(6,3) * t99 + t221 * t97) * t333 + (-t15 * t270 - t16 * t271 + t354) * mrSges(6,3) + t16 * (mrSges(6,2) * t99 + mrSges(6,3) * t298) + t15 * (-mrSges(6,1) * t99 + mrSges(6,3) * t294) - (Ifges(5,1) * t97 + t28 - t89) * t99 / 0.2e1 - (-Ifges(5,2) * t99 + t51 + t88) * t97 / 0.2e1 + (Ifges(6,5) * t187 + Ifges(6,6) * t191) * t338 + (Ifges(6,2) * t191 + t305) * t340 + (Ifges(6,1) * t187 + t304) * t341 + Ifges(5,5) * t41 + Ifges(5,6) * t42 - t13 * mrSges(5,2) + t14 * mrSges(5,1) + (t298 / 0.2e1 + t245) * t29 - t30 * t294 / 0.2e1 + t187 * t342 - t75 * (mrSges(5,1) * t99 + mrSges(5,2) * t97) - t97 * t375 + (t221 * t332 + t223 * t336 + t226 * t334 + t254 + t375) * qJD(5) + t12 * t135 + Ifges(5,3) * t180 + t49 * t309 + t48 * t310 + t3 * t325 + t50 * t329 - t181 * (Ifges(5,5) * t97 - Ifges(5,6) * t99) / 0.2e1;
t127 = qJ(3) * t192 + t188 * t195;
t194 = cos(qJ(1));
t178 = t194 * pkin(6);
t173 = t189 * qJ(3);
t243 = -pkin(1) - t173;
t356 = t190 * (t193 * t195 + t243) - pkin(7) * t194 + t178;
t10 = -mrSges(6,2) * t40 + mrSges(6,3) * t22;
t9 = mrSges(6,1) * t40 - mrSges(6,3) * t21;
t355 = t10 * t191 - t187 * t9;
t353 = g(1) * t194 + g(2) * t190;
t352 = qJD(3) + t167;
t43 = -mrSges(6,2) * t91 + mrSges(6,3) * t67;
t44 = mrSges(6,1) * t91 - mrSges(6,3) * t68;
t72 = -mrSges(5,2) * t181 + t310;
t351 = -t187 * t44 + t191 * t43 + t72;
t282 = t193 * t194;
t284 = t189 * t194;
t84 = -t188 * t284 - t192 * t282;
t85 = t188 * t282 - t192 * t284;
t350 = t84 * mrSges(5,2) + t385 * t85;
t82 = t212 * t190;
t83 = t111 * t190;
t349 = -t83 * mrSges(5,2) + t385 * t82;
t197 = (-t15 * t191 - t16 * t187) * qJD(5) + t354;
t345 = m(6) * t197 - t270 * t44 - t271 * t43 + t355;
t328 = pkin(6) - pkin(7);
t326 = t189 / 0.2e1;
t318 = pkin(6) * t189;
t317 = pkin(6) * t193;
t177 = t193 * pkin(2);
t307 = Ifges(3,4) * t189;
t306 = Ifges(3,4) * t193;
t303 = Ifges(4,5) * t189;
t293 = t193 * mrSges(4,3);
t288 = t212 * t187;
t287 = t212 * t191;
t274 = qJD(2) * t193;
t281 = qJ(3) * t274 + t172;
t280 = t177 + t173;
t279 = pkin(1) * t194 + pkin(6) * t190;
t276 = qJD(2) * t189;
t275 = qJD(2) * t192;
t269 = m(4) + m(5) + m(6);
t263 = Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t40;
t141 = t328 * t193;
t255 = mrSges(4,2) * t277;
t253 = pkin(3) * t193 + t280;
t241 = -t268 / 0.2e1;
t86 = -qJDD(2) * mrSges(4,1) + mrSges(4,2) * t125;
t238 = pkin(2) * t282 + qJ(3) * t284 + t279;
t237 = qJD(2) * t141;
t58 = pkin(4) * t99 - pkin(8) * t97;
t230 = mrSges(3,1) * t189 + mrSges(3,2) * t193;
t225 = t193 * Ifges(3,2) + t307;
t220 = -t15 * t187 + t16 * t191;
t211 = pkin(8) * t212 + t253;
t47 = pkin(4) * t111 + pkin(1) + t211;
t140 = t328 * t189;
t70 = t140 * t188 + t141 * t192;
t26 = -t187 * t70 + t191 * t47;
t27 = t187 * t47 + t191 * t70;
t217 = t187 * t83 - t191 * t194;
t216 = -t187 * t194 - t191 * t83;
t126 = -qJ(3) * t188 + t192 * t195;
t128 = -qJD(2) * pkin(2) + t352;
t132 = t168 + t184;
t214 = t128 * t193 - t132 * t189;
t213 = t140 * t192 - t141 * t188;
t159 = qJ(3) * t277;
t81 = qJD(1) * t259 + t159;
t208 = pkin(1) * t230;
t206 = -t191 * t63 - t212 * t271;
t205 = t101 * (t189 * mrSges(4,1) - t293);
t204 = t189 * (Ifges(3,1) * t193 - t307);
t203 = t193 * (Ifges(4,3) * t189 + t302);
t74 = t189 * t252 + t281;
t200 = pkin(3) * t282 - pkin(7) * t190 + t238;
t199 = t293 + (-m(4) * pkin(2) - mrSges(4,1)) * t189;
t164 = Ifges(4,5) * t278;
t152 = qJ(3) * t282;
t134 = qJD(2) * mrSges(4,3) + t255;
t133 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t277;
t129 = -pkin(1) - t280;
t123 = t328 * t276;
t120 = pkin(4) - t126;
t115 = pkin(2) * t278 - t159;
t114 = t229 * qJD(1);
t104 = pkin(1) + t253;
t98 = t187 * t278 + t191 * t275;
t96 = -t187 * t275 + t191 * t278;
t93 = Ifges(3,6) * qJD(2) + qJD(1) * t225;
t92 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t277 + t164;
t90 = pkin(2) * t276 - t281;
t87 = -mrSges(4,2) * t124 + qJDD(2) * mrSges(4,3);
t66 = -t187 * t190 - t191 * t84;
t65 = t187 * t84 - t190 * t191;
t62 = t188 * t274 - t193 * t273 + (t272 - t275) * t189;
t61 = t122 * t192 + t188 * t236;
t57 = t228 * t212;
t55 = -mrSges(5,1) * t97 + mrSges(5,2) * t99;
t37 = -t58 + t81;
t35 = -mrSges(5,2) * t180 + mrSges(5,3) * t42;
t31 = qJD(4) * t213 - t123 * t192 + t188 * t237;
t25 = t187 * t58 + t191 * t48;
t24 = -t187 * t48 + t191 * t58;
t23 = pkin(4) * t62 - pkin(8) * t63 + t74;
t18 = t187 * t37 + t191 * t61;
t17 = -t187 * t61 + t191 * t37;
t6 = -qJD(5) * t27 - t187 * t31 + t191 * t23;
t5 = qJD(5) * t26 + t187 * t23 + t191 * t31;
t4 = [(-t93 / 0.2e1 + t92 / 0.2e1) * t276 - t363 * (qJD(4) * t70 - t123 * t188 - t192 * t237) + m(6) * (t1 * t27 + t15 * t6 + t16 * t5 + t2 * t26) + (t14 * mrSges(5,3) - Ifges(5,1) * t41 - Ifges(5,4) * t42 - Ifges(5,5) * t180 - t221 * t338 - t223 * t340 - t226 * t341 - t245 * t30) * t212 + t231 * t185 + (t303 / 0.2e1 - pkin(1) * mrSges(3,1) + t129 * mrSges(4,1) - mrSges(3,3) * t317 - t225 / 0.2e1 + (Ifges(4,5) - Ifges(3,4)) * t326 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t193) * t124 + (t1 * t288 + t15 * t206 - t16 * t207 + t2 * t287) * mrSges(6,3) + (t189 * (Ifges(4,1) * t193 + t303) + t193 * (-Ifges(3,2) * t189 + t306) + t204) * t268 / 0.2e1 + (t189 * Ifges(3,1) + t227 + t306) * t125 / 0.2e1 + (-Ifges(6,5) * t206 - Ifges(6,6) * t207) * t332 + m(5) * (t104 * t39 + t13 * t70 + t31 * t49 + t74 * t75) + (-pkin(1) * t125 - qJDD(2) * t317) * mrSges(3,2) + (-Ifges(6,1) * t206 - Ifges(6,4) * t207) * t334 - t193 * (Ifges(4,5) * t125 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (-t48 * t63 - t49 * t62) * mrSges(5,3) - t12 * t57 + t5 * t43 + t6 * t44 - t208 * t268 + t26 * t9 + t27 * t10 + (-Ifges(6,4) * t206 - Ifges(6,2) * t207) * t336 + (t373 * t274 + (-t134 - t133) * t276 + m(4) * (qJD(2) * t214 + t368) + m(3) * t369) * pkin(6) + (Ifges(6,3) * t338 + Ifges(6,6) * t340 + Ifges(6,5) * t341 - Ifges(5,4) * t41 - Ifges(5,2) * t42 - Ifges(5,6) * t180 + t263 / 0.2e1 - t13 * mrSges(5,3) + t364) * t111 + t193 * (Ifges(3,4) * t125 + Ifges(3,6) * qJDD(2)) / 0.2e1 + m(4) * (t101 * t90 + t129 * t54) + (-qJDD(2) * mrSges(3,1) + t86) * t318 + t3 * t288 / 0.2e1 - t287 * t342 + t39 * t239 + t70 * t35 + t31 * t72 + t74 * t55 + (-t216 * mrSges(6,1) - t217 * mrSges(6,2) - (-pkin(4) * t83 + t356) * m(6) - t356 * m(5) + t83 * mrSges(5,1) + t365 * t82 + (-m(4) - m(3)) * t178 + (-m(4) * (t243 - t177) + m(3) * pkin(1) - t367) * t190 + t366 * t194) * g(1) + (-m(6) * (-pkin(4) * t84 + t200) - t66 * mrSges(6,1) - t65 * mrSges(6,2) - m(5) * t200 + t84 * mrSges(5,1) - m(3) * t279 - m(4) * t238 - t365 * t85 + t367 * t194 + t366 * t190) * g(2) + (t128 * t274 - t132 * t276 + t368) * mrSges(4,2) + (t125 * t318 + t369) * mrSges(3,3) - t207 * t29 / 0.2e1 + t104 * (-mrSges(5,1) * t42 + mrSges(5,2) * t41) + t378 * t274 / 0.2e1 + t45 * (mrSges(6,1) * t207 - mrSges(6,2) * t206) - t90 * t114 + ((Ifges(3,1) + Ifges(4,1)) * t125 + t380 * qJDD(2)) * t326 + (t380 * t189 - t379 * t193) * qJDD(2) / 0.2e1 + (t28 / 0.2e1 - t50 / 0.2e1 + Ifges(6,6) * t336 + t75 * mrSges(5,1) - t381 / 0.2e1 + t15 * mrSges(6,1) - t16 * mrSges(6,2) - Ifges(5,4) * t329 + Ifges(6,3) * t332 + Ifges(6,5) * t334 - t376 / 0.2e1) * t62 + (t51 / 0.2e1 + t75 * mrSges(5,2) + t88 / 0.2e1 + t254 + Ifges(5,1) * t329 + t377 / 0.2e1) * t63 + t203 * t241 + (m(5) * t14 - m(6) * t12 + t382) * t213 + t87 * t317 + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * qJDD(1) - t54 * t229 + (t205 + t372 * qJD(2) / 0.2e1) * qJD(2) - t129 * mrSges(4,3) * t125; (-m(5) * (t194 * t259 + t152) - m(4) * t152 - t199 * t194 - m(6) * (pkin(8) * t84 + t195 * t284 + t152) - t84 * mrSges(6,3) + t350) * g(1) + t345 * (-pkin(8) + t127) + (m(5) * t49 + m(6) * t220 + t351) * (qJD(3) * t192 + qJD(4) * t126) + t352 * t134 + t353 * t230 + (-m(5) * t384 - m(4) * t151 - t199 * t190 - m(6) * (-pkin(8) * t83 + t384) + t83 * mrSges(6,3) + t349) * g(2) - (Ifges(4,1) * t277 + t164 + t92) * t278 / 0.2e1 + t133 * t167 + t132 * t256 + t93 * t278 / 0.2e1 - t18 * t43 - t17 * t44 + (t126 * t14 + t127 * t13 - t49 * t61 - t75 * t81) * m(5) + (t12 * t120 - t15 * t17 - t16 * t18) * m(6) + (-pkin(2) * t80 + qJ(3) * t76 + qJD(3) * t132 - t101 * t115) * m(4) + (-pkin(6) * t214 * m(4) - t205 + (t208 + t203 / 0.2e1 - t204 / 0.2e1) * qJD(1)) * qJD(1) - t362 - t128 * t255 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) + t363 * (-qJD(4) * t127 + t188 * t383 - t192 * t236) - t61 * t72 + t76 * mrSges(4,3) - t80 * mrSges(4,1) - t81 * t55 + (-m(4) * t280 - m(5) * t253 - m(6) * t211 - mrSges(6,3) * t212 - t370 - t371) * g(3) + t372 * t241 - t373 * t168 - pkin(2) * t86 + qJ(3) * t87 + t107 * mrSges(3,2) - t108 * mrSges(3,1) - (-Ifges(3,2) * t278 + t165 + t378) * t277 / 0.2e1 + t379 * t124 + t115 * t114 + t380 * t125 + t120 * t7 + t126 * t34 + t127 * t35; -qJD(2) * t134 - t98 * t43 - t96 * t44 + t269 * t193 * g(3) + (-qJD(2) * t72 + qJD(4) * t351 + t382) * t192 + ((-t114 - t55) * qJD(1) - t353 * t269) * t189 + (t35 + (-t187 * t43 - t191 * t44) * qJD(5) + t181 * t361 + t355) * t188 + t86 + (-t15 * t96 - t16 * t98 + (qJD(4) * t220 - t12) * t192 + (t181 * t45 + t197) * t188) * m(6) + (t13 * t188 + t14 * t192 - t278 * t75 + t181 * (-t188 * t48 + t192 * t49)) * m(5) + (-qJD(2) * t132 + t101 * t278 + t80) * m(4); (t212 * t260 + t370) * g(3) - t25 * t43 - t24 * t44 - pkin(4) * t7 + t345 * pkin(8) + (t260 * t84 - t350) * g(1) + (-t260 * t83 - t349) * g(2) + (-pkin(4) * t12 - t15 * t24 - t16 * t25 - t45 * t49) * m(6) - t361 * t49 - t48 * t72 + t362; -t45 * (mrSges(6,1) * t68 + mrSges(6,2) * t67) + (Ifges(6,1) * t67 - t323) * t335 + t29 * t334 + (Ifges(6,5) * t67 - Ifges(6,6) * t68) * t333 - t15 * t43 + t16 * t44 - g(1) * (mrSges(6,1) * t65 - mrSges(6,2) * t66) - g(2) * (-mrSges(6,1) * t217 + mrSges(6,2) * t216) - g(3) * t57 + (t15 * t67 + t16 * t68) * mrSges(6,3) + t263 + (-Ifges(6,2) * t68 + t30 + t64) * t337 + t364;];
tau = t4;
