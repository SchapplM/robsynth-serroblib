% Calculate vector of inverse dynamics joint torques for
% S5RRPRR5
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:46
% DurationCPUTime: 3.96s
% Computational Cost: add. (5951->373), mult. (9273->494), div. (0->0), fcn. (6441->16), ass. (0->184)
t217 = cos(qJ(2));
t261 = qJD(1) * pkin(1);
t245 = t217 * t261;
t231 = qJD(3) - t245;
t208 = sin(pkin(9));
t212 = sin(qJ(4));
t209 = cos(pkin(9));
t216 = cos(qJ(4));
t251 = t209 * t216;
t155 = -t208 * t212 + t251;
t210 = -pkin(7) - qJ(3);
t167 = t210 * t208;
t197 = t209 * pkin(7);
t168 = qJ(3) * t209 + t197;
t248 = qJD(4) * t216;
t285 = t167 * t248 + qJD(3) * t251 + (-qJD(3) * t208 - qJD(4) * t168) * t212 - t155 * t245;
t113 = t212 * t167 + t216 * t168;
t156 = t208 * t216 + t209 * t212;
t284 = -t113 * qJD(4) - t231 * t156;
t142 = t155 * qJD(4);
t274 = pkin(8) * t142;
t303 = -t274 + t284;
t143 = t156 * qJD(4);
t140 = t143 * pkin(8);
t302 = t140 - t285;
t204 = pkin(9) + qJ(4);
t192 = sin(t204);
t193 = cos(t204);
t301 = -mrSges(5,1) * t193 + mrSges(5,2) * t192;
t230 = -mrSges(4,1) * t209 + mrSges(4,2) * t208;
t199 = qJDD(4) + qJDD(5);
t205 = qJD(4) + qJD(5);
t215 = cos(qJ(5));
t211 = sin(qJ(5));
t206 = qJD(1) + qJD(2);
t131 = t155 * t206;
t213 = sin(qJ(2));
t246 = t213 * t261;
t163 = qJ(3) * t206 + t246;
t241 = pkin(7) * t206 + t163;
t119 = t241 * t208;
t120 = t241 * t209;
t71 = -t119 * t212 + t120 * t216;
t53 = pkin(8) * t131 + t71;
t260 = t211 * t53;
t132 = t156 * t206;
t256 = t120 * t212;
t70 = -t216 * t119 - t256;
t52 = -pkin(8) * t132 + t70;
t50 = qJD(4) * pkin(4) + t52;
t22 = t215 * t50 - t260;
t259 = t215 * t53;
t23 = t211 * t50 + t259;
t234 = t215 * t131 - t132 * t211;
t86 = t131 * t211 + t132 * t215;
t277 = Ifges(6,4) * t86;
t243 = qJD(2) * t261;
t257 = qJDD(1) * pkin(1);
t151 = t213 * t257 + t217 * t243;
t200 = qJDD(1) + qJDD(2);
t118 = qJ(3) * t200 + qJD(3) * t206 + t151;
t242 = pkin(7) * t200 + t118;
t100 = t242 * t208;
t101 = t242 * t209;
t37 = -t71 * qJD(4) - t216 * t100 - t101 * t212;
t91 = t206 * t142 + t156 * t200;
t18 = qJDD(4) * pkin(4) - pkin(8) * t91 + t37;
t36 = -qJD(4) * t256 - t212 * t100 + t216 * t101 - t119 * t248;
t92 = -t206 * t143 + t155 * t200;
t19 = pkin(8) * t92 + t36;
t3 = t22 * qJD(5) + t18 * t211 + t19 * t215;
t34 = t234 * qJD(5) + t211 * t92 + t215 * t91;
t35 = -t86 * qJD(5) - t211 * t91 + t215 * t92;
t4 = -t23 * qJD(5) + t18 * t215 - t19 * t211;
t78 = Ifges(6,4) * t234;
t43 = Ifges(6,1) * t86 + Ifges(6,5) * t205 + t78;
t187 = pkin(3) * t209 + pkin(2);
t129 = -t187 * t206 + t231;
t93 = -pkin(4) * t131 + t129;
t300 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t34 + Ifges(6,6) * t35 + Ifges(6,3) * t199 - (Ifges(6,5) * t234 - Ifges(6,6) * t86) * t205 / 0.2e1 + (t22 * t234 + t23 * t86) * mrSges(6,3) - (-Ifges(6,2) * t86 + t43 + t78) * t234 / 0.2e1 - t93 * (mrSges(6,1) * t86 + mrSges(6,2) * t234) - (Ifges(6,1) * t234 - t277) * t86 / 0.2e1;
t249 = t208 ^ 2 + t209 ^ 2;
t240 = t249 * mrSges(4,3);
t232 = m(4) * qJ(3) + mrSges(5,3) + mrSges(6,3);
t299 = mrSges(4,3) + t232;
t298 = -mrSges(3,2) - m(6) * (-pkin(8) + t210) - m(5) * t210;
t194 = qJ(5) + t204;
t185 = cos(t194);
t174 = t185 * mrSges(6,1);
t297 = m(4) * pkin(2) + m(5) * t187 + m(6) * (pkin(4) * t193 + t187) + mrSges(3,1) + t174 - t230 - t301;
t42 = Ifges(6,2) * t234 + Ifges(6,6) * t205 + t277;
t295 = t42 / 0.2e1;
t112 = t216 * t167 - t168 * t212;
t273 = pkin(8) * t156;
t94 = t112 - t273;
t148 = t155 * pkin(8);
t95 = t148 + t113;
t45 = -t211 * t95 + t215 * t94;
t291 = t45 * qJD(5) + t303 * t211 - t302 * t215;
t46 = t211 * t94 + t215 * t95;
t290 = -t46 * qJD(5) + t302 * t211 + t303 * t215;
t186 = pkin(1) * t213 + qJ(3);
t145 = (-pkin(7) - t186) * t208;
t146 = t186 * t209 + t197;
t103 = t212 * t145 + t216 * t146;
t281 = t86 / 0.2e1;
t279 = t132 / 0.2e1;
t276 = pkin(1) * t217;
t275 = pkin(4) * t143;
t207 = qJ(1) + qJ(2);
t195 = sin(t207);
t272 = g(2) * t195;
t196 = cos(t207);
t271 = g(3) * t196;
t184 = sin(t194);
t266 = mrSges(6,2) * t184;
t265 = mrSges(6,2) * t185;
t263 = Ifges(5,4) * t132;
t262 = pkin(1) * qJD(2);
t258 = mrSges(5,1) * t131 - mrSges(5,2) * t132 - t230 * t206;
t254 = t184 * t195;
t253 = t200 * t208;
t252 = t200 * t209;
t250 = mrSges(6,1) * t254 + t195 * t265;
t247 = -m(4) - m(5) - m(6);
t244 = t213 * t262;
t47 = -t92 * mrSges(5,1) + t91 * mrSges(5,2);
t9 = -t35 * mrSges(6,1) + t34 * mrSges(6,2);
t239 = t174 - t266;
t238 = qJ(3) * t249;
t237 = t249 * t163;
t175 = t217 * t262 + qJD(3);
t236 = t249 * t175;
t235 = t249 * t186;
t102 = t216 * t145 - t146 * t212;
t233 = qJD(3) * t249;
t139 = -mrSges(4,1) * t252 + mrSges(4,2) * t253;
t150 = -t213 * t243 + t217 * t257;
t229 = -mrSges(5,1) * t192 - mrSges(5,2) * t193;
t228 = -mrSges(6,1) * t184 - t265;
t76 = t102 - t273;
t77 = t148 + t103;
t38 = -t211 * t77 + t215 * t76;
t39 = t211 * t76 + t215 * t77;
t105 = t155 * t215 - t156 * t211;
t106 = t155 * t211 + t156 * t215;
t128 = -pkin(4) * t155 - t187;
t227 = mrSges(2,1) + (m(3) - t247) * pkin(1);
t225 = qJDD(3) - t150;
t60 = t145 * t248 + t175 * t251 + (-qJD(4) * t146 - t175 * t208) * t212;
t110 = -t187 * t200 + t225;
t61 = -t103 * qJD(4) - t156 * t175;
t221 = t298 * t195 + (-t266 + t297) * t196;
t220 = -mrSges(6,2) * t254 + t297 * t195 + (-t298 - t299) * t196;
t130 = -pkin(2) * t200 + t225;
t56 = t105 * qJD(5) + t142 * t215 - t143 * t211;
t57 = -t106 * qJD(5) - t142 * t211 - t143 * t215;
t58 = -pkin(4) * t92 + t110;
t81 = t131 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t263;
t127 = Ifges(5,4) * t131;
t82 = t132 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t127;
t219 = t118 * t240 + (Ifges(4,4) * t208 + Ifges(4,2) * t209) * t252 + (Ifges(4,1) * t208 + Ifges(4,4) * t209) * t253 + t130 * t230 + (-t142 * t70 - t71 * t143) * mrSges(5,3) + qJD(4) * (Ifges(5,5) * t142 - Ifges(5,6) * t143) / 0.2e1 + t129 * (mrSges(5,1) * t143 + mrSges(5,2) * t142) + t131 * (Ifges(5,4) * t142 - Ifges(5,2) * t143) / 0.2e1 + (Ifges(5,1) * t142 - Ifges(5,4) * t143) * t279 + t234 * (Ifges(6,4) * t56 + Ifges(6,2) * t57) / 0.2e1 + (-t22 * t56 + t23 * t57) * mrSges(6,3) + (t110 * mrSges(5,2) - t37 * mrSges(5,3) + Ifges(5,1) * t91 + Ifges(5,4) * t92 + Ifges(5,5) * qJDD(4)) * t156 + (-t110 * mrSges(5,1) + t36 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t92 + Ifges(5,6) * qJDD(4)) * t155 + (t58 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t34 + Ifges(6,4) * t35 + Ifges(6,5) * t199) * t106 + (-t58 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t34 + Ifges(6,2) * t35 + Ifges(6,6) * t199) * t105 + t56 * t43 / 0.2e1 + t57 * t295 + t93 * (-mrSges(6,1) * t57 + mrSges(6,2) * t56) + t142 * t82 / 0.2e1 - t143 * t81 / 0.2e1 + t150 * mrSges(3,1) - t151 * mrSges(3,2) + Ifges(3,3) * t200 + t205 * (Ifges(6,5) * t56 + Ifges(6,6) * t57) / 0.2e1 + (Ifges(6,1) * t56 + Ifges(6,4) * t57) * t281;
t218 = cos(qJ(1));
t214 = sin(qJ(1));
t188 = -pkin(2) - t276;
t166 = -t187 - t276;
t157 = -pkin(2) * t206 + t231;
t121 = t244 + t275;
t117 = t128 - t276;
t116 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t132;
t115 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t131;
t80 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t92;
t79 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t91;
t73 = mrSges(6,1) * t205 - mrSges(6,3) * t86;
t72 = -mrSges(6,2) * t205 + mrSges(6,3) * t234;
t49 = t61 - t274;
t48 = -t140 + t60;
t44 = -mrSges(6,1) * t234 + mrSges(6,2) * t86;
t29 = -mrSges(6,2) * t199 + mrSges(6,3) * t35;
t28 = mrSges(6,1) * t199 - mrSges(6,3) * t34;
t25 = t215 * t52 - t260;
t24 = -t211 * t52 - t259;
t8 = -t39 * qJD(5) - t211 * t48 + t215 * t49;
t7 = t38 * qJD(5) + t211 * t49 + t215 * t48;
t1 = [t219 + m(6) * (t117 * t58 + t121 * t93 + t22 * t8 + t23 * t7 + t3 * t39 + t38 * t4) + m(5) * (t102 * t37 + t103 * t36 + t110 * t166 + t129 * t244 + t60 * t71 + t61 * t70) + m(4) * (t118 * t235 + t130 * t188 + t157 * t244 + t163 * t236) + t7 * t72 + t8 * t73 + t38 * t28 + t39 * t29 + (-t214 * mrSges(2,2) + t195 * t299 + t227 * t218 + t221) * g(2) + t102 * t79 + t103 * t80 + t60 * t115 + t61 * t116 + t117 * t9 + t121 * t44 + t166 * t47 + t188 * t139 + Ifges(2,3) * qJDD(1) + (t200 * t235 + t206 * t236) * mrSges(4,3) + (mrSges(2,2) * t218 + t227 * t214 + t220) * g(3) + (m(3) * (t150 * t217 + t151 * t213) + (mrSges(3,1) * t217 - mrSges(3,2) * t213) * t200 + (-mrSges(3,2) * t206 * t217 + (-mrSges(3,1) * t206 - t258) * t213) * qJD(2)) * pkin(1); t219 + ((-t44 + t258) * t213 + (mrSges(3,1) * t213 + (mrSges(3,2) - t240) * t217) * t206) * t261 + t290 * t73 + t291 * t72 + (t200 * t238 + t206 * t233 + t272) * mrSges(4,3) + t45 * t28 + t46 * t29 + (t232 * t195 + t221) * g(2) + t284 * t116 + t220 * g(3) + t285 * t115 + t44 * t275 + t112 * t79 + t113 * t80 + t128 * t9 - pkin(2) * t139 - t187 * t47 + (t128 * t58 + t3 * t46 + t4 * t45 + (-t246 + t275) * t93 + t291 * t23 + t290 * t22) * m(6) + (-t110 * t187 + t112 * t37 + t113 * t36 - t129 * t246 + t284 * t70 + t285 * t71) * m(5) + (-(t157 * t213 + t217 * t237) * t261 - pkin(2) * t130 + t118 * t238 + t163 * t233) * m(4); -t206 ^ 2 * t240 - t131 * t115 + t132 * t116 - t234 * t72 + t86 * t73 + t139 + t47 + t9 + (g(2) * t196 + g(3) * t195) * t247 + (t22 * t86 - t23 * t234 + t58) * m(6) + (-t131 * t71 + t132 * t70 + t110) * m(5) + (-t206 * t237 + t130) * m(4); -(-Ifges(5,2) * t132 + t127 + t82) * t131 / 0.2e1 + t86 * t295 + (-t239 + t301) * g(1) - t132 * (Ifges(5,1) * t131 - t263) / 0.2e1 + (t229 * t195 - t250) * g(2) + (-t132 * t44 + t211 * t29 + t215 * t28 + (-t132 * t93 + t211 * t3 + t215 * t4 - g(1) * t193 + (t271 - t272) * t192) * m(6) + (-t211 * t73 + t215 * t72 + (-t211 * t22 + t215 * t23) * m(6)) * qJD(5)) * pkin(4) + t300 - t25 * t72 - t24 * t73 + t37 * mrSges(5,1) - t36 * mrSges(5,2) - m(6) * (t22 * t24 + t23 * t25) + Ifges(5,5) * t91 + Ifges(5,6) * t92 - t70 * t115 + t71 * t116 - qJD(4) * (Ifges(5,5) * t131 - Ifges(5,6) * t132) / 0.2e1 - t129 * (mrSges(5,1) * t132 + mrSges(5,2) * t131) + (-t228 - t229) * t271 + t81 * t279 + Ifges(5,3) * qJDD(4) + (t131 * t70 + t132 * t71) * mrSges(5,3); -g(1) * t239 - g(2) * t250 - t22 * t72 - t228 * t271 + t23 * t73 + t42 * t281 + t300;];
tau = t1;
