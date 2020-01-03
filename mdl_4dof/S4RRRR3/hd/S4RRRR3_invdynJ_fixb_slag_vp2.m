% Calculate vector of inverse dynamics joint torques for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:26
% DurationCPUTime: 4.48s
% Computational Cost: add. (3289->361), mult. (7833->508), div. (0->0), fcn. (5278->12), ass. (0->174)
t161 = cos(qJ(3));
t142 = pkin(2) * t161 + pkin(3);
t156 = sin(qJ(4));
t157 = sin(qJ(3));
t160 = cos(qJ(4));
t193 = qJD(4) * t160;
t194 = qJD(4) * t156;
t200 = t156 * t157;
t158 = sin(qJ(2));
t162 = cos(qJ(2));
t117 = -t157 * t158 + t161 * t162;
t105 = t117 * qJD(1);
t221 = pkin(7) * t105;
t164 = -pkin(6) - pkin(5);
t135 = t164 * t162;
t124 = qJD(1) * t135;
t110 = t161 * t124;
t134 = t164 * t158;
t123 = qJD(1) * t134;
t81 = -t123 * t157 + t110;
t56 = t81 - t221;
t118 = t157 * t162 + t158 * t161;
t106 = t118 * qJD(1);
t100 = t106 * pkin(7);
t107 = t157 * t124;
t82 = t161 * t123 + t107;
t57 = -t100 + t82;
t244 = -t156 * t56 - t160 * t57 + t142 * t193 + (-t157 * t194 + (t160 * t161 - t200) * qJD(3)) * pkin(2);
t199 = t157 * t160;
t243 = t156 * t57 - t160 * t56 - t142 * t194 + (-t157 * t193 + (-t156 * t161 - t199) * qJD(3)) * pkin(2);
t132 = -mrSges(3,1) * t162 + mrSges(3,2) * t158;
t155 = qJ(2) + qJ(3);
t149 = sin(t155);
t150 = cos(t155);
t151 = qJ(4) + t155;
t140 = sin(t151);
t141 = cos(t151);
t181 = t141 * mrSges(5,1) - t140 * mrSges(5,2);
t170 = -t150 * mrSges(4,1) + mrSges(4,2) * t149 - t181;
t242 = t132 + t170;
t190 = qJD(1) * qJD(2);
t127 = qJDD(1) * t162 - t158 * t190;
t159 = sin(qJ(1));
t163 = cos(qJ(1));
t234 = g(1) * t163 + g(2) * t159;
t180 = t160 * t105 - t106 * t156;
t128 = qJDD(1) * t158 + t162 * t190;
t172 = t117 * qJD(3);
t54 = qJD(1) * t172 + t127 * t157 + t128 * t161;
t173 = t118 * qJD(3);
t55 = -qJD(1) * t173 + t127 * t161 - t128 * t157;
t13 = t180 * qJD(4) + t156 * t55 + t160 * t54;
t69 = t105 * t156 + t106 * t160;
t14 = -qJD(4) * t69 - t156 * t54 + t160 * t55;
t152 = qJDD(2) + qJDD(3);
t147 = qJDD(4) + t152;
t113 = qJD(2) * pkin(2) + t123;
t195 = qJD(3) * t161;
t196 = qJD(3) * t157;
t116 = t128 * pkin(5);
t86 = qJDD(2) * pkin(2) - pkin(6) * t128 - t116;
t115 = t127 * pkin(5);
t87 = pkin(6) * t127 + t115;
t28 = t113 * t195 + t124 * t196 + t157 * t86 + t161 * t87;
t10 = pkin(7) * t55 + t28;
t76 = t113 * t157 - t110;
t51 = t76 + t221;
t207 = t156 * t51;
t153 = qJD(2) + qJD(3);
t75 = t161 * t113 + t107;
t50 = -t100 + t75;
t46 = pkin(3) * t153 + t50;
t19 = t160 * t46 - t207;
t29 = -qJD(3) * t76 - t157 * t87 + t161 * t86;
t6 = pkin(3) * t152 - pkin(7) * t54 + t29;
t2 = qJD(4) * t19 + t10 * t160 + t156 * t6;
t206 = t160 * t51;
t20 = t156 * t46 + t206;
t3 = -qJD(4) * t20 - t10 * t156 + t160 * t6;
t148 = qJD(4) + t153;
t63 = Ifges(5,4) * t180;
t34 = Ifges(5,1) * t69 + Ifges(5,5) * t148 + t63;
t241 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t13 + Ifges(5,6) * t14 + Ifges(5,3) * t147 - (-Ifges(5,2) * t69 + t34 + t63) * t180 / 0.2e1;
t223 = Ifges(5,4) * t69;
t33 = Ifges(5,2) * t180 + Ifges(5,6) * t148 + t223;
t240 = t33 / 0.2e1;
t225 = t158 / 0.2e1;
t239 = mrSges(5,2) * t180;
t238 = Ifges(5,1) * t180;
t237 = Ifges(5,5) * t180;
t204 = qJDD(1) * pkin(1);
t90 = t157 * t134 - t161 * t135;
t236 = t115 * t162 + t116 * t158;
t139 = pkin(3) * t150;
t143 = pkin(2) * t162 + pkin(1);
t233 = mrSges(2,1) + m(5) * (t139 + t143) + m(4) * t143 + m(3) * pkin(1) - t242;
t232 = mrSges(2,2) + m(5) * (-pkin(7) + t164) - mrSges(5,3) + m(4) * t164 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t230 = -t69 / 0.2e1;
t229 = t69 / 0.2e1;
t227 = t106 / 0.2e1;
t226 = -t148 / 0.2e1;
t222 = pkin(2) * t158;
t218 = g(3) * t162;
t217 = t19 * mrSges(5,3);
t216 = t20 * mrSges(5,3);
t215 = mrSges(3,2) * t162;
t213 = mrSges(4,3) * t105;
t212 = Ifges(3,4) * t158;
t211 = Ifges(3,4) * t162;
t210 = t106 * mrSges(4,3);
t209 = t106 * Ifges(4,4);
t205 = t162 * Ifges(3,2);
t198 = qJD(2) * t158;
t197 = qJD(2) * t162;
t192 = t158 * qJD(1);
t191 = t162 * qJD(1);
t188 = pkin(2) * t198;
t186 = qJD(2) * t164;
t89 = t161 * t134 + t135 * t157;
t133 = t143 * qJD(1);
t179 = mrSges(5,1) * t140 + mrSges(5,2) * t141;
t178 = t205 + t212;
t177 = Ifges(3,5) * t162 - Ifges(3,6) * t158;
t98 = -pkin(2) * t127 - t204;
t64 = -pkin(7) * t118 + t89;
t65 = pkin(7) * t117 + t90;
t35 = -t156 * t65 + t160 * t64;
t36 = t156 * t64 + t160 * t65;
t79 = t117 * t160 - t118 * t156;
t80 = t117 * t156 + t118 * t160;
t175 = pkin(1) * (mrSges(3,1) * t158 + t215);
t174 = t158 * (Ifges(3,1) * t162 - t212);
t125 = t158 * t186;
t126 = t162 * t186;
t40 = t161 * t125 + t157 * t126 + t134 * t195 + t135 * t196;
t167 = mrSges(4,1) * t149 + mrSges(4,2) * t150 + t179;
t41 = -t90 * qJD(3) - t125 * t157 + t161 * t126;
t61 = t105 * Ifges(4,2) + t153 * Ifges(4,6) + t209;
t99 = Ifges(4,4) * t105;
t62 = t106 * Ifges(4,1) + t153 * Ifges(4,5) + t99;
t88 = -pkin(3) * t105 - t133;
t165 = -t28 * mrSges(4,2) + t133 * (mrSges(4,1) * t106 + mrSges(4,2) * t105) + t180 * t217 - t88 * t239 + Ifges(4,3) * t152 + t238 * t230 + t29 * mrSges(4,1) + t237 * t226 + t61 * t227 - t106 * (Ifges(4,1) * t105 - t209) / 0.2e1 + t75 * t213 + Ifges(4,6) * t55 + Ifges(4,5) * t54 - t153 * (Ifges(4,5) * t105 - Ifges(4,6) * t106) / 0.2e1 + (-t88 * mrSges(5,1) - Ifges(5,4) * t230 - Ifges(5,6) * t226 + t216 + t240) * t69 - (-Ifges(4,2) * t106 + t62 + t99) * t105 / 0.2e1 + t241;
t145 = Ifges(3,4) * t191;
t131 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t191;
t130 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t192;
t104 = Ifges(3,1) * t192 + Ifges(3,5) * qJD(2) + t145;
t103 = Ifges(3,6) * qJD(2) + qJD(1) * t178;
t102 = pkin(2) * t199 + t142 * t156;
t101 = -pkin(2) * t200 + t142 * t160;
t94 = -pkin(3) * t117 - t143;
t93 = mrSges(4,1) * t153 - t210;
t92 = -mrSges(4,2) * t153 + t213;
t91 = pkin(2) * t192 + pkin(3) * t106;
t84 = -qJD(2) * t118 - t173;
t83 = qJD(2) * t117 + t172;
t74 = -mrSges(4,1) * t105 + mrSges(4,2) * t106;
t70 = -pkin(3) * t84 + t188;
t59 = mrSges(5,1) * t148 - mrSges(5,3) * t69;
t58 = -mrSges(5,2) * t148 + mrSges(5,3) * t180;
t45 = -mrSges(4,2) * t152 + mrSges(4,3) * t55;
t44 = mrSges(4,1) * t152 - mrSges(4,3) * t54;
t39 = -pkin(3) * t55 + t98;
t38 = -mrSges(5,1) * t180 + mrSges(5,2) * t69;
t31 = -pkin(7) * t83 + t41;
t30 = pkin(7) * t84 + t40;
t26 = -qJD(4) * t80 - t156 * t83 + t160 * t84;
t25 = qJD(4) * t79 + t156 * t84 + t160 * t83;
t22 = t160 * t50 - t207;
t21 = -t156 * t50 - t206;
t9 = -mrSges(5,2) * t147 + mrSges(5,3) * t14;
t8 = mrSges(5,1) * t147 - mrSges(5,3) * t13;
t5 = -qJD(4) * t36 - t156 * t30 + t160 * t31;
t4 = qJD(4) * t35 + t156 * t31 + t160 * t30;
t1 = [(-mrSges(5,1) * t39 + mrSges(5,3) * t2 + Ifges(5,4) * t13 + Ifges(5,2) * t14 + Ifges(5,6) * t147) * t79 + (mrSges(4,2) * t98 - mrSges(4,3) * t29 + Ifges(4,1) * t54 + Ifges(4,4) * t55 + Ifges(4,5) * t152) * t118 + (-mrSges(4,1) * t98 + mrSges(4,3) * t28 + Ifges(4,4) * t54 + Ifges(4,2) * t55 + Ifges(4,6) * t152) * t117 + (-t75 * t83 + t76 * t84) * mrSges(4,3) + t153 * (Ifges(4,5) * t83 + Ifges(4,6) * t84) / 0.2e1 - t143 * (-mrSges(4,1) * t55 + mrSges(4,2) * t54) + t148 * (Ifges(5,5) * t25 + Ifges(5,6) * t26) / 0.2e1 - t133 * (-mrSges(4,1) * t84 + mrSges(4,2) * t83) + t105 * (Ifges(4,4) * t83 + Ifges(4,2) * t84) / 0.2e1 + (t162 * (-Ifges(3,2) * t158 + t211) + t174) * t190 / 0.2e1 + t88 * (-mrSges(5,1) * t26 + mrSges(5,2) * t25) + t89 * t44 + t90 * t45 + t180 * (Ifges(5,4) * t25 + Ifges(5,2) * t26) / 0.2e1 + t40 * t92 + t41 * t93 + t94 * (-mrSges(5,1) * t14 + mrSges(5,2) * t13) + m(4) * (-t133 * t188 - t143 * t98 + t28 * t90 + t29 * t89 + t40 * t76 + t41 * t75) + t83 * t62 / 0.2e1 + t84 * t61 / 0.2e1 + t70 * t38 + t4 * t58 + t5 * t59 + t25 * t34 / 0.2e1 + t35 * t8 + t36 * t9 + m(5) * (t19 * t5 + t2 * t36 + t20 * t4 + t3 * t35 + t39 * t94 + t70 * t88) + t74 * t188 + t162 * (Ifges(3,4) * t128 + Ifges(3,2) * t127) / 0.2e1 + (m(3) * t204 + mrSges(3,1) * t127 - mrSges(3,2) * t128) * pkin(1) + (mrSges(5,2) * t39 - mrSges(5,3) * t3 + Ifges(5,1) * t13 + Ifges(5,4) * t14 + Ifges(5,5) * t147) * t80 + (0.2e1 * Ifges(3,5) * t225 + Ifges(3,6) * t162) * qJDD(2) + (Ifges(4,1) * t83 + Ifges(4,4) * t84) * t227 + (Ifges(5,1) * t25 + Ifges(5,4) * t26) * t229 + t26 * t216 + qJD(2) ^ 2 * t177 / 0.2e1 + t127 * t178 / 0.2e1 + (t162 * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t127) - t158 * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t128) - t130 * t197 - t131 * t198 + m(3) * t236) * pkin(5) + t236 * mrSges(3,3) + (t233 * t159 + t232 * t163) * g(1) + (t232 * t159 - t233 * t163) * g(2) + (Ifges(3,1) * t128 + Ifges(3,4) * t127) * t225 + t26 * t240 - t175 * t190 + t104 * t197 / 0.2e1 - t103 * t198 / 0.2e1 - t132 * t204 + t128 * (t158 * Ifges(3,1) + t211) / 0.2e1 - t25 * t217 + Ifges(2,3) * qJDD(1); Ifges(3,6) * t127 + Ifges(3,5) * t128 - t115 * mrSges(3,2) - t116 * mrSges(3,1) + t101 * t8 + t102 * t9 - t91 * t38 - t82 * t92 - t81 * t93 - m(4) * (t75 * t81 + t76 * t82) + t165 + t76 * t210 + t243 * t59 + t244 * t58 + (t157 * t45 + t161 * t44 + (-t157 * t93 + t161 * t92) * qJD(3) + (t157 * t28 + t161 * t29 + t195 * t76 - t196 * t75 - t218) * m(4)) * pkin(2) + (-qJD(2) * t177 / 0.2e1 + t103 * t225 + (t205 * t225 - t174 / 0.2e1 + t175) * qJD(1) + (t162 * t130 + t158 * t131) * pkin(5) + (m(4) * t133 - t74) * t222 - (t145 + t104) * t162 / 0.2e1) * qJD(1) + Ifges(3,3) * qJDD(2) + t242 * g(3) + t234 * (-m(5) * (-pkin(3) * t149 - t222) + t215 + (m(4) * pkin(2) + mrSges(3,1)) * t158 + t167) + (-t218 * pkin(2) - t139 * g(3) + t101 * t3 + t102 * t2 + t243 * t19 + t244 * t20 - t88 * t91) * m(5); -t75 * t92 - t22 * t58 - t21 * t59 + t165 + t170 * g(3) - m(5) * (t19 * t21 + t20 * t22) + (-t106 * t38 + t156 * t9 + t160 * t8 + (-t156 * t59 + t160 * t58) * qJD(4) + (-g(3) * t150 - t106 * t88 + t234 * t149 + t156 * t2 + t160 * t3 - t19 * t194 + t20 * t193) * m(5)) * pkin(3) + (t93 + t210) * t76 + t234 * t167; -t88 * (mrSges(5,1) * t69 + t239) + (-t223 + t238) * t230 + t33 * t229 + (-Ifges(5,6) * t69 + t237) * t226 - t19 * t58 + t20 * t59 - g(3) * t181 + (t180 * t19 + t20 * t69) * mrSges(5,3) + t234 * t179 + t241;];
tau = t1;
