% Calculate vector of inverse dynamics joint torques for
% S4RRPR10
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:22
% DurationCPUTime: 6.31s
% Computational Cost: add. (1402->370), mult. (3077->522), div. (0->0), fcn. (1600->6), ass. (0->184)
t99 = sin(qJ(2));
t176 = t99 * qJD(1);
t181 = qJDD(1) * pkin(1);
t102 = cos(qJ(2));
t168 = qJD(1) * qJD(2);
t65 = qJDD(1) * t99 + t102 * t168;
t110 = -qJ(3) * t65 - qJD(3) * t176 - t181;
t218 = pkin(2) + pkin(6);
t157 = t99 * t168;
t167 = t102 * qJDD(1);
t64 = t157 - t167;
t10 = t218 * t64 + t110;
t101 = cos(qJ(4));
t92 = t99 * qJ(3);
t154 = -pkin(1) - t92;
t31 = (-t102 * t218 + t154) * qJD(1);
t86 = pkin(5) * t176;
t60 = -pkin(3) * t176 - t86;
t231 = -t60 + qJD(3);
t33 = -qJD(2) * t218 + t231;
t98 = sin(qJ(4));
t11 = t101 * t33 - t31 * t98;
t54 = t65 * pkin(5);
t158 = qJDD(3) + t54;
t23 = pkin(3) * t65 - qJDD(2) * t218 + t158;
t1 = qJD(4) * t11 + t10 * t101 + t23 * t98;
t249 = t1 * mrSges(5,2);
t12 = t101 * t31 + t33 * t98;
t2 = -qJD(4) * t12 - t10 * t98 + t101 * t23;
t248 = t2 * mrSges(5,1);
t244 = -Ifges(4,4) + Ifges(3,5);
t243 = Ifges(4,5) - Ifges(3,6);
t100 = sin(qJ(1));
t208 = g(2) * t100;
t136 = t102 * mrSges(4,2) - t99 * mrSges(4,3);
t139 = mrSges(3,1) * t102 - mrSges(3,2) * t99;
t247 = t136 - t139;
t103 = cos(qJ(1));
t232 = g(1) * t103 + t208;
t169 = t102 * qJD(1);
t56 = -qJD(2) * t98 - t101 * t169;
t81 = qJD(4) + t176;
t28 = -mrSges(5,2) * t81 + mrSges(5,3) * t56;
t175 = qJD(2) * t101;
t114 = t169 * t98 - t175;
t29 = mrSges(5,1) * t81 + mrSges(5,3) * t114;
t122 = t101 * t28 - t98 * t29;
t21 = qJD(4) * t56 + qJDD(2) * t101 + t64 * t98;
t55 = qJDD(4) + t65;
t8 = mrSges(5,1) * t55 - mrSges(5,3) * t21;
t22 = qJD(4) * t114 - qJDD(2) * t98 + t101 * t64;
t9 = -mrSges(5,2) * t55 + mrSges(5,3) * t22;
t246 = t122 * qJD(4) + t101 * t8 + t98 * t9;
t245 = m(4) + m(5);
t85 = Ifges(3,4) * t169;
t242 = Ifges(3,1) * t176 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t114 + t56 * Ifges(5,6) + t81 * Ifges(5,3) + t85;
t27 = -mrSges(5,1) * t56 - mrSges(5,2) * t114;
t160 = mrSges(4,1) * t169;
t71 = -qJD(2) * mrSges(4,3) - t160;
t241 = -t71 + t27;
t240 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t169 + t71;
t162 = mrSges(4,1) * t176;
t239 = mrSges(3,3) * t176 + t162 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t195 = Ifges(4,6) * t102;
t211 = Ifges(4,6) * t99;
t237 = t99 * (-Ifges(4,2) * t102 + t211) + t102 * (Ifges(4,3) * t99 - t195);
t236 = t102 * t244 + t243 * t99;
t84 = pkin(5) * t167;
t53 = -pkin(5) * t157 + t84;
t235 = t102 * t53 + t54 * t99;
t32 = qJD(2) * (-qJD(3) + t86) - qJDD(2) * qJ(3) - t84;
t34 = -qJDD(2) * pkin(2) + t158;
t234 = -t102 * t32 + t34 * t99;
t230 = -mrSges(2,1) + t247;
t125 = -t102 * Ifges(4,3) - t211;
t213 = Ifges(5,4) * t114;
t14 = Ifges(5,2) * t56 + Ifges(5,6) * t81 - t213;
t52 = Ifges(5,4) * t56;
t15 = -Ifges(5,1) * t114 + Ifges(5,5) * t81 + t52;
t229 = Ifges(4,5) * qJD(2) + qJD(1) * t125 + t101 * t14 + t98 * t15;
t228 = -m(5) * pkin(3) + mrSges(2,2) - mrSges(3,3);
t177 = t102 * t103;
t184 = qJ(3) * t102;
t227 = -g(1) * qJ(3) * t177 - t184 * t208;
t225 = -t21 * Ifges(5,4) / 0.2e1 - t22 * Ifges(5,2) / 0.2e1 - t55 * Ifges(5,6) / 0.2e1;
t224 = t21 / 0.2e1;
t223 = t22 / 0.2e1;
t222 = t55 / 0.2e1;
t220 = -t114 / 0.2e1;
t217 = pkin(3) + pkin(5);
t216 = pkin(5) * t99;
t214 = Ifges(3,4) * t99;
t212 = Ifges(5,4) * t98;
t210 = pkin(5) * t102;
t94 = t102 * pkin(2);
t204 = t98 * mrSges(5,3);
t199 = t94 + t92;
t198 = t103 * pkin(1) + t100 * pkin(5);
t197 = Ifges(3,4) * t102;
t196 = Ifges(5,4) * t101;
t194 = t100 * t98;
t193 = t101 * mrSges(5,3);
t190 = t102 * mrSges(4,3);
t187 = t102 * t98;
t186 = t103 * t98;
t183 = qJD(2) * t99;
t182 = qJD(4) * t98;
t180 = t100 * t101;
t179 = t101 * t102;
t178 = t101 * t103;
t174 = qJD(2) * t102;
t173 = qJD(4) * t101;
t172 = qJD(4) * t102;
t164 = Ifges(5,5) * t21 + Ifges(5,6) * t22 + Ifges(5,3) * t55;
t88 = pkin(5) * t169;
t75 = t217 * t102;
t159 = t98 * t172;
t37 = t65 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t151 = -t173 / 0.2e1;
t149 = -t168 / 0.2e1;
t147 = pkin(2) * t177 + t103 * t92 + t198;
t146 = pkin(2) * t183 - qJD(3) * t99;
t69 = -qJD(2) * qJ(3) - t88;
t144 = pkin(6) * t102 + t199;
t143 = -m(5) * t218 - mrSges(5,3);
t140 = mrSges(3,1) * t99 + mrSges(3,2) * t102;
t138 = mrSges(5,1) * t98 + mrSges(5,2) * t101;
t137 = mrSges(5,1) * t101 - mrSges(5,2) * t98;
t135 = Ifges(5,1) * t98 + t196;
t134 = Ifges(5,1) * t101 - t212;
t133 = t102 * Ifges(3,2) + t214;
t131 = Ifges(5,2) * t101 + t212;
t130 = -Ifges(5,2) * t98 + t196;
t128 = Ifges(5,5) * t98 + Ifges(5,6) * t101;
t127 = Ifges(5,5) * t101 - Ifges(5,6) * t98;
t126 = -t99 * Ifges(4,2) - t195;
t124 = pkin(6) * t99 - t184;
t123 = t12 * t101 - t11 * t98;
t51 = -pkin(1) - t144;
t74 = t217 * t99;
t25 = t101 * t74 - t51 * t98;
t26 = t101 * t51 + t74 * t98;
t66 = -qJD(2) * pkin(2) + qJD(3) + t86;
t121 = t66 * t102 + t69 * t99;
t120 = t154 - t94;
t119 = pkin(1) * t140;
t44 = t120 * qJD(1);
t118 = t44 * (-t99 * mrSges(4,2) - t190);
t117 = t99 * (Ifges(3,1) * t102 - t214);
t113 = -t101 * t172 + t183 * t98;
t112 = t175 * t99 + t159;
t109 = Ifges(5,5) * t102 + t135 * t99;
t108 = Ifges(5,6) * t102 + t131 * t99;
t107 = Ifges(5,3) * t102 + t128 * t99;
t106 = qJD(4) * t123 + t1 * t98 + t101 * t2;
t89 = pkin(3) * t169;
t87 = pkin(2) * t176;
t67 = -pkin(1) - t199;
t63 = qJD(2) * t75;
t62 = t217 * t183;
t61 = t88 + t89;
t59 = -qJ(3) * t169 + t87;
t58 = t136 * qJD(1);
t50 = t137 * t102;
t48 = -t194 * t99 + t178;
t47 = t180 * t99 + t186;
t46 = t186 * t99 + t180;
t45 = t178 * t99 - t194;
t43 = -t69 + t89;
t42 = Ifges(4,4) * qJD(2) + qJD(1) * t126;
t39 = Ifges(3,6) * qJD(2) + qJD(1) * t133;
t38 = -qJ(3) * t174 + t146;
t36 = mrSges(4,1) * t64 - qJDD(2) * mrSges(4,3);
t35 = qJD(1) * t124 + t87;
t30 = qJD(2) * t124 + t146;
t24 = -pkin(3) * t64 - t32;
t20 = pkin(2) * t64 + t110;
t19 = t101 * t35 + t61 * t98;
t18 = t101 * t61 - t35 * t98;
t7 = -qJD(4) * t26 + t101 * t63 - t30 * t98;
t6 = qJD(4) * t25 + t101 * t30 + t63 * t98;
t5 = -mrSges(5,1) * t22 + mrSges(5,2) * t21;
t4 = Ifges(5,1) * t21 + Ifges(5,4) * t22 + Ifges(5,5) * t55;
t3 = [(t239 * pkin(5) - t42 / 0.2e1 - t12 * mrSges(5,2) + t11 * mrSges(5,1) + t66 * mrSges(4,1) + t242 / 0.2e1) * t174 + (-m(4) * t147 - m(3) * t198 - m(5) * (pkin(6) * t177 + t147) - t46 * mrSges(5,1) - t45 * mrSges(5,2) + t228 * t100 + t230 * t103) * g(2) + (-g(2) * t177 - t1 * t179 - t11 * t113 + t112 * t12 + t187 * t2) * mrSges(5,3) + (-qJDD(2) * mrSges(3,2) - t36) * t210 + t65 * t197 / 0.2e1 + (t229 / 0.2e1 + t240 * pkin(5) - t39 / 0.2e1 + t69 * mrSges(4,1)) * t183 + t139 * t181 + (-t48 * mrSges(5,1) + t47 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * t120 - m(5) * t154 - t102 * t143 - t230) * t100 + ((-m(3) - t245) * pkin(5) + t228) * t103) * g(1) + t43 * (-mrSges(5,1) * t112 + mrSges(5,2) * t113) + t99 * t248 + t65 * t99 * Ifges(3,1) + (-Ifges(3,4) * t64 + Ifges(3,5) * qJDD(2) + t164) * t99 / 0.2e1 - t102 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t65 + Ifges(4,3) * t64) / 0.2e1 + t102 * (Ifges(3,4) * t65 - Ifges(3,2) * t64 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t99 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t65 + Ifges(4,6) * t64) / 0.2e1 - t4 * t187 / 0.2e1 + t67 * (-mrSges(4,2) * t64 - mrSges(4,3) * t65) + t75 * t5 - t62 * t27 - pkin(1) * (mrSges(3,1) * t64 + mrSges(3,2) * t65) + t38 * t58 + t24 * t50 + t56 * (qJD(2) * t108 - t130 * t172) / 0.2e1 + t81 * (qJD(2) * t107 - t127 * t172) / 0.2e1 + t6 * t28 + t7 * t29 + t25 * t8 + t26 * t9 - t119 * t168 + (-t232 + t234) * mrSges(4,1) + t102 * t15 * t151 + (-qJDD(2) * mrSges(3,1) + t37) * t216 - t99 * t249 + m(5) * (t1 * t26 + t11 * t7 + t12 * t6 + t2 * t25 + t24 * t75 - t43 * t62) + (t102 * (-Ifges(3,2) * t99 + t197) + t117) * t168 / 0.2e1 + (-t210 * t64 + t216 * t65 + t235) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t235) + t236 * qJD(2) ^ 2 / 0.2e1 + t237 * t149 + m(4) * (t20 * t67 + t38 * t44 + (qJD(2) * t121 + t234) * pkin(5)) + t64 * t125 / 0.2e1 - t65 * t126 / 0.2e1 - t64 * t133 / 0.2e1 + t20 * t136 + qJD(2) * t118 + (qJD(2) * t109 - t134 * t172) * t220 + (Ifges(5,3) * t99 - t102 * t128) * t222 + (Ifges(5,6) * t99 - t102 * t131) * t223 + (Ifges(5,5) * t99 - t102 * t135) * t224 + t179 * t225 + (-t102 * t243 + t244 * t99) * qJDD(2) / 0.2e1 + Ifges(2,3) * qJDD(1) + t14 * t159 / 0.2e1; (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (-pkin(2) * t34 - qJ(3) * t32 - qJD(3) * t69 - t44 * t59 + t227) * m(4) + (-t36 + t5) * qJ(3) + (t24 * qJ(3) - t11 * t18 - t12 * t19 + t231 * t43 + t227) * m(5) - t229 * t176 / 0.2e1 + (-t102 * t138 + t140 - t190 + (m(4) * pkin(2) - mrSges(4,2) - t143) * t99) * t232 + (-t117 / 0.2e1 + t119 + t237 / 0.2e1) * qJD(1) ^ 2 - (-t114 * t135 + t128 * t81 + t131 * t56) * qJD(4) / 0.2e1 - (t107 * t81 + t108 * t56 - t109 * t114) * qJD(1) / 0.2e1 + (-t11 * (mrSges(5,1) * t102 - t204 * t99) - t12 * (-mrSges(5,2) * t102 + t193 * t99) - t118 - m(4) * pkin(5) * t121) * qJD(1) + t81 * t43 * t137 + (-m(5) * t106 - t246) * t218 - t1 * t204 + t101 * t4 / 0.2e1 - t2 * t193 - t15 * t182 / 0.2e1 - t59 * t58 - t60 * t27 + t39 * t176 / 0.2e1 - t53 * mrSges(3,2) - t54 * mrSges(3,1) + t34 * mrSges(4,2) - pkin(2) * t37 - t19 * t28 - t18 * t29 + t42 * t169 / 0.2e1 - t32 * mrSges(4,3) - t69 * t162 - t66 * t160 + (-m(4) * t199 - m(5) * t144 - t102 * mrSges(5,3) - t138 * t99 + t247) * g(3) + t236 * t149 + (t11 * t182 - t12 * t173) * mrSges(5,3) + t14 * t151 + t24 * t138 + t127 * t222 + t130 * t223 + t134 * t224 + t98 * t225 - t239 * t88 - t240 * t86 + t241 * qJD(3) - (-Ifges(3,2) * t176 + t242 + t85) * t169 / 0.2e1 + t243 * t64 + t244 * t65; -t241 * qJD(2) + t245 * t102 * g(3) + ((t122 + t58) * qJD(1) - t245 * t232) * t99 + t37 + (-qJD(2) * t43 + t123 * t176 + t106) * m(5) + (qJD(2) * t69 + t176 * t44 + t34) * m(4) + t246; -t249 + t248 - t43 * (-mrSges(5,1) * t114 + mrSges(5,2) * t56) + t114 * (Ifges(5,1) * t56 + t213) / 0.2e1 + t14 * t220 - t81 * (Ifges(5,5) * t56 + Ifges(5,6) * t114) / 0.2e1 - t11 * t28 + t12 * t29 - g(1) * (mrSges(5,1) * t45 - mrSges(5,2) * t46) - g(2) * (mrSges(5,1) * t47 + mrSges(5,2) * t48) + g(3) * t50 + (t11 * t56 - t114 * t12) * mrSges(5,3) + t164 - (Ifges(5,2) * t114 + t15 + t52) * t56 / 0.2e1;];
tau = t3;
