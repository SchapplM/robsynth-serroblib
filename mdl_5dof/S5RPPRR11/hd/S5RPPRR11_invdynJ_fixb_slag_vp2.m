% Calculate vector of inverse dynamics joint torques for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:32
% EndTime: 2019-12-31 18:05:41
% DurationCPUTime: 4.26s
% Computational Cost: add. (1584->343), mult. (2932->478), div. (0->0), fcn. (1422->6), ass. (0->161)
t80 = sin(qJ(4));
t210 = -t80 / 0.2e1;
t83 = cos(qJ(4));
t141 = t83 * qJD(1);
t82 = cos(qJ(5));
t148 = qJD(4) * t82;
t79 = sin(qJ(5));
t49 = -t141 * t79 + t148;
t209 = -t49 / 0.2e1;
t150 = qJD(4) * t79;
t50 = t141 * t82 + t150;
t208 = -t50 / 0.2e1;
t142 = t80 * qJD(1);
t65 = qJD(5) + t142;
t207 = -t65 / 0.2e1;
t176 = Ifges(5,4) * t83;
t206 = Ifges(5,2) * t210 + t176 / 0.2e1;
t154 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t49 + mrSges(6,2) * t50 + mrSges(5,3) * t141;
t139 = qJD(1) * qJD(4);
t55 = qJDD(1) * t83 - t139 * t80;
t18 = qJD(5) * t49 + qJDD(4) * t79 + t55 * t82;
t188 = t18 / 0.2e1;
t19 = -qJD(5) * t50 + qJDD(4) * t82 - t55 * t79;
t187 = t19 / 0.2e1;
t56 = -qJDD(1) * t80 - t139 * t83;
t46 = qJDD(5) - t56;
t186 = t46 / 0.2e1;
t183 = -m(5) - m(4);
t205 = -m(5) - m(6);
t5 = -mrSges(6,1) * t19 + mrSges(6,2) * t18;
t204 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t55 - t5;
t203 = mrSges(3,2) - mrSges(4,3);
t202 = -mrSges(4,2) - mrSges(3,3);
t73 = qJD(1) * qJD(2);
t201 = qJDD(1) * qJ(2) + t73;
t177 = Ifges(5,4) * t80;
t112 = t83 * Ifges(5,1) - t177;
t45 = Ifges(6,4) * t49;
t15 = t50 * Ifges(6,1) + t65 * Ifges(6,5) + t45;
t200 = Ifges(5,5) * qJD(4) + qJD(1) * t112 + t82 * t15;
t121 = pkin(4) * t80 - pkin(7) * t83;
t78 = pkin(1) + qJ(3);
t57 = t121 + t78;
t31 = qJD(1) * t57 - qJD(2);
t66 = qJD(1) * qJ(2) + qJD(3);
t62 = -qJD(1) * pkin(6) + t66;
t170 = t62 * t80;
t37 = qJD(4) * pkin(7) + t170;
t10 = t31 * t82 - t37 * t79;
t11 = t31 * t79 + t37 * t82;
t145 = qJD(5) * t82;
t146 = qJD(5) * t79;
t199 = -t10 * t145 - t11 * t146;
t147 = qJD(4) * t83;
t77 = qJ(2) - pkin(6);
t198 = qJD(2) * t80 + t77 * t147;
t149 = qJD(4) * t80;
t61 = qJDD(3) + t201;
t51 = -qJDD(1) * pkin(6) + t61;
t25 = -t149 * t62 + t51 * t83;
t26 = t62 * t147 + t80 * t51;
t196 = t25 * t83 + t26 * t80;
t8 = mrSges(6,1) * t46 - mrSges(6,3) * t18;
t9 = -mrSges(6,2) * t46 + mrSges(6,3) * t19;
t195 = -t79 * t8 + t82 * t9;
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t194 = -g(1) * t84 - g(2) * t81;
t193 = Ifges(5,6) * qJD(4) / 0.2e1 + qJD(1) * t206 + Ifges(6,5) * t208 + Ifges(6,6) * t209 + Ifges(6,3) * t207;
t192 = mrSges(2,2) + mrSges(5,3) + t202;
t113 = mrSges(6,1) * t79 + mrSges(6,2) * t82;
t175 = Ifges(6,4) * t50;
t14 = Ifges(6,2) * t49 + Ifges(6,6) * t65 + t175;
t169 = t62 * t83;
t38 = -qJD(4) * pkin(4) - t169;
t191 = t113 * t38 - t79 * t14 / 0.2e1;
t115 = t80 * mrSges(5,1) + t83 * mrSges(5,2);
t157 = t83 * mrSges(6,3);
t190 = mrSges(2,1) + t115 - t157 - t203;
t85 = qJD(1) ^ 2;
t189 = Ifges(6,1) * t188 + Ifges(6,4) * t187 + Ifges(6,5) * t186;
t184 = t50 / 0.2e1;
t140 = qJD(1) * qJD(3);
t52 = qJDD(1) * t78 - qJDD(2) + t140;
t12 = -pkin(4) * t56 - pkin(7) * t55 + t52;
t21 = qJDD(4) * pkin(7) + t26;
t2 = -qJD(5) * t11 + t12 * t82 - t21 * t79;
t180 = t2 * t79;
t174 = Ifges(6,4) * t79;
t173 = Ifges(6,4) * t82;
t168 = t77 * t80;
t166 = t79 * t80;
t165 = t79 * t81;
t164 = t79 * t83;
t163 = t79 * t84;
t162 = t81 * t82;
t161 = t82 * mrSges(6,3);
t159 = t82 * t83;
t158 = t82 * t84;
t153 = t84 * pkin(1) + t81 * qJ(2);
t144 = qJD(5) * t83;
t143 = qJDD(1) * pkin(1);
t138 = Ifges(6,5) * t18 + Ifges(6,6) * t19 + Ifges(6,3) * t46;
t137 = t84 * qJ(3) + t153;
t134 = t79 * t149;
t129 = m(6) * pkin(7) + mrSges(6,3);
t127 = (t80 ^ 2 + t83 ^ 2) * t62;
t123 = -t139 / 0.2e1;
t122 = pkin(4) * t83 + pkin(7) * t80;
t119 = qJD(4) * t122 - qJD(5) * t168 + qJD(3);
t1 = qJD(5) * t10 + t12 * t79 + t21 * t82;
t118 = t1 * t82 - t180;
t117 = -t56 * mrSges(5,1) + t55 * mrSges(5,2);
t116 = mrSges(5,1) * t83 - mrSges(5,2) * t80;
t114 = -mrSges(6,1) * t82 + mrSges(6,2) * t79;
t111 = Ifges(6,1) * t82 - t174;
t110 = Ifges(6,1) * t79 + t173;
t108 = -Ifges(6,2) * t79 + t173;
t107 = Ifges(6,2) * t82 + t174;
t106 = -Ifges(5,5) * t80 - Ifges(5,6) * t83;
t105 = Ifges(6,5) * t82 - Ifges(6,6) * t79;
t104 = Ifges(6,5) * t79 + Ifges(6,6) * t82;
t103 = t10 * t82 + t11 * t79;
t102 = t10 * t79 - t11 * t82;
t29 = -mrSges(6,2) * t65 + mrSges(6,3) * t49;
t30 = mrSges(6,1) * t65 - mrSges(6,3) * t50;
t101 = -t82 * t29 + t79 * t30;
t100 = -t79 * t29 - t82 * t30;
t63 = qJD(1) * t78 - qJD(2);
t98 = qJD(3) * t63 + t52 * t78;
t58 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t142;
t97 = -t101 + t58;
t96 = t63 * t116;
t95 = t80 * (-Ifges(5,2) * t83 - t177);
t94 = t83 * (-Ifges(5,1) * t80 - t176);
t93 = m(6) * pkin(4) - t114;
t92 = t144 * t79 + t148 * t80;
t91 = -t144 * t82 + t134;
t90 = qJD(5) * t57 + t198;
t89 = Ifges(6,5) * t83 - t111 * t80;
t88 = Ifges(6,6) * t83 - t108 * t80;
t87 = Ifges(6,3) * t83 - t105 * t80;
t70 = t84 * qJ(2);
t67 = qJDD(2) - t143;
t54 = t122 * qJD(1);
t53 = t115 * qJD(1);
t44 = t113 * t83;
t42 = t158 * t80 - t165;
t41 = -t163 * t80 - t162;
t40 = -t162 * t80 - t163;
t39 = t165 * t80 - t158;
t33 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t56;
t28 = t168 * t82 + t57 * t79;
t27 = -t166 * t77 + t57 * t82;
t23 = t159 * t62 + t54 * t79;
t22 = -t164 * t62 + t54 * t82;
t20 = -qJDD(4) * pkin(4) - t25;
t7 = t119 * t82 - t79 * t90;
t6 = t119 * t79 + t82 * t90;
t3 = t18 * Ifges(6,4) + t19 * Ifges(6,2) + t46 * Ifges(6,6);
t4 = [(t67 - t143) * mrSges(3,2) + (mrSges(4,3) * t78 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + t83 * (t55 * Ifges(5,1) + t56 * Ifges(5,4)) / 0.2e1 + (t10 * t92 + t11 * t91) * mrSges(6,3) + (Ifges(5,4) * t55 + Ifges(5,2) * t56) * t210 - (t82 * t14 + t79 * t15) * t144 / 0.2e1 + m(6) * (t1 * t28 + t10 * t7 + t11 * t6 + t2 * t27) + qJDD(4) * (Ifges(5,5) * t83 - Ifges(5,6) * t80) + qJD(3) * t53 + t20 * t44 + t78 * t117 + (m(6) * t38 + t154) * (-qJD(2) * t83 + t77 * t149) + t38 * (-mrSges(6,1) * t91 - mrSges(6,2) * t92) + t27 * t8 + t28 * t9 + t6 * t29 + t7 * t30 + t52 * t115 + m(3) * (-pkin(1) * t67 + (t201 + t73) * qJ(2)) + (-m(6) * t20 + t204) * t77 * t83 + m(4) * (qJ(2) * t61 + qJD(2) * t66 + t98) + qJD(4) ^ 2 * t106 / 0.2e1 + t55 * t112 / 0.2e1 + t14 * t134 / 0.2e1 + t80 * t138 / 0.2e1 + t94 * t139 / 0.2e1 + t95 * t123 + t49 * (qJD(4) * t88 - t107 * t144) / 0.2e1 + t65 * (qJD(4) * t87 - t104 * t144) / 0.2e1 + qJD(4) * t96 + t33 * t168 + (qJD(4) * t89 - t110 * t144) * t184 + (Ifges(6,3) * t80 + t105 * t83) * t186 + (Ifges(6,6) * t80 + t108 * t83) * t187 + (Ifges(6,5) * t80 + t111 * t83) * t188 + t159 * t189 + t56 * t206 + (t10 * mrSges(6,1) - t11 * mrSges(6,2) - t193) * t147 + m(5) * (qJD(2) * t127 + t196 * t77 + t98) - t196 * mrSges(5,3) + t198 * t58 + t2 * (mrSges(6,1) * t80 - t157 * t82) + t1 * (-mrSges(6,2) * t80 - t157 * t79) - t200 * t149 / 0.2e1 + (t61 + t201) * mrSges(4,2) - t3 * t164 / 0.2e1 + (-m(3) * t153 - m(4) * t137 - t42 * mrSges(6,1) - t41 * mrSges(6,2) + t205 * (-pkin(6) * t81 + t137) + t192 * t81 + (-m(6) * t121 - t190) * t84) * g(2) + (-t40 * mrSges(6,1) - t39 * mrSges(6,2) + (-m(3) - m(4)) * t70 + t205 * (-pkin(6) * t84 + t70) + t192 * t84 + (m(3) * pkin(1) + m(6) * t57 - t183 * t78 + t190) * t81) * g(1) + 0.2e1 * t201 * mrSges(3,3) + (t52 + t140) * mrSges(4,3); -t79 * t9 - t82 * t8 + t203 * qJDD(1) + t101 * qJD(5) + (-m(3) * qJ(2) + t202) * t85 + m(6) * (qJD(5) * t102 - t1 * t79 - t2 * t82) + m(3) * t67 - m(5) * t52 + (t154 * t83 - t97 * t80 - m(6) * (-t102 * t80 - t38 * t83) - m(5) * t127) * qJD(1) - t117 + (-g(1) * t81 + g(2) * t84) * (m(3) + m(6) - t183) + (-qJD(1) * t66 - t52) * m(4); qJDD(1) * mrSges(4,2) - t85 * mrSges(4,3) + m(4) * t61 + (t97 * qJD(4) + m(6) * (-t10 * t150 + t11 * t148 - t20) + m(5) * t25 + t204) * t83 + (t33 + t100 * qJD(5) + t154 * qJD(4) + m(6) * (qJD(4) * t38 + t118 + t199) + m(5) * t26 + t195) * t80 + (-m(6) * t103 + t183 * t63 + t100 - t53) * qJD(1) + (m(4) - t205) * t194; (-t94 / 0.2e1 + t95 / 0.2e1) * t85 + (-t10 * (mrSges(6,1) * t83 + t161 * t80) - t11 * (-mrSges(6,2) * t83 + mrSges(6,3) * t166) - t96) * qJD(1) + (t191 + t200 / 0.2e1) * t142 + (-t20 * pkin(4) - t10 * t22 - t11 * t23 - t170 * t38) * m(6) + (t105 * t65 + t108 * t49 + t111 * t50) * qJD(5) / 0.2e1 - (t49 * t88 + t50 * t89 + t65 * t87) * qJD(1) / 0.2e1 + t82 * t3 / 0.2e1 + Ifges(5,5) * t55 + Ifges(5,6) * t56 - t23 * t29 - t22 * t30 + t25 * mrSges(5,1) - t26 * mrSges(5,2) - pkin(4) * t5 + (-t129 * t83 + t80 * t93 + t115) * g(3) + t20 * t114 + (t129 * t80 + t83 * t93 + t116) * t194 + Ifges(5,3) * qJDD(4) + t106 * t123 + t15 * t145 / 0.2e1 + t191 * qJD(5) + t1 * t161 + t104 * t186 + t107 * t187 + t110 * t188 + t79 * t189 + t193 * t141 + (m(6) * (-qJD(5) * t103 + t118) - t30 * t145 - t29 * t146 + t195) * pkin(7) + (-t180 + t199) * mrSges(6,3) - t154 * t170 - t58 * t169; -t1 * mrSges(6,2) + t2 * mrSges(6,1) - t38 * (mrSges(6,1) * t50 + mrSges(6,2) * t49) + (Ifges(6,1) * t49 - t175) * t208 + t14 * t184 + (Ifges(6,5) * t49 - Ifges(6,6) * t50) * t207 - t10 * t29 + t11 * t30 - g(1) * (mrSges(6,1) * t41 - mrSges(6,2) * t42) - g(2) * (-mrSges(6,1) * t39 + mrSges(6,2) * t40) + g(3) * t44 + (t10 * t49 + t11 * t50) * mrSges(6,3) + t138 + (-Ifges(6,2) * t50 + t15 + t45) * t209;];
tau = t4;
