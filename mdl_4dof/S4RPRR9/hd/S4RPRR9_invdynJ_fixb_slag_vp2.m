% Calculate vector of inverse dynamics joint torques for
% S4RPRR9
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:08
% EndTime: 2019-12-31 16:56:15
% DurationCPUTime: 3.57s
% Computational Cost: add. (1208->292), mult. (2441->421), div. (0->0), fcn. (1262->6), ass. (0->147)
t188 = -Ifges(4,2) / 0.2e1;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t105 = pkin(3) * t70 - pkin(6) * t73;
t55 = qJ(2) + t105;
t35 = t55 * qJD(1);
t75 = -pkin(1) - pkin(5);
t60 = qJD(1) * t75 + qJD(2);
t152 = t60 * t70;
t40 = qJD(3) * pkin(6) + t152;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t13 = t35 * t72 - t40 * t69;
t122 = qJD(1) * qJD(3);
t53 = qJDD(1) * t73 - t122 * t70;
t54 = -qJDD(1) * t70 - t122 * t73;
t123 = qJD(1) * qJD(2);
t61 = qJDD(1) * qJ(2) + t123;
t19 = -pkin(3) * t54 - pkin(6) * t53 + t61;
t130 = qJD(3) * t73;
t59 = qJDD(1) * t75 + qJDD(2);
t26 = t60 * t130 + t70 * t59;
t24 = qJDD(3) * pkin(6) + t26;
t1 = qJD(4) * t13 + t19 * t69 + t24 * t72;
t14 = t35 * t69 + t40 * t72;
t2 = -qJD(4) * t14 + t19 * t72 - t24 * t69;
t187 = t2 * mrSges(5,1) - t1 * mrSges(5,2);
t131 = qJD(3) * t72;
t134 = qJD(1) * t73;
t47 = -t134 * t69 + t131;
t186 = -t47 / 0.2e1;
t133 = qJD(3) * t69;
t48 = t134 * t72 + t133;
t185 = -t48 / 0.2e1;
t125 = t70 * qJD(1);
t62 = qJD(4) + t125;
t184 = -t62 / 0.2e1;
t183 = t73 / 0.2e1;
t158 = Ifges(4,4) * t73;
t182 = t70 * t188 + t158 / 0.2e1;
t71 = sin(qJ(1));
t74 = cos(qJ(1));
t104 = -g(1) * t71 + g(2) * t74;
t17 = qJD(4) * t47 + qJDD(3) * t69 + t53 * t72;
t167 = t17 / 0.2e1;
t18 = -qJD(4) * t48 + qJDD(3) * t72 - t53 * t69;
t166 = t18 / 0.2e1;
t45 = qJDD(4) - t54;
t165 = t45 / 0.2e1;
t181 = -m(4) - m(5);
t5 = -mrSges(5,1) * t18 + mrSges(5,2) * t17;
t180 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t53 - t5;
t179 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t47 - mrSges(5,2) * t48 - mrSges(4,3) * t134;
t44 = Ifges(5,4) * t47;
t12 = Ifges(5,1) * t48 + Ifges(5,5) * t62 + t44;
t146 = t70 * Ifges(4,4);
t98 = t73 * Ifges(4,1) - t146;
t178 = Ifges(4,5) * qJD(3) + qJD(1) * t98 + t72 * t12;
t127 = qJD(4) * t72;
t128 = qJD(4) * t69;
t177 = -t13 * t127 - t14 * t128;
t132 = qJD(3) * t70;
t25 = -t132 * t60 + t59 * t73;
t176 = t25 * t73 + t26 * t70;
t6 = mrSges(5,1) * t45 - mrSges(5,3) * t17;
t7 = -mrSges(5,2) * t45 + mrSges(5,3) * t18;
t175 = -t69 * t6 + t72 * t7;
t102 = mrSges(4,1) * t73 - mrSges(4,2) * t70;
t174 = qJ(2) * t102 + (-Ifges(4,1) * t70 - t158) * t183;
t173 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t172 = Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t182 + Ifges(5,5) * t185 + Ifges(5,6) * t186 + Ifges(5,3) * t184;
t157 = Ifges(5,4) * t48;
t11 = Ifges(5,2) * t47 + Ifges(5,6) * t62 + t157;
t151 = t60 * t73;
t41 = -qJD(3) * pkin(3) - t151;
t99 = mrSges(5,1) * t69 + mrSges(5,2) * t72;
t171 = -t69 * t11 / 0.2e1 + t41 * t99;
t101 = t70 * mrSges(4,1) + mrSges(4,2) * t73;
t139 = t73 * mrSges(5,3);
t170 = -m(5) * t105 + mrSges(2,2) - mrSges(3,3) - t101 + t139;
t76 = qJD(1) ^ 2;
t168 = Ifges(5,1) * t167 + Ifges(5,4) * t166 + Ifges(5,5) * t165;
t163 = t48 / 0.2e1;
t156 = Ifges(5,4) * t69;
t155 = Ifges(5,4) * t72;
t150 = t69 * mrSges(5,3);
t149 = t69 * t71;
t148 = t69 * t73;
t147 = t69 * t74;
t145 = t70 * t75;
t144 = t71 * t72;
t143 = t72 * mrSges(5,3);
t141 = t72 * t73;
t140 = t72 * t74;
t135 = t74 * pkin(1) + t71 * qJ(2);
t129 = qJD(3) * t75;
t126 = qJD(4) * t73;
t124 = qJDD(1) * mrSges(3,2);
t120 = Ifges(5,5) * t17 + Ifges(5,6) * t18 + Ifges(5,3) * t45;
t117 = t73 * t129;
t114 = t69 * t132;
t112 = m(5) * pkin(6) + mrSges(5,3);
t110 = -t126 / 0.2e1;
t108 = (t123 + t61) * qJ(2);
t107 = -t122 / 0.2e1;
t106 = pkin(3) * t73 + pkin(6) * t70;
t103 = t1 * t72 - t2 * t69;
t100 = -t72 * mrSges(5,1) + t69 * mrSges(5,2);
t97 = Ifges(5,1) * t72 - t156;
t96 = Ifges(5,1) * t69 + t155;
t94 = -Ifges(5,2) * t69 + t155;
t93 = Ifges(5,2) * t72 + t156;
t92 = -Ifges(4,5) * t70 - Ifges(4,6) * t73;
t91 = Ifges(5,5) * t72 - Ifges(5,6) * t69;
t90 = Ifges(5,5) * t69 + Ifges(5,6) * t72;
t89 = t13 * t72 + t14 * t69;
t27 = -mrSges(5,2) * t62 + mrSges(5,3) * t47;
t28 = mrSges(5,1) * t62 - mrSges(5,3) * t48;
t88 = -t69 * t27 - t72 * t28;
t30 = t145 * t72 + t55 * t69;
t29 = -t145 * t69 + t55 * t72;
t87 = t70 * (-Ifges(4,2) * t73 - t146);
t84 = m(5) * pkin(3) - t100;
t83 = -t76 * qJ(2) + t104;
t82 = t126 * t69 + t131 * t70;
t81 = -t126 * t72 + t114;
t79 = Ifges(5,5) * t73 - t70 * t97;
t78 = Ifges(5,6) * t73 - t70 * t94;
t77 = Ifges(5,3) * t73 - t70 * t91;
t63 = -qJDD(1) * pkin(1) + qJDD(2);
t56 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t125;
t51 = t106 * qJD(1);
t50 = t101 * qJD(1);
t46 = qJD(3) * t106 + qJD(2);
t43 = t99 * t73;
t39 = t140 * t70 - t149;
t38 = t147 * t70 + t144;
t37 = t144 * t70 + t147;
t36 = -t149 * t70 + t140;
t32 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t54;
t23 = -qJDD(3) * pkin(3) - t25;
t21 = t141 * t60 + t51 * t69;
t20 = -t148 * t60 + t51 * t72;
t9 = -qJD(4) * t30 - t117 * t69 + t46 * t72;
t8 = qJD(4) * t29 + t117 * t72 + t46 * t69;
t3 = Ifges(5,4) * t17 + Ifges(5,2) * t18 + Ifges(5,6) * t45;
t4 = [(-Ifges(4,6) * qJDD(3) - t179 * t129 - Ifges(4,4) * t53 / 0.2e1 + t54 * t188 + Ifges(5,3) * t165 + Ifges(5,6) * t166 + Ifges(5,5) * t167 + t120 / 0.2e1 + t187) * t70 + (t72 * t110 + t114 / 0.2e1) * t11 + (-t1 * t69 - t2 * t72) * t139 + (t101 + 0.2e1 * mrSges(3,3)) * t61 + (-m(3) * t135 - t37 * mrSges(5,1) - t36 * mrSges(5,2) + t181 * (t74 * pkin(5) + t135) + t173 * t74 + t170 * t71) * g(2) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + m(4) * (t176 * t75 + t108) - t176 * mrSges(4,3) - t178 * t132 / 0.2e1 + (t13 * mrSges(5,1) - t14 * mrSges(5,2) - t172) * t130 + t174 * t122 + (Ifges(4,5) * qJDD(3) + t91 * t165 + t94 * t166 + t97 * t167 + (-m(5) * t23 + t180) * t75) * t73 + t63 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t54 + mrSges(4,2) * t53) + qJD(2) * t50 + t23 * t43 + t8 * t27 + t9 * t28 + t29 * t6 + t30 * t7 + (t13 * t82 + t14 * t81) * mrSges(5,3) + t32 * t145 + (qJD(3) * t79 - t126 * t96) * t163 + (Ifges(4,1) * t53 + Ifges(4,4) * t54) * t183 + t54 * t182 - t3 * t148 / 0.2e1 + t56 * t117 + t141 * t168 + t69 * t12 * t110 + (-t39 * mrSges(5,1) + t38 * mrSges(5,2) + (m(3) * pkin(1) + t181 * t75 - t173) * t71 + ((-m(3) + t181) * qJ(2) + t170) * t74) * g(1) + m(5) * (t132 * t41 * t75 + t1 * t30 + t13 * t9 + t14 * t8 + t2 * t29) + t41 * (-mrSges(5,1) * t81 - mrSges(5,2) * t82) + t87 * t107 + qJD(3) ^ 2 * t92 / 0.2e1 + t53 * t98 / 0.2e1 + m(3) * (-pkin(1) * t63 + t108) - pkin(1) * t124 + t47 * (qJD(3) * t78 - t126 * t93) / 0.2e1 + t62 * (qJD(3) * t77 - t126 * t90) / 0.2e1; t124 - t76 * mrSges(3,3) + t104 * m(5) + t83 * m(4) + (t63 + t83) * m(3) + (-m(5) * t89 - t50 + t88) * qJD(1) + ((t72 * t27 - t69 * t28 + t56) * qJD(3) + m(4) * t25 + m(5) * (-t13 * t133 + t131 * t14 - t23) + t180) * t73 + (t32 + t88 * qJD(4) - t179 * qJD(3) + m(4) * t26 + m(5) * (qJD(3) * t41 + t103 + t177) + t175) * t70; (-pkin(3) * t23 - t13 * t20 - t14 * t21 - t41 * t152) * m(5) + (t47 * t94 + t48 * t97 + t62 * t91) * qJD(4) / 0.2e1 - (t47 * t78 + t48 * t79 + t62 * t77) * qJD(1) / 0.2e1 + (-t112 * t73 + t70 * t84 + t101) * g(3) + (-t13 * (mrSges(5,1) * t73 + t143 * t70) - t14 * (-mrSges(5,2) * t73 + t150 * t70)) * qJD(1) + t179 * t152 + (-t28 * t127 - t27 * t128 + m(5) * (-qJD(4) * t89 + t103) + t175) * pkin(6) + t177 * mrSges(5,3) + t171 * qJD(4) + t172 * t134 + (t87 / 0.2e1 - t174) * t76 + t72 * t3 / 0.2e1 + Ifges(4,5) * t53 + Ifges(4,6) * t54 - t21 * t27 - t20 * t28 + t1 * t143 + t25 * mrSges(4,1) - t26 * mrSges(4,2) - t2 * t150 - pkin(3) * t5 - t56 * t151 + t90 * t165 + t93 * t166 + t96 * t167 + t69 * t168 + t92 * t107 + t23 * t100 + (t112 * t70 + t73 * t84 + t102) * t104 + (t178 / 0.2e1 + t171) * t125 + Ifges(4,3) * qJDD(3) + t12 * t127 / 0.2e1; -t41 * (mrSges(5,1) * t48 + mrSges(5,2) * t47) + (Ifges(5,1) * t47 - t157) * t185 + t11 * t163 + (Ifges(5,5) * t47 - Ifges(5,6) * t48) * t184 - t13 * t27 + t14 * t28 - g(1) * (mrSges(5,1) * t36 - mrSges(5,2) * t37) - g(2) * (mrSges(5,1) * t38 + mrSges(5,2) * t39) + g(3) * t43 + (t13 * t47 + t14 * t48) * mrSges(5,3) + t120 + (-Ifges(5,2) * t48 + t12 + t44) * t186 + t187;];
tau = t4;
