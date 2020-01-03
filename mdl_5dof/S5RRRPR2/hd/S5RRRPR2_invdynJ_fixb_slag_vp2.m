% Calculate vector of inverse dynamics joint torques for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:26
% DurationCPUTime: 1.89s
% Computational Cost: add. (3151->293), mult. (5136->389), div. (0->0), fcn. (2732->16), ass. (0->146)
t132 = sin(qJ(5));
t200 = t132 / 0.2e1;
t204 = mrSges(5,2) - mrSges(6,3);
t130 = sin(pkin(9));
t131 = cos(pkin(9));
t137 = cos(qJ(3));
t133 = sin(qJ(3));
t177 = t130 * t133;
t186 = pkin(2) * qJD(3);
t138 = cos(qJ(2));
t134 = sin(qJ(2));
t174 = t134 * t137;
t152 = -t133 * t138 - t174;
t185 = qJD(1) * pkin(1);
t70 = t152 * t185;
t175 = t133 * t134;
t151 = t137 * t138 - t175;
t71 = t151 * t185;
t192 = t130 * t70 + t131 * t71 - (t131 * t137 - t177) * t186;
t136 = cos(qJ(5));
t191 = mrSges(6,1) * t136;
t203 = -mrSges(5,1) - t191;
t157 = mrSges(6,1) * t132 + mrSges(6,2) * t136;
t128 = qJD(1) + qJD(2);
t121 = qJD(3) + t128;
t165 = t134 * t185;
t164 = t138 * t185;
t84 = pkin(2) * t128 + t164;
t50 = t133 * t84 + t137 * t165;
t183 = t130 * t50;
t49 = -t133 * t165 + t137 * t84;
t42 = pkin(3) * t121 + t49;
t23 = t131 * t42 - t183;
t18 = -pkin(4) * t121 - t23;
t202 = t18 * t157 + qJD(5) * (Ifges(6,5) * t136 - Ifges(6,6) * t132) / 0.2e1;
t45 = t131 * t50;
t24 = t130 * t42 + t45;
t19 = pkin(8) * t121 + t24;
t14 = qJD(4) * t132 + t136 * t19;
t181 = qJD(5) * t14;
t127 = qJDD(1) + qJDD(2);
t120 = qJDD(3) + t127;
t197 = pkin(1) * t138;
t77 = -qJD(2) * t165 + qJDD(1) * t197;
t62 = pkin(2) * t127 + t77;
t171 = qJD(2) * t138;
t78 = (qJD(1) * t171 + qJDD(1) * t134) * pkin(1);
t22 = -qJD(3) * t50 - t133 * t78 + t137 * t62;
t16 = pkin(3) * t120 + t22;
t21 = t49 * qJD(3) + t133 * t62 + t137 * t78;
t9 = t130 * t16 + t131 * t21;
t6 = pkin(8) * t120 + t9;
t3 = qJDD(4) * t136 - t132 * t6 - t181;
t201 = -t3 - t181;
t176 = t131 * t133;
t193 = t130 * t71 - t131 * t70 - (t130 * t137 + t176) * t186;
t167 = qJD(5) * t136;
t61 = t120 * t132 + t121 * t167;
t48 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t61;
t178 = t121 * t136;
t80 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t178;
t162 = -qJD(5) * t80 - t48;
t13 = qJD(4) * t136 - t132 * t19;
t182 = t136 * t14;
t153 = -t13 * t132 + t182;
t198 = mrSges(6,3) * t13;
t129 = qJ(1) + qJ(2);
t122 = sin(t129);
t111 = pkin(2) * t122;
t123 = cos(t129);
t112 = pkin(2) * t123;
t124 = qJ(3) + t129;
t114 = sin(t124);
t101 = pkin(3) * t114;
t115 = cos(t124);
t102 = pkin(3) * t115;
t196 = pkin(3) * t130;
t195 = pkin(3) * t131;
t2 = qJD(5) * t13 + qJDD(4) * t132 + t136 * t6;
t194 = t136 * t2;
t117 = pkin(2) + t197;
t72 = -pkin(1) * t175 + t137 * t117;
t67 = pkin(3) + t72;
t73 = pkin(1) * t174 + t117 * t133;
t35 = t130 * t67 + t131 * t73;
t190 = mrSges(6,2) * t132;
t189 = Ifges(6,4) * t132;
t188 = Ifges(6,4) * t136;
t187 = Ifges(6,2) * t136;
t116 = pkin(2) * t137 + pkin(3);
t69 = pkin(2) * t176 + t130 * t116;
t179 = t121 * t132;
t135 = sin(qJ(1));
t125 = t135 * pkin(1);
t173 = t111 + t125;
t139 = cos(qJ(1));
t126 = t139 * pkin(1);
t172 = t112 + t126;
t170 = qJD(3) * t133;
t169 = qJD(3) * t137;
t168 = qJD(5) * t132;
t113 = pkin(9) + t124;
t100 = cos(t113);
t99 = sin(t113);
t166 = t100 * pkin(4) + t99 * pkin(8) + t102;
t161 = t112 + t166;
t159 = t99 * pkin(4) - pkin(8) * t100 + t101;
t158 = t201 * mrSges(6,3);
t85 = t190 - t191;
t156 = t187 + t189;
t154 = t13 * t136 + t132 * t14;
t8 = -t130 * t21 + t131 * t16;
t34 = -t130 * t73 + t131 * t67;
t150 = t111 + t159;
t68 = -pkin(2) * t177 + t116 * t131;
t148 = t132 * (Ifges(6,1) * t136 - t189);
t147 = -t114 * mrSges(4,1) - t115 * mrSges(4,2) - t204 * t100 + t203 * t99;
t146 = -t115 * mrSges(4,1) + t114 * mrSges(4,2) + t203 * t100 + t204 * t99;
t145 = -t122 * mrSges(3,1) - t123 * mrSges(3,2) + t147;
t144 = -t123 * mrSges(3,1) + t122 * mrSges(3,2) + t146;
t143 = t158 + (g(2) * t100 + g(3) * t99) * mrSges(6,2);
t142 = m(6) * (-t154 * qJD(5) - t132 * t3 + t194);
t5 = -pkin(4) * t120 - t8;
t51 = Ifges(6,6) * qJD(5) + t121 * t156;
t87 = Ifges(6,4) * t178;
t52 = Ifges(6,1) * t179 + Ifges(6,5) * qJD(5) + t87;
t60 = t120 * t136 - t121 * t168;
t141 = t22 * mrSges(4,1) + t8 * mrSges(5,1) - t21 * mrSges(4,2) + mrSges(6,3) * t194 + t5 * t85 + (Ifges(6,1) * t61 + Ifges(6,4) * t60) * t200 + t136 * (Ifges(6,4) * t61 + Ifges(6,2) * t60) / 0.2e1 + t60 * t156 / 0.2e1 + t61 * (Ifges(6,1) * t132 + t188) / 0.2e1 - t51 * t168 / 0.2e1 + (t52 + t121 * (-Ifges(6,2) * t132 + t188)) * t167 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * t120 + (0.2e1 * Ifges(6,5) * t200 + Ifges(6,6) * t136) * qJDD(5) + (t148 * t121 / 0.2e1 + t202) * qJD(5);
t140 = t77 * mrSges(3,1) - t9 * mrSges(5,2) + Ifges(3,3) * t127 + t141;
t110 = -pkin(4) - t195;
t109 = pkin(8) + t196;
t79 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t179;
t64 = pkin(8) + t69;
t63 = -pkin(4) - t68;
t59 = t85 * t121;
t47 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t60;
t38 = -t117 * t170 + (qJD(2) * t152 - t134 * t169) * pkin(1);
t37 = t117 * t169 + (qJD(2) * t151 - t134 * t170) * pkin(1);
t33 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t30 = pkin(8) + t35;
t29 = -pkin(4) - t34;
t26 = t131 * t49 - t183;
t25 = t130 * t49 + t45;
t12 = t130 * t38 + t131 * t37;
t11 = t130 * t37 - t131 * t38;
t1 = [(t158 + (-m(6) * t13 - t79) * t12 + (m(6) * t201 + t162) * t30) * t132 + ((-t127 * t134 - t128 * t171) * mrSges(3,2) + (-qJD(2) * t128 * t134 + t127 * t138) * mrSges(3,1) + (-g(2) * t139 - g(3) * t135 + t134 * t78 + t138 * t77) * m(3)) * pkin(1) + t140 + (mrSges(4,1) * t38 - mrSges(5,1) * t11 - mrSges(4,2) * t37 - mrSges(5,2) * t12) * t121 + (mrSges(4,1) * t72 + mrSges(5,1) * t34 - mrSges(4,2) * t73 - mrSges(5,2) * t35) * t120 + (-t135 * mrSges(2,1) - mrSges(2,2) * t139 - m(6) * (t125 + t150) + t99 * t190 - m(5) * (t101 + t173) - m(4) * t173 + t145) * g(3) + (-mrSges(2,1) * t139 + t135 * mrSges(2,2) - m(6) * (t126 + t161) + t100 * t190 - m(5) * (t102 + t172) - m(4) * t172 + t144) * g(2) - t78 * mrSges(3,2) + t11 * t59 + m(6) * (t11 * t18 + t12 * t182 + t29 * t5 + (-t13 * t167 + t194) * t30) + t29 * t33 + m(4) * (t21 * t73 + t22 * t72 + t37 * t50 + t38 * t49) + m(5) * (-t11 * t23 + t12 * t24 + t34 * t8 + t35 * t9) + Ifges(2,3) * qJDD(1) + (t12 * t80 + t30 * t47 + (-t30 * t79 - t198) * qJD(5)) * t136; t140 + t63 * t33 + (t128 * t164 - t78) * mrSges(3,2) + t145 * g(3) + t144 * g(2) + (mrSges(5,1) * t68 - mrSges(5,2) * t69 + (mrSges(4,1) * t137 - mrSges(4,2) * t133) * pkin(2)) * t120 + t128 * mrSges(3,1) * t165 - t193 * t59 + (t64 * t47 - t192 * t80 + (-t64 * t79 - t198) * qJD(5)) * t136 + (t162 * t64 + t192 * t79 + t143) * t132 + (t192 * mrSges(5,2) + (-pkin(2) * t169 + t71) * mrSges(4,2) + t193 * mrSges(5,1) + (-pkin(2) * t170 - t70) * mrSges(4,1)) * t121 + t64 * t142 + (-t161 * g(2) - t150 * g(3) - t192 * t153 - t193 * t18 + t5 * t63) * m(6) + ((-t101 - t111) * g(3) + (-t102 - t112) * g(2) + t68 * t8 + t69 * t9 - t192 * t24 + t193 * t23) * m(5) + (-t111 * g(3) - t112 * g(2) - t49 * t70 - t50 * t71 + (t133 * t21 + t137 * t22 + (-t133 * t49 + t137 * t50) * qJD(3)) * pkin(2)) * m(4); t147 * g(3) + t146 * g(2) + t141 + (mrSges(4,1) * t50 + mrSges(5,1) * t25 + mrSges(4,2) * t49 + mrSges(5,2) * t26) * t121 + t110 * t33 - t25 * t59 + t120 * mrSges(5,1) * t195 + t109 * t142 + (-t120 * t196 - t9) * mrSges(5,2) + (t109 * t47 - t26 * t80 + (-t109 * t79 - t198) * qJD(5)) * t136 + (t109 * t162 + t26 * t79 + t143) * t132 + (-t166 * g(2) - t159 * g(3) + t110 * t5 - t153 * t26 - t18 * t25) * m(6) + (-t101 * g(3) - t102 * g(2) + t23 * t25 - t24 * t26 + (t130 * t9 + t131 * t8) * pkin(3)) * m(5); t132 * t47 + t136 * t48 + (-t132 * t79 + t136 * t80) * qJD(5) + (qJD(5) * t153 + t132 * t2 + t136 * t3 - g(1)) * m(6) + (qJDD(4) - g(1)) * m(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * qJDD(5) + g(1) * t85 - t13 * t80 + t14 * t79 + (t51 * t200 + (-t148 / 0.2e1 + t187 * t200) * t121 + t154 * mrSges(6,3) - (t52 + t87) * t136 / 0.2e1 - t202) * t121 + (g(2) * t99 - g(3) * t100) * t157;];
tau = t1;
