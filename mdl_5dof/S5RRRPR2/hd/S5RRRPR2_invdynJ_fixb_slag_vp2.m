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
% m [6x1]
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:37
% DurationCPUTime: 2.02s
% Computational Cost: add. (3151->287), mult. (5136->381), div. (0->0), fcn. (2732->16), ass. (0->144)
t124 = sin(qJ(5));
t192 = t124 / 0.2e1;
t198 = mrSges(5,2) - mrSges(6,3);
t122 = sin(pkin(9));
t123 = cos(pkin(9));
t129 = cos(qJ(3));
t125 = sin(qJ(3));
t169 = t122 * t125;
t176 = pkin(2) * qJD(3);
t130 = cos(qJ(2));
t126 = sin(qJ(2));
t166 = t126 * t129;
t143 = -t125 * t130 - t166;
t177 = pkin(1) * qJD(1);
t70 = t143 * t177;
t167 = t125 * t126;
t142 = t129 * t130 - t167;
t71 = t142 * t177;
t183 = t122 * t70 + t123 * t71 - (t123 * t129 - t169) * t176;
t128 = cos(qJ(5));
t182 = mrSges(6,1) * t128;
t197 = -mrSges(5,1) - t182;
t148 = mrSges(6,1) * t124 + mrSges(6,2) * t128;
t120 = qJD(1) + qJD(2);
t114 = qJD(3) + t120;
t156 = t126 * t177;
t155 = t130 * t177;
t84 = pkin(2) * t120 + t155;
t50 = t125 * t84 + t129 * t156;
t175 = t122 * t50;
t49 = -t125 * t156 + t129 * t84;
t42 = pkin(3) * t114 + t49;
t23 = t123 * t42 - t175;
t18 = -pkin(4) * t114 - t23;
t196 = t18 * t148 + qJD(5) * (Ifges(6,5) * t128 - Ifges(6,6) * t124) / 0.2e1;
t168 = t123 * t125;
t184 = t122 * t71 - t123 * t70 - (t122 * t129 + t168) * t176;
t45 = t123 * t50;
t24 = t122 * t42 + t45;
t19 = pkin(8) * t114 + t24;
t13 = qJD(4) * t128 - t124 * t19;
t14 = qJD(4) * t124 + t128 * t19;
t144 = -t124 * t13 + t128 * t14;
t121 = qJ(1) + qJ(2);
t117 = qJ(3) + t121;
t106 = pkin(9) + t117;
t97 = cos(t106);
t93 = t97 * pkin(8);
t96 = sin(t106);
t195 = t96 * mrSges(5,1) - m(6) * (-pkin(4) * t96 + t93);
t194 = m(3) * pkin(1);
t193 = m(5) + m(6);
t190 = mrSges(6,3) * t13;
t189 = pkin(1) * t130;
t116 = cos(t121);
t105 = pkin(2) * t116;
t108 = cos(t117);
t98 = pkin(3) * t108;
t188 = pkin(3) * t122;
t187 = pkin(3) * t123;
t119 = qJDD(1) + qJDD(2);
t113 = qJDD(3) + t119;
t77 = -qJD(2) * t156 + qJDD(1) * t189;
t62 = pkin(2) * t119 + t77;
t164 = qJD(2) * t130;
t78 = (qJD(1) * t164 + qJDD(1) * t126) * pkin(1);
t22 = -t50 * qJD(3) - t125 * t78 + t129 * t62;
t16 = pkin(3) * t113 + t22;
t21 = qJD(3) * t49 + t125 * t62 + t129 * t78;
t9 = t122 * t16 + t123 * t21;
t6 = pkin(8) * t113 + t9;
t2 = t13 * qJD(5) + qJDD(4) * t124 + t128 * t6;
t186 = t128 * t2;
t110 = pkin(2) + t189;
t72 = -pkin(1) * t167 + t129 * t110;
t67 = pkin(3) + t72;
t73 = pkin(1) * t166 + t110 * t125;
t35 = t122 * t67 + t123 * t73;
t109 = pkin(2) * t129 + pkin(3);
t69 = pkin(2) * t168 + t122 * t109;
t181 = mrSges(6,2) * t124;
t180 = Ifges(6,4) * t124;
t179 = Ifges(6,4) * t128;
t178 = Ifges(6,2) * t128;
t172 = qJD(5) * t14;
t171 = t114 * t124;
t170 = t114 * t128;
t131 = cos(qJ(1));
t118 = t131 * pkin(1);
t165 = t105 + t118;
t163 = qJD(3) * t125;
t162 = qJD(3) * t129;
t161 = qJD(5) * t124;
t160 = qJD(5) * t128;
t159 = m(4) + t193;
t158 = g(1) * mrSges(6,1) * t96;
t157 = t97 * pkin(4) + t96 * pkin(8) + t98;
t61 = t113 * t124 + t114 * t160;
t48 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t61;
t80 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t170;
t153 = -qJD(5) * t80 - t48;
t152 = t105 + t157;
t3 = qJDD(4) * t128 - t124 * t6 - t172;
t149 = (-t3 - t172) * mrSges(6,3);
t85 = t181 - t182;
t147 = t178 + t180;
t8 = -t122 * t21 + t123 * t16;
t34 = -t122 * t73 + t123 * t67;
t145 = t124 * t14 + t128 * t13;
t68 = -pkin(2) * t169 + t109 * t123;
t140 = t124 * (Ifges(6,1) * t128 - t180);
t107 = sin(t117);
t139 = -t108 * mrSges(4,1) + t107 * mrSges(4,2) + t197 * t97 + t198 * t96;
t138 = g(2) * mrSges(6,2) * t97 + t149;
t115 = sin(t121);
t137 = -t116 * mrSges(3,1) + t115 * mrSges(3,2) + t139;
t136 = t108 * mrSges(4,2) + (t193 * pkin(3) + mrSges(4,1)) * t107 - t96 * t181 + t198 * t97;
t135 = m(6) * (-t145 * qJD(5) - t124 * t3 + t186);
t5 = -pkin(4) * t113 - t8;
t51 = Ifges(6,6) * qJD(5) + t147 * t114;
t87 = Ifges(6,4) * t170;
t52 = Ifges(6,1) * t171 + Ifges(6,5) * qJD(5) + t87;
t60 = t113 * t128 - t114 * t161;
t134 = t22 * mrSges(4,1) + t8 * mrSges(5,1) - t21 * mrSges(4,2) + mrSges(6,3) * t186 + t5 * t85 + (Ifges(6,1) * t61 + Ifges(6,4) * t60) * t192 + t128 * (Ifges(6,4) * t61 + Ifges(6,2) * t60) / 0.2e1 + t60 * t147 / 0.2e1 + t61 * (Ifges(6,1) * t124 + t179) / 0.2e1 - t51 * t161 / 0.2e1 + (t52 + t114 * (-Ifges(6,2) * t124 + t179)) * t160 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * t113 + (0.2e1 * Ifges(6,5) * t192 + Ifges(6,6) * t128) * qJDD(5) + (t140 * t114 / 0.2e1 + t196) * qJD(5);
t133 = t77 * mrSges(3,1) - t9 * mrSges(5,2) + Ifges(3,3) * t119 + t134;
t132 = mrSges(3,2) * t116 + (t159 * pkin(2) + mrSges(3,1)) * t115 + t136;
t127 = sin(qJ(1));
t104 = -pkin(4) - t187;
t103 = pkin(8) + t188;
t79 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t171;
t64 = pkin(8) + t69;
t63 = -pkin(4) - t68;
t59 = t85 * t114;
t47 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t60;
t38 = -t110 * t163 + (t143 * qJD(2) - t126 * t162) * pkin(1);
t37 = t110 * t162 + (t142 * qJD(2) - t126 * t163) * pkin(1);
t33 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t30 = pkin(8) + t35;
t29 = -pkin(4) - t34;
t26 = t123 * t49 - t175;
t25 = t122 * t49 + t45;
t12 = t122 * t38 + t123 * t37;
t11 = t122 * t37 - t123 * t38;
t1 = [(-m(6) * t93 + mrSges(2,2) * t131 + (m(6) * pkin(4) - t197) * t96 + (mrSges(2,1) + (m(3) + t159) * pkin(1)) * t127 + t132) * g(1) + t133 + m(4) * (t21 * t73 + t22 * t72 + t37 * t50 + t38 * t49) + m(5) * (-t11 * t23 + t12 * t24 + t34 * t8 + t35 * t9) - t78 * mrSges(3,2) + t11 * t59 + t29 * t33 + (t127 * mrSges(2,2) - m(6) * (t118 + t152) + t97 * t181 - m(4) * t165 - m(5) * (t98 + t165) + (-mrSges(2,1) - t194) * t131 + t137) * g(2) + (-t12 * t79 + t153 * t30 + t149) * t124 + ((-t119 * t126 - t120 * t164) * mrSges(3,2) + (-qJD(2) * t120 * t126 + t119 * t130) * mrSges(3,1)) * pkin(1) + (t12 * t80 + t30 * t47 + (-t30 * t79 - t190) * qJD(5)) * t128 + (mrSges(4,1) * t38 - mrSges(5,1) * t11 - mrSges(4,2) * t37 - mrSges(5,2) * t12) * t114 + (mrSges(4,1) * t72 + mrSges(5,1) * t34 - mrSges(4,2) * t73 - mrSges(5,2) * t35) * t113 + (t126 * t78 + t130 * t77) * t194 + Ifges(2,3) * qJDD(1) + t30 * t135 + m(6) * (t11 * t18 + t144 * t12 + t29 * t5); (t132 + t195) * g(1) + t133 + t137 * g(2) - t184 * t59 + t63 * t33 + (t183 * mrSges(5,2) + (-pkin(2) * t162 + t71) * mrSges(4,2) + t184 * mrSges(5,1) + (-pkin(2) * t163 - t70) * mrSges(4,1)) * t114 + (t153 * t64 + t183 * t79 + t138) * t124 + t120 * mrSges(3,1) * t156 + t64 * t135 + (t120 * t155 - t78) * mrSges(3,2) + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + (mrSges(4,1) * t129 - mrSges(4,2) * t125) * pkin(2)) * t113 + (t158 + t64 * t47 - t183 * t80 + (-t64 * t79 - t190) * qJD(5)) * t128 + (-t152 * g(2) - t183 * t144 - t18 * t184 + t5 * t63) * m(6) + ((-t98 - t105) * g(2) + t68 * t8 + t69 * t9 - t183 * t24 + t184 * t23) * m(5) + (-t105 * g(2) - t49 * t70 - t50 * t71 + (t125 * t21 + t129 * t22 + (-t125 * t49 + t129 * t50) * qJD(3)) * pkin(2)) * m(4); t134 + t139 * g(2) + (mrSges(4,1) * t50 + mrSges(5,1) * t25 + mrSges(4,2) * t49 + mrSges(5,2) * t26) * t114 + t104 * t33 - t25 * t59 + t113 * mrSges(5,1) * t187 + t103 * t135 + (-t113 * t188 - t9) * mrSges(5,2) + (t136 + t195) * g(1) + (t158 + t103 * t47 - t26 * t80 + (-t103 * t79 - t190) * qJD(5)) * t128 + (t153 * t103 + t26 * t79 + t138) * t124 + (-t157 * g(2) + t104 * t5 - t144 * t26 - t18 * t25) * m(6) + (-t98 * g(2) + t23 * t25 - t24 * t26 + (t122 * t9 + t123 * t8) * pkin(3)) * m(5); t124 * t47 + t128 * t48 + (-t124 * t79 + t128 * t80) * qJD(5) + (t144 * qJD(5) + t124 * t2 + t128 * t3 - g(3)) * m(6) + (qJDD(4) - g(3)) * m(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * qJDD(5) + g(3) * t85 - t13 * t80 + t14 * t79 + (t51 * t192 + (-t140 / 0.2e1 + t178 * t192) * t114 + t145 * mrSges(6,3) - (t52 + t87) * t128 / 0.2e1 - t196) * t114 + (g(1) * t97 + g(2) * t96) * t148;];
tau = t1;
