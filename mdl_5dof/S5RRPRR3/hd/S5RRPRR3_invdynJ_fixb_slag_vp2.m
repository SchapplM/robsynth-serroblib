% Calculate vector of inverse dynamics joint torques for
% S5RRPRR3
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:25
% DurationCPUTime: 1.78s
% Computational Cost: add. (2960->279), mult. (4825->363), div. (0->0), fcn. (2646->16), ass. (0->141)
t126 = sin(qJ(5));
t193 = t126 / 0.2e1;
t127 = sin(qJ(4));
t131 = cos(qJ(4));
t124 = sin(pkin(9));
t132 = cos(qJ(2));
t125 = cos(pkin(9));
t128 = sin(qJ(2));
t164 = t125 * t128;
t141 = pkin(1) * (-t124 * t132 - t164);
t62 = qJD(1) * t141;
t165 = t124 * t128;
t140 = pkin(1) * (t125 * t132 - t165);
t64 = qJD(1) * t140;
t188 = pkin(2) * t125;
t106 = pkin(3) + t188;
t189 = pkin(2) * t124;
t67 = t106 * t131 - t127 * t189;
t184 = -t67 * qJD(4) + t127 * t62 + t131 * t64;
t130 = cos(qJ(5));
t181 = mrSges(6,1) * t130;
t202 = -mrSges(5,1) - t181;
t122 = qJD(1) + qJD(2);
t176 = qJD(1) * pkin(1);
t154 = t128 * t176;
t80 = pkin(2) * t122 + t132 * t176;
t42 = -t124 * t154 + t125 * t80;
t37 = pkin(3) * t122 + t42;
t43 = t124 * t80 + t125 * t154;
t17 = -t127 * t43 + t131 * t37;
t149 = mrSges(6,1) * t126 + mrSges(6,2) * t130;
t115 = qJD(4) + t122;
t15 = -pkin(4) * t115 - t17;
t201 = t15 * t149 + qJD(5) * (Ifges(6,5) * t130 - Ifges(6,6) * t126) / 0.2e1;
t18 = t127 * t37 + t131 * t43;
t16 = pkin(8) * t115 + t18;
t14 = qJD(3) * t126 + t130 * t16;
t169 = qJD(5) * t14;
t121 = qJDD(1) + qJDD(2);
t114 = qJDD(4) + t121;
t190 = pkin(1) * t132;
t73 = -qJD(2) * t154 + qJDD(1) * t190;
t58 = pkin(2) * t121 + t73;
t161 = qJD(2) * t132;
t74 = (qJD(1) * t161 + qJDD(1) * t128) * pkin(1);
t32 = -t124 * t74 + t125 * t58;
t22 = pkin(3) * t121 + t32;
t33 = t124 * t58 + t125 * t74;
t8 = t17 * qJD(4) + t127 * t22 + t131 * t33;
t5 = pkin(8) * t114 + t8;
t3 = qJDD(3) * t130 - t126 * t5 - t169;
t200 = -t3 - t169;
t69 = t127 * t106 + t131 * t189;
t183 = t69 * qJD(4) - t127 * t64 + t131 * t62;
t159 = qJD(5) * t130;
t57 = t114 * t126 + t115 * t159;
t41 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t57;
t166 = t115 * t130;
t76 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t166;
t199 = -qJD(5) * t76 - t41;
t13 = qJD(3) * t130 - t126 * t16;
t2 = qJD(5) * t13 + qJDD(3) * t126 + t130 * t5;
t187 = t130 * t2;
t198 = -t13 * t159 + t187;
t180 = mrSges(6,2) * t126;
t123 = qJ(1) + qJ(2);
t116 = pkin(9) + t123;
t109 = qJ(4) + t116;
t97 = sin(t109);
t98 = cos(t109);
t197 = t97 * mrSges(5,2) + t98 * t180;
t173 = t130 * t14;
t145 = -t126 * t13 + t173;
t9 = -t18 * qJD(4) - t127 * t33 + t131 * t22;
t196 = m(4) + m(5);
t195 = g(2) * t97;
t194 = g(3) * t98;
t191 = mrSges(6,3) * t13;
t117 = sin(t123);
t107 = pkin(2) * t117;
t118 = cos(t123);
t108 = pkin(2) * t118;
t186 = t3 * t126;
t110 = pkin(2) + t190;
t66 = -pkin(1) * t165 + t125 * t110;
t59 = pkin(3) + t66;
t68 = pkin(1) * t164 + t110 * t124;
t29 = t127 * t59 + t131 * t68;
t182 = t98 * pkin(4) + t97 * pkin(8);
t179 = Ifges(6,4) * t126;
t178 = Ifges(6,4) * t130;
t177 = Ifges(6,2) * t130;
t104 = sin(t116);
t95 = pkin(3) * t104;
t171 = t95 + t107;
t105 = cos(t116);
t96 = pkin(3) * t105;
t170 = t96 + t108;
t167 = t115 * t126;
t129 = sin(qJ(1));
t119 = t129 * pkin(1);
t163 = t107 + t119;
t133 = cos(qJ(1));
t120 = t133 * pkin(1);
t162 = t108 + t120;
t160 = qJD(5) * t126;
t158 = -t98 * mrSges(5,2) + t202 * t97;
t157 = -t97 * mrSges(6,3) + t202 * t98;
t156 = t97 * t180;
t151 = t170 + t182;
t150 = t200 * mrSges(6,3);
t81 = t180 - t181;
t148 = t177 + t179;
t146 = t126 * t14 + t13 * t130;
t28 = -t127 * t68 + t131 * t59;
t90 = t97 * pkin(4);
t143 = -pkin(8) * t98 + t171 + t90;
t139 = t126 * (Ifges(6,1) * t130 - t179);
t138 = -t117 * mrSges(3,1) - t104 * mrSges(4,1) - t118 * mrSges(3,2) - t105 * mrSges(4,2) + t98 * mrSges(6,3) + t158;
t137 = -t118 * mrSges(3,1) - t105 * mrSges(4,1) + t117 * mrSges(3,2) + t104 * mrSges(4,2) + t157;
t136 = -t146 * qJD(5) - t186;
t45 = Ifges(6,6) * qJD(5) + t148 * t115;
t84 = Ifges(6,4) * t166;
t46 = Ifges(6,1) * t167 + Ifges(6,5) * qJD(5) + t84;
t56 = t114 * t130 - t115 * t160;
t6 = -pkin(4) * t114 - t9;
t135 = t9 * mrSges(5,1) + mrSges(6,3) * t187 + t6 * t81 + Ifges(5,3) * t114 + (Ifges(6,1) * t57 + Ifges(6,4) * t56) * t193 + t130 * (Ifges(6,4) * t57 + Ifges(6,2) * t56) / 0.2e1 + t56 * t148 / 0.2e1 + t57 * (Ifges(6,1) * t126 + t178) / 0.2e1 - t45 * t160 / 0.2e1 + (t46 + t115 * (-Ifges(6,2) * t126 + t178)) * t159 / 0.2e1 + (0.2e1 * Ifges(6,5) * t193 + Ifges(6,6) * t130) * qJDD(5) + (t139 * t115 / 0.2e1 + t201) * qJD(5);
t134 = t73 * mrSges(3,1) + t32 * mrSges(4,1) - t74 * mrSges(3,2) + t135 + (Ifges(3,3) + Ifges(4,3)) * t121;
t75 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t167;
t65 = qJD(2) * t140;
t63 = qJD(2) * t141;
t61 = pkin(8) + t69;
t60 = -pkin(4) - t67;
t55 = t81 * t115;
t40 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t56;
t30 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t24 = pkin(8) + t29;
t23 = -pkin(4) - t28;
t11 = t29 * qJD(4) + t127 * t65 - t131 * t63;
t10 = t28 * qJD(4) + t127 * t63 + t131 * t65;
t1 = [(-t11 * t115 + t114 * t28) * mrSges(5,1) + (t121 * t66 + t122 * t63) * mrSges(4,1) + (t10 * t76 + t24 * t40 + (-t24 * t75 - t191) * qJD(5)) * t130 + m(4) * (t32 * t66 + t33 * t68 + t42 * t63 + t43 * t65) + m(5) * (t10 * t18 - t11 * t17 + t28 * t9 + t29 * t8) + (-t10 * t115 - t114 * t29 - t8) * mrSges(5,2) + (-t121 * t68 - t122 * t65 - t33) * mrSges(4,2) + m(6) * (t10 * t173 + t11 * t15 + t198 * t24 + t23 * t6) + t11 * t55 + t23 * t30 + (-t129 * mrSges(2,1) - t133 * mrSges(2,2) - m(6) * (t119 + t143) + t156 - m(4) * t163 - m(5) * (t95 + t163) + t138) * g(3) + (-t133 * mrSges(2,1) + t129 * mrSges(2,2) - m(6) * (t120 + t151) - m(4) * t162 - m(5) * (t96 + t162) + t137 + t197) * g(2) + t134 + ((-t121 * t128 - t122 * t161) * mrSges(3,2) + (-qJD(2) * t122 * t128 + t121 * t132) * mrSges(3,1) + (-g(2) * t133 - g(3) * t129 + t128 * t74 + t132 * t73) * m(3)) * pkin(1) + (t150 + (-m(6) * t13 - t75) * t10 + (m(6) * t200 + t199) * t24) * t126 + Ifges(2,3) * qJDD(1); t137 * g(2) + t183 * t55 + t121 * mrSges(4,1) * t188 + (t184 * t75 + t199 * t61 + t150 + (g(2) * t98 + g(3) * t97) * mrSges(6,2)) * t126 + (t61 * t40 - t184 * t76 + (-t61 * t75 - t191) * qJD(5)) * t130 + (-t69 * t114 + t184 * t115 + t195 - t8) * mrSges(5,2) + (-t62 * mrSges(4,1) + t64 * mrSges(4,2) + (mrSges(3,1) * t128 + mrSges(3,2) * t132) * t176) * t122 + (t67 * t114 - t183 * t115) * mrSges(5,1) + (-t121 * t189 - t33) * mrSges(4,2) + t138 * g(3) + t60 * t30 + t134 + (-t170 * g(2) - t171 * g(3) - t183 * t17 - t184 * t18 + t67 * t9 + t69 * t8) * m(5) + (-t108 * g(2) + (t124 * t33 + t125 * t32) * pkin(2) - t107 * g(3) - t42 * t62 - t43 * t64) * m(4) + (-t151 * g(2) + (t136 + t187) * t61 - t143 * g(3) + t6 * t60 + t183 * t15 - t184 * t145) * m(6); m(6) * (t145 * qJD(5) + t126 * t2 + t130 * t3) + t76 * t159 + t126 * t40 - t75 * t160 + t130 * t41 + t196 * qJDD(3) + (-m(6) - t196) * g(1); (t136 + t194) * mrSges(6,3) + (t115 * mrSges(5,1) - t55) * t18 + (t115 * mrSges(5,2) + t126 * t75 - t130 * t76) * t17 + (t156 + t158) * g(3) + (t157 + t197) * g(2) - pkin(4) * t30 - t8 * mrSges(5,2) + t135 + (-t126 * t41 + t130 * t40 + (-t126 * t76 - t130 * t75) * qJD(5)) * pkin(8) + (-pkin(4) * t6 - t90 * g(3) - t182 * g(2) - t145 * t17 - t15 * t18 + (-t14 * t160 - t186 + t194 + t198) * pkin(8)) * m(6); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t57 + Ifges(6,6) * t56 + Ifges(6,3) * qJDD(5) + g(1) * t81 - t13 * t76 + t14 * t75 + (t45 * t193 + (-t139 / 0.2e1 + t177 * t193) * t115 + t146 * mrSges(6,3) - (t46 + t84) * t130 / 0.2e1 - t201) * t115 + (-t194 + t195) * t149;];
tau = t1;
