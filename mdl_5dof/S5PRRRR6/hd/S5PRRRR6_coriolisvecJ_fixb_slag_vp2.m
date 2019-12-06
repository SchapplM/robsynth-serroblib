% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:34
% EndTime: 2019-12-05 17:09:39
% DurationCPUTime: 1.59s
% Computational Cost: add. (2683->251), mult. (5124->378), div. (0->0), fcn. (3387->8), ass. (0->139)
t120 = qJD(2) + qJD(3);
t123 = sin(qJ(5));
t124 = sin(qJ(4));
t127 = cos(qJ(5));
t128 = cos(qJ(4));
t141 = t123 * t124 - t127 * t128;
t89 = t141 * t120;
t196 = t89 / 0.2e1;
t195 = -mrSges(5,1) * t128 + mrSges(5,2) * t124 - mrSges(4,1);
t188 = -pkin(8) - pkin(7);
t156 = qJD(4) * t188;
t104 = t124 * t156;
t105 = t128 * t156;
t110 = t188 * t124;
t118 = t128 * pkin(8);
t111 = pkin(7) * t128 + t118;
t74 = t110 * t127 - t111 * t123;
t130 = cos(qJ(2));
t166 = qJD(1) * t130;
t112 = qJD(2) * pkin(2) + t166;
t125 = sin(qJ(3));
t126 = sin(qJ(2));
t167 = qJD(1) * t126;
t114 = t125 * t167;
t129 = cos(qJ(3));
t85 = t112 * t129 - t114;
t194 = qJD(5) * t74 + t104 * t127 + t105 * t123 + t141 * t85;
t102 = t123 * t128 + t124 * t127;
t75 = t110 * t123 + t111 * t127;
t193 = -qJD(5) * t75 + t102 * t85 - t104 * t123 + t105 * t127;
t115 = pkin(2) * t125 + pkin(7);
t182 = -pkin(8) - t115;
t98 = t182 * t124;
t170 = t115 * t128;
t99 = t118 + t170;
t64 = t123 * t98 + t127 * t99;
t150 = qJD(4) * t182;
t164 = qJD(3) * t129;
t158 = pkin(2) * t164;
t76 = t124 * t150 + t128 * t158;
t77 = -t124 * t158 + t128 * t150;
t97 = t129 * t166 - t114;
t192 = -qJD(5) * t64 + t102 * t97 - t123 * t76 + t127 * t77;
t63 = -t123 * t99 + t127 * t98;
t191 = qJD(5) * t63 + t123 * t77 + t127 * t76 + t141 * t97;
t103 = t125 * t130 + t126 * t129;
t59 = t141 * t103;
t119 = qJD(4) + qJD(5);
t86 = t112 * t125 + t129 * t167;
t79 = pkin(7) * t120 + t86;
t149 = (t124 ^ 2 + t128 ^ 2) * t79;
t90 = t102 * t120;
t189 = t90 / 0.2e1;
t187 = mrSges(6,3) * t89;
t186 = Ifges(6,4) * t90;
t185 = pkin(2) * t129;
t161 = qJD(4) * t128;
t101 = t125 * t126 - t129 * t130;
t137 = t101 * qJD(2);
t165 = qJD(3) * t125;
t55 = t112 * t164 + (-t126 * t165 - t137) * qJD(1);
t176 = t124 * t55;
t30 = -t161 * t79 - t176;
t184 = t30 * mrSges(5,3);
t183 = t90 * mrSges(6,3);
t181 = Ifges(5,4) * t124;
t138 = t103 * qJD(2);
t56 = t112 * t165 + (t126 * t164 + t138) * qJD(1);
t179 = t101 * t56;
t154 = pkin(8) * t120 + t79;
t62 = t154 * t128;
t178 = t123 * t62;
t177 = t124 * t30;
t175 = t127 * t62;
t162 = qJD(4) * t124;
t53 = t128 * t55;
t29 = -t162 * t79 + t53;
t174 = t128 * t29;
t172 = Ifges(5,5) * qJD(4);
t171 = Ifges(5,6) * qJD(4);
t169 = t120 * t124;
t168 = t120 * t128;
t163 = qJD(4) * t115;
t160 = qJD(5) * t123;
t159 = qJD(5) * t127;
t157 = pkin(4) * t162;
t117 = -pkin(4) * t128 - pkin(3);
t155 = t120 * t162;
t151 = t195 * t120;
t52 = mrSges(6,1) * t89 + mrSges(6,2) * t90;
t148 = -t151 - t52;
t61 = t154 * t124;
t57 = qJD(4) * pkin(4) - t61;
t17 = t127 * t57 - t178;
t18 = t123 * t57 + t175;
t145 = qJD(4) * t154;
t19 = -t124 * t145 + t53;
t20 = -t128 * t145 - t176;
t4 = -qJD(5) * t18 - t123 * t19 + t127 * t20;
t66 = t119 * t141;
t146 = -t4 * t102 + t17 * t66;
t143 = t174 - t177;
t107 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t169;
t108 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t168;
t142 = t107 * t128 + t108 * t124;
t140 = (Ifges(5,2) * t128 + t181) * t120;
t139 = (mrSges(5,1) * t124 + mrSges(5,2) * t128) * qJD(4);
t136 = -t120 * mrSges(4,2) - t124 * t107 + t128 * t108;
t3 = qJD(5) * t17 + t123 * t20 + t127 * t19;
t43 = -Ifges(6,2) * t89 + Ifges(6,6) * t119 + t186;
t82 = Ifges(6,4) * t89;
t44 = Ifges(6,1) * t90 + Ifges(6,5) * t119 - t82;
t47 = t66 * t120;
t67 = t119 * t102;
t48 = t67 * t120;
t65 = t117 * t120 - t85;
t133 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t17 * t187 + t43 * t189 - t65 * (mrSges(6,1) * t90 - mrSges(6,2) * t89) - t90 * (-Ifges(6,1) * t89 - t186) / 0.2e1 - t119 * (-Ifges(6,5) * t89 - Ifges(6,6) * t90) / 0.2e1 - Ifges(6,6) * t48 - Ifges(6,5) * t47 + (-Ifges(6,2) * t90 + t44 - t82) * t196;
t41 = pkin(4) * t155 + t56;
t78 = -pkin(3) * t120 - t85;
t87 = t140 + t171;
t113 = Ifges(5,4) * t168;
t88 = Ifges(5,1) * t169 + t113 + t172;
t132 = mrSges(5,3) * t174 + t78 * t139 + t119 * (-Ifges(6,5) * t66 - Ifges(6,6) * t67) / 0.2e1 + t41 * (mrSges(6,1) * t141 + mrSges(6,2) * t102) + qJD(4) ^ 2 * (Ifges(5,5) * t128 - Ifges(5,6) * t124) / 0.2e1 - t66 * t44 / 0.2e1 - t67 * t43 / 0.2e1 + t65 * (mrSges(6,1) * t67 - mrSges(6,2) * t66) + (Ifges(5,1) * t128 - t181) * t155 + t195 * t56 - (t140 + t87) * t162 / 0.2e1 + (t141 * t48 + t67 * t196) * Ifges(6,2) + (-t102 * t47 - t66 * t189) * Ifges(6,1) + (-t141 * t3 - t18 * t67) * mrSges(6,3) + (-t102 * t48 + t141 * t47 - t67 * t189 + t66 * t196) * Ifges(6,4) + (t88 + (0.3e1 * Ifges(5,4) * t128 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t124) * t120) * t161 / 0.2e1;
t116 = -pkin(3) - t185;
t109 = t117 - t185;
t106 = pkin(2) * t165 + t157;
t96 = t103 * qJD(1);
t91 = t120 * t139;
t71 = mrSges(6,1) * t119 - t183;
t70 = -mrSges(6,2) * t119 - t187;
t69 = qJD(3) * t103 + t138;
t68 = -qJD(3) * t101 - t137;
t58 = t102 * t103;
t24 = -t127 * t61 - t178;
t23 = t123 * t61 - t175;
t14 = mrSges(6,1) * t48 - mrSges(6,2) * t47;
t6 = -t102 * t68 + t119 * t59;
t5 = -t103 * t67 - t141 * t68;
t1 = [t5 * t70 + t6 * t71 + (-mrSges(3,1) * t126 - mrSges(3,2) * t130) * qJD(2) ^ 2 + (t14 + t91) * t101 - t142 * t103 * qJD(4) + (-t47 * t58 + t48 * t59) * mrSges(6,3) - t148 * t69 + t136 * t68 + m(4) * (t103 * t55 + t68 * t86 - t69 * t85 + t179) + m(5) * (t103 * t143 + t149 * t68 + t69 * t78 + t179) + m(6) * (t101 * t41 + t17 * t6 + t18 * t5 - t3 * t59 - t4 * t58 + t65 * t69); -m(5) * (t97 * t149 + t78 * t96) + t148 * t96 + (-t107 * t163 - t97 * t108) * t128 + (t97 * t107 - t108 * t163 - t184) * t124 + (t63 * t47 - t64 * t48 + t146) * mrSges(6,3) + t116 * t91 + t106 * t52 + t109 * t14 + m(5) * (-t115 * t177 + t116 * t56 + t29 * t170) + (m(4) * (t125 * t55 - t129 * t56) + (t151 * t125 + t136 * t129 + m(5) * (t125 * t78 + t129 * t149) + m(4) * (-t125 * t85 + t129 * t86)) * qJD(3)) * pkin(2) + (t97 * t120 - t55) * mrSges(4,2) - m(4) * (-t85 * t96 + t86 * t97) + t192 * t71 + t191 * t70 + t132 + (t109 * t41 + t3 * t64 + t4 * t63 + (t106 - t96) * t65 + t191 * t18 + t192 * t17) * m(6); (t74 * t47 - t75 * t48 + t146) * mrSges(6,3) + (t85 * t120 - t55) * mrSges(4,2) + t117 * t14 - pkin(3) * t91 + (-t184 + t85 * t107 + (pkin(4) * t52 - pkin(7) * t108) * qJD(4)) * t124 + t132 + t194 * t70 + t193 * t71 + (-qJD(4) * pkin(7) * t107 - t85 * t108) * t128 + t148 * t86 + (t117 * t41 + t3 * t75 + t4 * t74 + (t157 - t86) * t65 + t194 * t18 + t193 * t17) * m(6) + (-pkin(3) * t56 + pkin(7) * t143 - t149 * t85 - t78 * t86) * m(5); -m(6) * (t17 * t23 + t18 * t24) + t142 * t79 - t24 * t70 - t23 * t71 - t29 * mrSges(5,2) + t30 * mrSges(5,1) + t18 * t183 + (t70 * t159 - t71 * t160 + m(6) * (t123 * t3 + t127 * t4 + t159 * t18 - t160 * t17) + (-t123 * t48 + t127 * t47) * mrSges(6,3)) * pkin(4) + ((t172 / 0.2e1 - t78 * mrSges(5,2) - t113 / 0.2e1 - t88 / 0.2e1) * t128 + (-t171 / 0.2e1 - t78 * mrSges(5,1) + t87 / 0.2e1 + (t181 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t128) * t120 + (-m(6) * t65 - t52) * pkin(4)) * t124) * t120 + t133; -t17 * t70 + (t71 + t183) * t18 + t133;];
tauc = t1(:);
