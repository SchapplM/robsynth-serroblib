% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:48
% EndTime: 2019-03-08 19:17:53
% DurationCPUTime: 2.47s
% Computational Cost: add. (2286->299), mult. (5500->416), div. (0->0), fcn. (3722->10), ass. (0->159)
t83 = sin(pkin(6));
t144 = qJD(1) * t83;
t91 = cos(qJ(2));
t130 = t91 * t144;
t84 = cos(pkin(11));
t121 = t84 * t130;
t88 = sin(qJ(2));
t131 = t88 * t144;
t82 = sin(pkin(11));
t73 = t82 * t131;
t58 = -t73 + t121;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t117 = pkin(5) * t90 + pkin(9) * t87;
t66 = t117 * qJD(5) + qJD(4);
t194 = -t58 + t66;
t89 = cos(qJ(6));
t160 = Ifges(7,4) * t89;
t86 = sin(qJ(6));
t109 = -Ifges(7,2) * t86 + t160;
t161 = Ifges(7,4) * t86;
t111 = Ifges(7,1) * t89 - t161;
t112 = mrSges(7,1) * t86 + mrSges(7,2) * t89;
t72 = qJD(2) * pkin(2) + t130;
t43 = t72 * t84 - t73;
t114 = qJD(4) - t43;
t33 = (-pkin(3) - pkin(8)) * qJD(2) + t114;
t85 = cos(pkin(6));
t77 = qJD(1) * t85 + qJD(3);
t25 = t33 * t87 + t77 * t90;
t23 = qJD(5) * pkin(9) + t25;
t100 = pkin(5) * t87 - pkin(9) * t90 + qJ(4);
t44 = t84 * t131 + t72 * t82;
t27 = t100 * qJD(2) + t44;
t5 = -t23 * t86 + t27 * t89;
t6 = t23 * t89 + t27 * t86;
t115 = t5 * t89 + t6 * t86;
t158 = Ifges(7,6) * t86;
t159 = Ifges(7,5) * t89;
t170 = t89 / 0.2e1;
t171 = -t86 / 0.2e1;
t140 = qJD(5) * t86;
t141 = qJD(2) * t90;
t68 = t89 * t141 + t140;
t174 = t68 / 0.2e1;
t24 = t33 * t90 - t77 * t87;
t22 = -qJD(5) * pkin(5) - t24;
t162 = Ifges(7,4) * t68;
t139 = qJD(5) * t89;
t67 = -t86 * t141 + t139;
t142 = qJD(2) * t87;
t80 = qJD(6) + t142;
t29 = Ifges(7,2) * t67 + Ifges(7,6) * t80 + t162;
t65 = Ifges(7,4) * t67;
t30 = Ifges(7,1) * t68 + Ifges(7,5) * t80 + t65;
t193 = -t115 * mrSges(7,3) + (-t158 + t159) * t80 / 0.2e1 + t109 * t67 / 0.2e1 + t111 * t174 + t22 * t112 + t30 * t170 + t29 * t171;
t52 = -mrSges(7,2) * t80 + mrSges(7,3) * t67;
t53 = mrSges(7,1) * t80 - mrSges(7,3) * t68;
t184 = -m(7) * t115 - t86 * t52 - t89 * t53;
t134 = qJD(2) * qJD(5);
t125 = t90 * t134;
t133 = qJD(5) * qJD(6);
t138 = qJD(6) * t90;
t50 = t89 * t133 + (-t86 * t138 - t87 * t139) * qJD(2);
t31 = mrSges(7,1) * t125 - mrSges(7,3) * t50;
t51 = -t86 * t133 + (-t89 * t138 + t87 * t140) * qJD(2);
t32 = -mrSges(7,2) * t125 + mrSges(7,3) * t51;
t192 = t184 * qJD(6) - t86 * t31 + t89 * t32;
t126 = -pkin(2) * t84 - pkin(3);
t79 = -pkin(8) + t126;
t129 = qJD(5) * t79 * t90;
t150 = t86 * t87;
t149 = t87 * t89;
t169 = pkin(2) * t82;
t64 = t100 + t169;
t35 = t79 * t149 + t64 * t86;
t60 = (t82 * t91 + t84 * t88) * t83;
t55 = qJD(1) * t60;
t191 = -t35 * qJD(6) - t86 * t129 + t55 * t150 + t194 * t89;
t34 = -t79 * t150 + t64 * t89;
t190 = t34 * qJD(6) + t89 * t129 - t55 * t149 + t194 * t86;
t146 = Ifges(6,5) * qJD(5);
t163 = Ifges(6,4) * t87;
t187 = qJD(2) / 0.2e1;
t37 = qJD(2) * qJ(4) + t44;
t189 = t37 * mrSges(6,2) + t146 / 0.2e1 + (t90 * Ifges(6,1) - t163) * t187 - t24 * mrSges(6,3) + t193;
t143 = qJD(2) * t83;
t152 = t82 * t88;
t119 = t143 * t152;
t70 = qJD(1) * t119;
t26 = -t70 + (t66 + t121) * qJD(2);
t137 = t24 * qJD(5);
t56 = qJD(2) * t60;
t48 = qJD(1) * t56;
t7 = t48 * t87 + t137;
t1 = t5 * qJD(6) + t26 * t86 + t7 * t89;
t2 = -t6 * qJD(6) + t26 * t89 - t7 * t86;
t116 = t1 * t89 - t2 * t86;
t186 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t50 + Ifges(7,6) * t51;
t127 = mrSges(6,3) * t141;
t147 = qJD(5) * mrSges(6,1) + mrSges(7,1) * t67 - mrSges(7,2) * t68 - t127;
t185 = m(7) * t22 - t147;
t136 = t25 * qJD(5);
t21 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t8 = -t90 * t48 + t136;
t183 = -m(6) * (-t8 + t136) + m(7) * (-t6 * t139 + t5 * t140 + t8) + t21;
t181 = 0.2e1 * m(6);
t180 = t50 / 0.2e1;
t179 = t51 / 0.2e1;
t178 = -t55 / 0.2e1;
t177 = -t67 / 0.2e1;
t175 = -t68 / 0.2e1;
t173 = -t80 / 0.2e1;
t151 = t84 * t91;
t59 = (-t151 + t152) * t83;
t102 = t59 * t90 - t85 * t87;
t166 = t102 * t8;
t157 = t37 * t58;
t156 = t48 * t59;
t148 = mrSges(4,1) - mrSges(5,2);
t145 = Ifges(6,6) * qJD(5);
t135 = -t58 + qJD(4);
t132 = -t79 * t8 / 0.2e1;
t128 = mrSges(6,3) * t142;
t124 = -t146 / 0.2e1;
t123 = -t145 / 0.2e1;
t120 = t143 * t151;
t113 = mrSges(7,1) * t89 - mrSges(7,2) * t86;
t110 = Ifges(7,1) * t86 + t160;
t108 = Ifges(7,2) * t89 + t161;
t106 = Ifges(7,5) * t86 + Ifges(7,6) * t89;
t41 = -t70 + (qJD(4) + t121) * qJD(2);
t57 = -t119 + t120;
t104 = t37 * t57 + t41 * t60;
t40 = t59 * t87 + t85 * t90;
t18 = t40 * t89 + t60 * t86;
t17 = -t40 * t86 + t60 * t89;
t81 = qJ(4) + t169;
t101 = t37 * qJD(4) + t41 * t81;
t75 = -qJD(5) * mrSges(6,2) - t128;
t99 = t89 * t52 - t86 * t53 + t75;
t95 = t6 * mrSges(7,2) - t80 * Ifges(7,3) - t68 * Ifges(7,5) - t67 * Ifges(7,6) + t145 / 0.2e1 + (Ifges(6,4) * t90 - t87 * Ifges(6,2)) * t187 - t37 * mrSges(6,1) - t5 * mrSges(7,1);
t94 = m(6) * (t7 - t137) + m(7) * (qJD(5) * t22 + t116) + t192;
t92 = qJD(2) ^ 2;
t78 = Ifges(7,3) * t125;
t71 = t117 * qJD(2);
t69 = (mrSges(6,1) * t87 + mrSges(6,2) * t90) * qJD(2);
t63 = (mrSges(6,1) * t90 - mrSges(6,2) * t87) * t134;
t49 = qJD(1) * t120 - t70;
t36 = -qJD(2) * pkin(3) + t114;
t16 = t40 * qJD(5) - t56 * t90;
t15 = t102 * qJD(5) + t56 * t87;
t14 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + Ifges(7,5) * t125;
t13 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + Ifges(7,6) * t125;
t10 = t24 * t89 + t71 * t86;
t9 = -t24 * t86 + t71 * t89;
t4 = t17 * qJD(6) + t15 * t89 + t57 * t86;
t3 = -t18 * qJD(6) - t15 * t86 + t57 * t89;
t11 = [t15 * t75 + t17 * t31 + t18 * t32 - t102 * t21 + t3 * t53 + t4 * t52 + t57 * t69 + t60 * t63 + (-mrSges(3,1) * t88 - mrSges(3,2) * t91) * t92 * t83 - t147 * t16 + m(4) * (-t43 * t56 + t44 * t57 + t49 * t60 + t156) + m(6) * (t15 * t25 - t16 * t24 + t40 * t7 + t104 - t166) + m(7) * (t1 * t18 + t16 * t22 + t17 * t2 + t3 * t5 + t4 * t6 - t166) + m(5) * (t36 * t56 + t104 + t156) + ((-mrSges(4,2) + mrSges(5,3)) * t57 - t148 * t56 + (t102 * t87 - t40 * t90) * qJD(5) * mrSges(6,3)) * qJD(2); -t49 * mrSges(4,2) + t41 * mrSges(5,3) + t34 * t31 + t35 * t32 + t81 * t63 + t135 * t69 + t191 * t53 + t190 * t52 - t148 * t48 + (t58 * mrSges(4,2) + t135 * mrSges(5,3) + t148 * t55) * qJD(2) + m(6) * t101 - m(6) * t157 + (t78 / 0.2e1 - t55 * t75 - t7 * mrSges(6,3) + t41 * mrSges(6,1) + (t7 * t79 / 0.2e1 + t25 * t178) * t181 + (0.3e1 / 0.2e1 * Ifges(6,4) * t142 + t124 + (-m(6) * t24 + t185) * t79 - t189) * qJD(5) + t186) * t87 + (t109 * t179 + t41 * mrSges(6,2) + t14 * t170 + t13 * t171 - t79 * t21 + t111 * t180 + (mrSges(6,3) + t112) * t8 + (-t1 * t86 - t2 * t89) * mrSges(7,3) + (t24 * t178 + t132) * t181 + (t22 * t113 + t106 * t173 + t110 * t175 + t108 * t177 - t89 * t29 / 0.2e1 + t30 * t171 + (t5 * t86 - t6 * t89) * mrSges(7,3)) * qJD(6) + (t79 * t75 + t123 + (m(6) * t79 - mrSges(6,3)) * t25 + ((-0.3e1 / 0.2e1 * Ifges(6,4) + t159 / 0.2e1 - t158 / 0.2e1) * t90 + (-0.3e1 / 0.2e1 * Ifges(6,1) + Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(6,2)) * t87) * qJD(2) - t95) * qJD(5) + t185 * t55) * t90 + (t1 * t35 + 0.2e1 * t132 * t90 + t190 * t6 + t191 * t5 + t2 * t34) * m(7) + (t48 * t126 - t36 * t55 + t101 - t157) * m(5) + ((-t48 * t84 + t49 * t82) * pkin(2) + t43 * t55 - t44 * t58) * m(4); ((-t99 - t128) * qJD(5) + t183) * t87 + ((-t127 - t147) * qJD(5) + t94) * t90; m(5) * t48 - t92 * mrSges(5,3) + (t99 * qJD(5) - t183) * t90 + (-t147 * qJD(5) + t94) * t87 + (-t69 + (-m(5) - m(6)) * t37 + t184) * qJD(2); t110 * t180 + t108 * t179 + t13 * t170 + t86 * t14 / 0.2e1 - t24 * t75 - t9 * t53 - t10 * t52 - pkin(5) * t21 - t7 * mrSges(6,2) + (-mrSges(6,1) - t113) * t8 + t147 * t25 + t116 * mrSges(7,3) + t193 * qJD(6) + ((t123 + qJD(5) * t106 / 0.2e1 + t25 * mrSges(6,3) + Ifges(6,4) * t141 / 0.2e1 + t95) * t90 + (t124 + (-t163 / 0.2e1 + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t90) * qJD(2) + t189) * t87) * qJD(2) + (-pkin(5) * t8 - t10 * t6 - t22 * t25 - t5 * t9) * m(7) + (m(7) * t116 + t192) * pkin(9); t78 - t22 * (mrSges(7,1) * t68 + mrSges(7,2) * t67) + (Ifges(7,1) * t67 - t162) * t175 + t29 * t174 + (Ifges(7,5) * t67 - Ifges(7,6) * t68) * t173 - t5 * t52 + t6 * t53 + (t5 * t67 + t6 * t68) * mrSges(7,3) + (-Ifges(7,2) * t68 + t30 + t65) * t177 + t186;];
tauc  = t11(:);
