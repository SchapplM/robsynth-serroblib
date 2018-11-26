% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:54:02
% EndTime: 2018-11-23 14:54:04
% DurationCPUTime: 2.39s
% Computational Cost: add. (2286->299), mult. (5500->415), div. (0->0), fcn. (3722->10), ass. (0->156)
t83 = sin(pkin(6));
t143 = qJD(1) * t83;
t91 = cos(qJ(2));
t128 = t91 * t143;
t84 = cos(pkin(11));
t120 = t84 * t128;
t88 = sin(qJ(2));
t129 = t88 * t143;
t82 = sin(pkin(11));
t73 = t82 * t129;
t57 = -t73 + t120;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t117 = pkin(5) * t90 + pkin(9) * t87;
t66 = t117 * qJD(5) + qJD(4);
t192 = -t57 + t66;
t89 = cos(qJ(6));
t159 = Ifges(7,4) * t89;
t86 = sin(qJ(6));
t109 = -Ifges(7,2) * t86 + t159;
t160 = Ifges(7,4) * t86;
t111 = Ifges(7,1) * t89 - t160;
t112 = mrSges(7,1) * t86 + mrSges(7,2) * t89;
t72 = qJD(2) * pkin(2) + t128;
t43 = t72 * t84 - t73;
t114 = qJD(4) - t43;
t33 = (-pkin(3) - pkin(8)) * qJD(2) + t114;
t85 = cos(pkin(6));
t77 = qJD(1) * t85 + qJD(3);
t25 = t33 * t87 + t77 * t90;
t23 = qJD(5) * pkin(9) + t25;
t100 = pkin(5) * t87 - pkin(9) * t90 + qJ(4);
t44 = t84 * t129 + t72 * t82;
t27 = t100 * qJD(2) + t44;
t5 = -t23 * t86 + t27 * t89;
t6 = t23 * t89 + t27 * t86;
t115 = t5 * t89 + t6 * t86;
t157 = Ifges(7,6) * t86;
t158 = Ifges(7,5) * t89;
t169 = t89 / 0.2e1;
t170 = -t86 / 0.2e1;
t139 = qJD(5) * t86;
t140 = qJD(2) * t90;
t68 = t89 * t140 + t139;
t173 = t68 / 0.2e1;
t24 = t33 * t90 - t77 * t87;
t22 = -qJD(5) * pkin(5) - t24;
t161 = Ifges(7,4) * t68;
t138 = qJD(5) * t89;
t67 = -t86 * t140 + t138;
t141 = qJD(2) * t87;
t80 = qJD(6) + t141;
t29 = Ifges(7,2) * t67 + Ifges(7,6) * t80 + t161;
t65 = Ifges(7,4) * t67;
t30 = Ifges(7,1) * t68 + Ifges(7,5) * t80 + t65;
t191 = -t115 * mrSges(7,3) + (-t157 + t158) * t80 / 0.2e1 + t109 * t67 / 0.2e1 + t111 * t173 + t22 * t112 + t30 * t169 + t29 * t170;
t52 = -mrSges(7,2) * t80 + mrSges(7,3) * t67;
t53 = mrSges(7,1) * t80 - mrSges(7,3) * t68;
t183 = -m(7) * t115 - t86 * t52 - t89 * t53;
t133 = qJD(5) * qJD(2);
t124 = t90 * t133;
t132 = qJD(5) * qJD(6);
t137 = qJD(6) * t90;
t50 = t89 * t132 + (-t86 * t137 - t87 * t138) * qJD(2);
t31 = mrSges(7,1) * t124 - mrSges(7,3) * t50;
t51 = -t86 * t132 + (-t89 * t137 + t87 * t139) * qJD(2);
t32 = -mrSges(7,2) * t124 + mrSges(7,3) * t51;
t190 = t183 * qJD(6) - t86 * t31 + t89 * t32;
t125 = -pkin(2) * t84 - pkin(3);
t79 = -pkin(8) + t125;
t127 = qJD(5) * t79 * t90;
t148 = t87 * t89;
t149 = t86 * t87;
t168 = pkin(2) * t82;
t64 = t100 + t168;
t34 = -t79 * t149 + t64 * t89;
t60 = (t82 * t91 + t84 * t88) * t83;
t55 = qJD(1) * t60;
t189 = t34 * qJD(6) + t89 * t127 - t55 * t148 + t192 * t86;
t35 = t79 * t148 + t64 * t86;
t188 = -t35 * qJD(6) - t86 * t127 + t55 * t149 + t192 * t89;
t145 = Ifges(6,5) * qJD(5);
t162 = Ifges(6,4) * t87;
t37 = qJD(2) * qJ(4) + t44;
t187 = t37 * mrSges(6,2) + t145 / 0.2e1 + (t90 * Ifges(6,1) - t162) * qJD(2) / 0.2e1 - t24 * mrSges(6,3) + t191;
t122 = -Ifges(6,6) * qJD(5) / 0.2e1;
t70 = qJD(2) * t73;
t26 = -t70 + (t66 + t120) * qJD(2);
t136 = t24 * qJD(5);
t56 = qJD(2) * t60;
t48 = qJD(1) * t56;
t7 = t48 * t87 + t136;
t1 = t5 * qJD(6) + t26 * t86 + t7 * t89;
t2 = -t6 * qJD(6) + t26 * t89 - t7 * t86;
t116 = t1 * t89 - t2 * t86;
t185 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t50 + Ifges(7,6) * t51;
t130 = mrSges(6,3) * t140;
t146 = qJD(5) * mrSges(6,1) + mrSges(7,1) * t67 - mrSges(7,2) * t68 - t130;
t184 = m(7) * t22 - t146;
t135 = t25 * qJD(5);
t21 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t8 = -t90 * t48 + t135;
t182 = -m(6) * (-t8 + t135) + m(7) * (-t6 * t138 + t5 * t139 + t8) + t21;
t180 = 0.2e1 * m(6);
t179 = t50 / 0.2e1;
t178 = t51 / 0.2e1;
t177 = -t55 / 0.2e1;
t176 = -t67 / 0.2e1;
t174 = -t68 / 0.2e1;
t172 = -t80 / 0.2e1;
t150 = t84 * t91;
t151 = t82 * t88;
t59 = (-t150 + t151) * t83;
t102 = t90 * t59 - t85 * t87;
t165 = t102 * t8;
t156 = t37 * t57;
t155 = t48 * t59;
t147 = mrSges(4,1) - mrSges(5,2);
t142 = qJD(2) * t83;
t134 = -t57 + qJD(4);
t131 = -t79 * t8 / 0.2e1;
t126 = mrSges(6,3) * t141;
t123 = -t145 / 0.2e1;
t119 = t142 * t150;
t113 = mrSges(7,1) * t89 - mrSges(7,2) * t86;
t110 = Ifges(7,1) * t86 + t159;
t108 = Ifges(7,2) * t89 + t160;
t106 = Ifges(7,5) * t86 + Ifges(7,6) * t89;
t41 = -t70 + (qJD(4) + t120) * qJD(2);
t58 = -t142 * t151 + t119;
t104 = t37 * t58 + t41 * t60;
t40 = t59 * t87 + t85 * t90;
t18 = t40 * t89 + t60 * t86;
t17 = -t40 * t86 + t60 * t89;
t81 = qJ(4) + t168;
t101 = t37 * qJD(4) + t41 * t81;
t75 = -qJD(5) * mrSges(6,2) - t126;
t99 = t89 * t52 - t86 * t53 + t75;
t95 = t37 * mrSges(6,1) + t5 * mrSges(7,1) + t80 * Ifges(7,3) + t68 * Ifges(7,5) + t67 * Ifges(7,6) + t122 - (Ifges(6,4) * t90 - t87 * Ifges(6,2)) * qJD(2) / 0.2e1 - t6 * mrSges(7,2);
t94 = m(6) * (t7 - t136) + m(7) * (qJD(5) * t22 + t116) + t190;
t92 = qJD(2) ^ 2;
t78 = Ifges(7,3) * t124;
t71 = t117 * qJD(2);
t69 = (mrSges(6,1) * t87 + mrSges(6,2) * t90) * qJD(2);
t63 = (mrSges(6,1) * t90 - mrSges(6,2) * t87) * t133;
t49 = qJD(1) * t119 - t70;
t36 = -qJD(2) * pkin(3) + t114;
t16 = t40 * qJD(5) - t90 * t56;
t15 = t102 * qJD(5) + t56 * t87;
t14 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + Ifges(7,5) * t124;
t13 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + Ifges(7,6) * t124;
t10 = t24 * t89 + t71 * t86;
t9 = -t24 * t86 + t71 * t89;
t4 = -t18 * qJD(6) - t15 * t86 + t58 * t89;
t3 = t17 * qJD(6) + t15 * t89 + t58 * t86;
t11 = [t15 * t75 + t17 * t31 + t18 * t32 - t102 * t21 + t3 * t52 + t4 * t53 + t58 * t69 + t60 * t63 + (-mrSges(3,1) * t88 - mrSges(3,2) * t91) * t92 * t83 - t146 * t16 + m(7) * (t1 * t18 + t16 * t22 + t17 * t2 + t3 * t6 + t4 * t5 - t165) + m(4) * (-t43 * t56 + t44 * t58 + t49 * t60 + t155) + m(6) * (t15 * t25 - t16 * t24 + t40 * t7 + t104 - t165) + m(5) * (t36 * t56 + t104 + t155) + ((-mrSges(4,2) + mrSges(5,3)) * t58 - t147 * t56 + (t102 * t87 - t40 * t90) * qJD(5) * mrSges(6,3)) * qJD(2); -t49 * mrSges(4,2) + t41 * mrSges(5,3) + t34 * t31 + t35 * t32 + t81 * t63 + t134 * t69 + t188 * t53 + t189 * t52 - t147 * t48 + (t57 * mrSges(4,2) + t134 * mrSges(5,3) + t147 * t55) * qJD(2) + m(6) * t101 - m(6) * t156 + (t78 / 0.2e1 + t41 * mrSges(6,1) - t55 * t75 - t7 * mrSges(6,3) + (t7 * t79 / 0.2e1 + t25 * t177) * t180 + (0.3e1 / 0.2e1 * Ifges(6,4) * t141 + t123 + (-m(6) * t24 + t184) * t79 - t187) * qJD(5) + t185) * t87 + (t14 * t169 + t13 * t170 - t79 * t21 + t111 * t179 + t109 * t178 + t41 * mrSges(6,2) + (mrSges(6,3) + t112) * t8 + (-t1 * t86 - t2 * t89) * mrSges(7,3) + (t24 * t177 + t131) * t180 + (t106 * t172 + t110 * t174 + t108 * t176 + t22 * t113 - t89 * t29 / 0.2e1 + t30 * t170 + (t5 * t86 - t6 * t89) * mrSges(7,3)) * qJD(6) + (t79 * t75 + t122 + (m(6) * t79 - mrSges(6,3)) * t25 + ((-0.3e1 / 0.2e1 * Ifges(6,4) + t158 / 0.2e1 - t157 / 0.2e1) * t90 + (-0.3e1 / 0.2e1 * Ifges(6,1) + Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(6,2)) * t87) * qJD(2) + t95) * qJD(5) + t184 * t55) * t90 + (t1 * t35 + 0.2e1 * t131 * t90 + t188 * t5 + t189 * t6 + t2 * t34) * m(7) + (t48 * t125 - t36 * t55 + t101 - t156) * m(5) + ((-t48 * t84 + t49 * t82) * pkin(2) + t43 * t55 - t44 * t57) * m(4); ((-t99 - t126) * qJD(5) + t182) * t87 + ((-t130 - t146) * qJD(5) + t94) * t90; m(5) * t48 - t92 * mrSges(5,3) + (t99 * qJD(5) - t182) * t90 + (-t146 * qJD(5) + t94) * t87 + (-t69 + (-m(6) - m(5)) * t37 + t183) * qJD(2); t110 * t179 + t108 * t178 + t13 * t169 + t86 * t14 / 0.2e1 - t24 * t75 - t10 * t52 - t9 * t53 - pkin(5) * t21 - t7 * mrSges(6,2) + (-mrSges(6,1) - t113) * t8 + t146 * t25 + t116 * mrSges(7,3) + t191 * qJD(6) + ((t122 + qJD(5) * t106 / 0.2e1 + t25 * mrSges(6,3) + Ifges(6,4) * t140 / 0.2e1 - t95) * t90 + (t123 + (-t162 / 0.2e1 + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t90) * qJD(2) + t187) * t87) * qJD(2) + (-pkin(5) * t8 - t10 * t6 - t22 * t25 - t5 * t9) * m(7) + (m(7) * t116 + t190) * pkin(9); t78 - t22 * (mrSges(7,1) * t68 + mrSges(7,2) * t67) + (Ifges(7,1) * t67 - t161) * t174 + t29 * t173 + (Ifges(7,5) * t67 - Ifges(7,6) * t68) * t172 - t5 * t52 + t6 * t53 + (t5 * t67 + t6 * t68) * mrSges(7,3) + (-Ifges(7,2) * t68 + t30 + t65) * t176 + t185;];
tauc  = t11(:);
