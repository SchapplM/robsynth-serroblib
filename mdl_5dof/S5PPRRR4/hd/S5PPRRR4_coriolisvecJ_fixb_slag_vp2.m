% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:37
% EndTime: 2019-12-05 15:18:45
% DurationCPUTime: 2.12s
% Computational Cost: add. (2195->284), mult. (6151->437), div. (0->0), fcn. (4967->12), ass. (0->152)
t93 = cos(qJ(4));
t132 = qJD(3) * t93;
t82 = Ifges(5,4) * t132;
t179 = -t82 / 0.2e1;
t178 = qJD(4) / 0.2e1;
t90 = sin(qJ(4));
t134 = qJD(3) * t90;
t177 = t134 / 0.2e1;
t116 = Ifges(5,5) * t178;
t86 = cos(pkin(11));
t87 = cos(pkin(6));
t141 = t86 * t87;
t83 = sin(pkin(11));
t85 = sin(pkin(5));
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t100 = (t141 * t91 + t83 * t94) * t85;
t84 = sin(pkin(6));
t143 = t84 * t91;
t88 = cos(pkin(5));
t78 = qJD(1) * t88 + qJD(2);
t37 = qJD(1) * t100 + t143 * t78;
t35 = qJD(3) * pkin(8) + t37;
t135 = qJD(1) * t85;
t120 = t86 * t135;
t58 = -t120 * t84 + t78 * t87;
t21 = -t35 * t90 + t58 * t93;
t144 = t83 * t91;
t36 = (t120 * t87 + t78 * t84) * t94 - t135 * t144;
t34 = -qJD(3) * pkin(3) - t36;
t92 = cos(qJ(5));
t153 = Ifges(6,4) * t92;
t89 = sin(qJ(5));
t107 = -Ifges(6,2) * t89 + t153;
t154 = Ifges(6,4) * t89;
t109 = Ifges(6,1) * t92 - t154;
t110 = mrSges(6,1) * t89 + mrSges(6,2) * t92;
t22 = t35 * t93 + t58 * t90;
t20 = qJD(4) * pkin(9) + t22;
t75 = -pkin(4) * t93 - pkin(9) * t90 - pkin(3);
t29 = qJD(3) * t75 - t36;
t5 = -t20 * t89 + t29 * t92;
t6 = t20 * t92 + t29 * t89;
t112 = t5 * t92 + t6 * t89;
t151 = Ifges(6,6) * t89;
t152 = Ifges(6,5) * t92;
t160 = t92 / 0.2e1;
t161 = -t89 / 0.2e1;
t131 = qJD(4) * t89;
t70 = t134 * t92 + t131;
t163 = t70 / 0.2e1;
t19 = -qJD(4) * pkin(4) - t21;
t155 = Ifges(6,4) * t70;
t130 = qJD(4) * t92;
t69 = -t134 * t89 + t130;
t80 = qJD(5) - t132;
t39 = Ifges(6,2) * t69 + Ifges(6,6) * t80 + t155;
t68 = Ifges(6,4) * t69;
t40 = Ifges(6,1) * t70 + Ifges(6,5) * t80 + t68;
t96 = -t112 * mrSges(6,3) + t40 * t160 + t39 * t161 + t80 * (-t151 + t152) / 0.2e1 + t69 * t107 / 0.2e1 + t109 * t163 + t19 * t110;
t176 = -t34 * mrSges(5,2) + t21 * mrSges(5,3) - Ifges(5,1) * t177 - t116 + t179 - t96;
t101 = (t141 * t94 - t144) * t85;
t142 = t84 * t94;
t175 = t88 * t142 + t101;
t32 = (qJD(1) * t101 + t142 * t78) * qJD(3);
t174 = t36 * qJD(3) - t32;
t113 = pkin(4) * t90 - pkin(9) * t93;
t73 = t113 * qJD(4);
t25 = (t73 + t37) * qJD(3);
t7 = qJD(4) * t21 + t32 * t93;
t1 = qJD(5) * t5 + t25 * t89 + t7 * t92;
t2 = -qJD(5) * t6 + t25 * t92 - t7 * t89;
t173 = t1 * t92 - t2 * t89;
t123 = qJD(4) * qJD(5);
t128 = qJD(5) * t89;
t129 = qJD(4) * t93;
t54 = t92 * t123 + (-t128 * t90 + t129 * t92) * qJD(3);
t127 = qJD(5) * t92;
t55 = -t89 * t123 + (-t127 * t90 - t129 * t89) * qJD(3);
t172 = -t2 * mrSges(6,1) + t1 * mrSges(6,2) - Ifges(6,5) * t54 - Ifges(6,6) * t55;
t170 = 2 * m(5);
t169 = pkin(8) / 0.2e1;
t168 = -t36 / 0.2e1;
t167 = t54 / 0.2e1;
t166 = t55 / 0.2e1;
t165 = -t69 / 0.2e1;
t164 = -t70 / 0.2e1;
t162 = -t80 / 0.2e1;
t46 = t143 * t88 + t100;
t62 = -t84 * t85 * t86 + t87 * t88;
t26 = t46 * t90 - t62 * t93;
t8 = qJD(4) * t22 + t32 * t90;
t157 = t26 * t8;
t63 = t143 * t90 - t93 * t87;
t156 = t63 * t8;
t33 = t37 * qJD(3);
t149 = t33 * t175;
t140 = t89 * t93;
t139 = t92 * t93;
t138 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t69 - mrSges(6,2) * t70 - mrSges(5,3) * t134;
t136 = Ifges(5,6) * qJD(4);
t133 = qJD(3) * t91;
t126 = qJD(5) * t93;
t124 = qJD(4) * qJD(3);
t121 = t8 * t169;
t119 = t84 * t133;
t118 = qJD(3) * t142;
t117 = t90 * t124;
t115 = -t136 / 0.2e1;
t111 = mrSges(6,1) * t92 - mrSges(6,2) * t89;
t108 = Ifges(6,1) * t89 + t153;
t106 = Ifges(6,2) * t92 + t154;
t105 = Ifges(6,5) * t89 + Ifges(6,6) * t92;
t27 = t46 * t93 + t62 * t90;
t12 = -t175 * t89 + t27 * t92;
t11 = -t175 * t92 - t27 * t89;
t64 = t143 * t93 + t87 * t90;
t50 = -t142 * t92 - t64 * t89;
t104 = t142 * t89 - t64 * t92;
t98 = t22 * mrSges(5,3) + t6 * mrSges(6,2) - t80 * Ifges(6,3) - t70 * Ifges(6,5) - t69 * Ifges(6,6) + t136 / 0.2e1 + (Ifges(5,4) * t90 + t93 * Ifges(5,2)) * qJD(3) / 0.2e1 - t34 * mrSges(5,1) - t5 * mrSges(6,1);
t95 = qJD(3) ^ 2;
t79 = Ifges(6,3) * t117;
t77 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t132;
t72 = t113 * qJD(3);
t71 = (-mrSges(5,1) * t93 + mrSges(5,2) * t90) * qJD(3);
t67 = (mrSges(5,1) * t90 + mrSges(5,2) * t93) * t124;
t60 = pkin(8) * t139 + t75 * t89;
t59 = -pkin(8) * t140 + t75 * t92;
t57 = mrSges(6,1) * t80 - mrSges(6,3) * t70;
t56 = -mrSges(6,2) * t80 + mrSges(6,3) * t69;
t49 = qJD(4) * t64 + t118 * t90;
t48 = -qJD(4) * t63 + t118 * t93;
t44 = -mrSges(6,2) * t117 + mrSges(6,3) * t55;
t43 = mrSges(6,1) * t117 - mrSges(6,3) * t54;
t42 = t46 * qJD(3);
t41 = t175 * qJD(3);
t31 = -t75 * t128 + t73 * t92 + (-t126 * t92 + t131 * t90) * pkin(8);
t30 = t75 * t127 + t73 * t89 + (-t126 * t89 - t130 * t90) * pkin(8);
t28 = -mrSges(6,1) * t55 + mrSges(6,2) * t54;
t24 = t54 * Ifges(6,1) + t55 * Ifges(6,4) + Ifges(6,5) * t117;
t23 = t54 * Ifges(6,4) + t55 * Ifges(6,2) + Ifges(6,6) * t117;
t18 = qJD(5) * t104 + t119 * t92 - t48 * t89;
t17 = qJD(5) * t50 + t119 * t89 + t48 * t92;
t16 = t139 * t36 + t37 * t89;
t15 = -t140 * t36 + t37 * t92;
t14 = t21 * t92 + t72 * t89;
t13 = -t21 * t89 + t72 * t92;
t10 = -qJD(4) * t26 + t41 * t93;
t9 = qJD(4) * t27 + t41 * t90;
t4 = qJD(5) * t11 + t10 * t92 + t42 * t89;
t3 = -qJD(5) * t12 - t10 * t89 + t42 * t92;
t38 = [t10 * t77 + t11 * t43 + t12 * t44 + t26 * t28 + t3 * t57 + t4 * t56 + t42 * t71 - t175 * t67 - t138 * t9 + m(4) * (t32 * t46 - t36 * t42 + t37 * t41 - t149) + m(5) * (t10 * t22 - t21 * t9 + t27 * t7 + t34 * t42 - t149 + t157) + m(6) * (t1 * t12 + t11 * t2 + t19 * t9 + t3 * t5 + t4 * t6 + t157) + (-t42 * mrSges(4,1) - t41 * mrSges(4,2) + (t26 * t93 - t27 * t90) * qJD(4) * mrSges(5,3)) * qJD(3); t17 * t56 + t18 * t57 + t63 * t28 + t50 * t43 - t104 * t44 + t48 * t77 - t138 * t49 + (t63 * t93 - t64 * t90) * mrSges(5,3) * t124 + m(5) * (-t21 * t49 + t22 * t48 + t64 * t7 + t156) + m(6) * (-t1 * t104 + t17 * t6 + t18 * t5 + t19 * t49 + t2 * t50 + t156) + ((-t95 * mrSges(4,2) - t67) * t94 + m(5) * (t133 * t34 - t33 * t94) + (-m(4) * t174 - t95 * mrSges(4,1) + qJD(3) * t71) * t91) * t84; -pkin(3) * t67 - t37 * t71 + t59 * t43 + t60 * t44 + (t31 - t15) * t57 + (t30 - t16) * t56 + t174 * mrSges(4,2) - m(6) * (t15 * t5 + t16 * t6) + m(6) * (t1 * t60 + t2 * t59 + t30 * t6 + t31 * t5) + (-t34 * t37 / 0.2e1 - pkin(3) * t33 / 0.2e1) * t170 + (-t36 * t77 + t7 * mrSges(5,3) - t79 / 0.2e1 - t33 * mrSges(5,1) + (t168 * t22 + t169 * t7) * t170 + (0.3e1 / 0.2e1 * t82 + t116 + (-m(5) * t21 + m(6) * t19 - t138) * pkin(8) - t176) * qJD(4) + t172) * t93 + (t24 * t160 + t23 * t161 + pkin(8) * t28 + t109 * t167 + t107 * t166 + t33 * mrSges(5,2) + (mrSges(5,3) + t110) * t8 + t138 * t36 + (-t1 * t89 - t2 * t92) * mrSges(6,3) + 0.2e1 * (t19 * t168 + t121) * m(6) + (t121 + t21 * t36 / 0.2e1) * t170 + (t105 * t162 + t106 * t165 + t108 * t164 + t19 * t111 - t92 * t39 / 0.2e1 + t40 * t161 + (t5 * t89 - t6 * t92) * mrSges(6,3)) * qJD(5) + (t115 + (-m(5) * t22 - t77) * pkin(8) + ((t152 / 0.2e1 - t151 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t90 + (-Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t93) * qJD(3) - t98) * qJD(4)) * t90; t108 * t167 + t106 * t166 + t23 * t160 + t89 * t24 / 0.2e1 - t21 * t77 - t14 * t56 - t13 * t57 - pkin(4) * t28 - t7 * mrSges(5,2) + (-mrSges(5,1) - t111) * t8 + t138 * t22 + t173 * mrSges(6,3) + t96 * qJD(5) + ((Ifges(5,4) * t177 + t105 * t178 + t115 + t98) * t90 + (t116 + t179 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t134 + t176) * t93) * qJD(3) + (-pkin(4) * t8 - t13 * t5 - t14 * t6 - t19 * t22) * m(6) + (-t89 * t43 + t92 * t44 + m(6) * t173 + (-m(6) * t112 - t89 * t56 - t92 * t57) * qJD(5)) * pkin(9); t79 - t19 * (mrSges(6,1) * t70 + mrSges(6,2) * t69) + (Ifges(6,1) * t69 - t155) * t164 + t39 * t163 + (Ifges(6,5) * t69 - Ifges(6,6) * t70) * t162 - t5 * t56 + t6 * t57 + (t5 * t69 + t6 * t70) * mrSges(6,3) + (-Ifges(6,2) * t70 + t40 + t68) * t165 - t172;];
tauc = t38(:);
