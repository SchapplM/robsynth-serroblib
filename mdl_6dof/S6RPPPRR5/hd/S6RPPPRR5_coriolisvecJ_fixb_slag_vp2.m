% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:23
% EndTime: 2019-03-09 01:37:26
% DurationCPUTime: 2.15s
% Computational Cost: add. (2446->306), mult. (4458->425), div. (0->0), fcn. (2209->6), ass. (0->142)
t88 = cos(qJ(5));
t125 = t88 * qJD(1);
t78 = Ifges(6,4) * t125;
t183 = -t78 / 0.2e1;
t182 = qJD(5) / 0.2e1;
t86 = sin(qJ(5));
t134 = qJD(1) * t86;
t118 = mrSges(6,3) * t134;
t87 = cos(qJ(6));
t129 = qJD(5) * t87;
t85 = sin(qJ(6));
t60 = -t85 * t134 + t129;
t131 = qJD(5) * t85;
t61 = t87 * t134 + t131;
t138 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t60 + mrSges(7,2) * t61 + t118;
t79 = qJD(1) * qJ(2) + qJD(3);
t71 = qJD(1) * pkin(3) + t79;
t84 = -pkin(1) - qJ(3);
t72 = t84 * qJD(1) + qJD(2);
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t37 = t81 * t71 + t82 * t72;
t32 = qJD(1) * pkin(7) + t37;
t22 = qJD(4) * t88 - t32 * t86;
t16 = -qJD(5) * pkin(5) - t22;
t181 = -m(7) * t16 - t138;
t180 = t134 / 0.2e1;
t116 = Ifges(6,5) * t182;
t119 = t86 * t129;
t141 = t85 * t88;
t48 = -t81 * t141 - t82 * t87;
t140 = t87 * t88;
t51 = t82 * t140 + t81 * t85;
t179 = t51 * qJD(1) + qJD(6) * t48 - t119 * t81;
t50 = t81 * t140 - t82 * t85;
t97 = t82 * t141 - t81 * t87;
t178 = -t50 * qJD(1) - qJD(6) * t97 - t119 * t82;
t130 = qJD(5) * t86;
t120 = t85 * t130;
t177 = -t97 * qJD(1) - qJD(6) * t50 + t120 * t81;
t176 = -t48 * qJD(1) - qJD(6) * t51 + t120 * t82;
t36 = t71 * t82 - t81 * t72;
t31 = -qJD(1) * pkin(4) - t36;
t155 = Ifges(7,4) * t87;
t103 = -Ifges(7,2) * t85 + t155;
t156 = Ifges(7,4) * t85;
t105 = Ifges(7,1) * t87 - t156;
t106 = mrSges(7,1) * t85 + mrSges(7,2) * t87;
t23 = qJD(4) * t86 + t32 * t88;
t17 = qJD(5) * pkin(8) + t23;
t98 = -pkin(5) * t88 - pkin(8) * t86 - pkin(4);
t18 = qJD(1) * t98 - t36;
t5 = -t17 * t85 + t18 * t87;
t6 = t17 * t87 + t18 * t85;
t108 = t5 * t87 + t6 * t85;
t153 = Ifges(7,6) * t85;
t154 = Ifges(7,5) * t87;
t161 = t87 / 0.2e1;
t162 = -t85 / 0.2e1;
t164 = t61 / 0.2e1;
t157 = Ifges(7,4) * t61;
t74 = qJD(6) - t125;
t20 = Ifges(7,2) * t60 + Ifges(7,6) * t74 + t157;
t57 = Ifges(7,4) * t60;
t21 = Ifges(7,1) * t61 + Ifges(7,5) * t74 + t57;
t90 = -t108 * mrSges(7,3) + t21 * t161 + t16 * t106 + t74 * (-t153 + t154) / 0.2e1 + t60 * t103 / 0.2e1 + t105 * t164 + t20 * t162;
t175 = -t31 * mrSges(6,2) + t22 * mrSges(6,3) - Ifges(6,1) * t180 - t116 + t183 - t90;
t121 = mrSges(6,3) * t125;
t70 = -qJD(5) * mrSges(6,2) + t121;
t174 = m(6) * (-t22 * t86 + t23 * t88) + m(5) * t37 + t88 * t70 - t181 * t86;
t115 = -Ifges(6,6) * qJD(5) / 0.2e1;
t111 = pkin(5) * t86 - pkin(8) * t88;
t67 = qJD(2) * t82 + qJD(3) * t81;
t43 = t111 * qJD(5) - t67;
t33 = t43 * qJD(1);
t133 = qJD(5) * t22;
t68 = qJD(2) * t81 - qJD(3) * t82;
t59 = t68 * qJD(1);
t9 = t59 * t88 + t133;
t1 = qJD(6) * t5 + t33 * t85 + t87 * t9;
t2 = -qJD(6) * t6 + t33 * t87 - t85 * t9;
t109 = t1 * t87 - t2 * t85;
t123 = qJD(5) * qJD(6);
t127 = qJD(6) * t85;
t128 = qJD(5) * t88;
t38 = t87 * t123 + (-t86 * t127 + t87 * t128) * qJD(1);
t126 = qJD(6) * t87;
t39 = -t85 * t123 + (-t86 * t126 - t85 * t128) * qJD(1);
t172 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t38 - Ifges(7,6) * t39;
t62 = (-mrSges(6,1) * t88 + mrSges(6,2) * t86) * qJD(1);
t170 = -m(5) * t36 + m(6) * t31 + t62;
t168 = t38 / 0.2e1;
t167 = t39 / 0.2e1;
t166 = -t60 / 0.2e1;
t165 = -t61 / 0.2e1;
t163 = -t74 / 0.2e1;
t132 = qJD(5) * t23;
t10 = t59 * t86 + t132;
t83 = qJ(2) + pkin(3);
t137 = t81 * t83 + t82 * t84;
t56 = pkin(7) + t137;
t152 = t10 * t56;
t151 = t10 * t86;
t58 = t67 * qJD(1);
t146 = t58 * t81;
t145 = t58 * t82;
t124 = qJD(1) * qJD(5);
t117 = t86 * t124;
t114 = -t81 * t84 + t82 * t83;
t113 = m(3) * qJ(2) + mrSges(4,1) + mrSges(3,3);
t110 = -qJD(6) * t56 * t88 + t43;
t107 = mrSges(7,1) * t87 - mrSges(7,2) * t85;
t104 = Ifges(7,1) * t85 + t155;
t102 = Ifges(7,2) * t87 + t156;
t101 = Ifges(7,5) * t85 + Ifges(7,6) * t87;
t24 = mrSges(7,1) * t117 - mrSges(7,3) * t38;
t25 = -mrSges(7,2) * t117 + mrSges(7,3) * t39;
t100 = -t85 * t24 + t87 * t25;
t40 = -mrSges(7,2) * t74 + mrSges(7,3) * t60;
t41 = mrSges(7,1) * t74 - mrSges(7,3) * t61;
t99 = -t85 * t40 - t87 * t41;
t96 = t16 * t128 + t151;
t42 = -t114 + t98;
t95 = -qJD(6) * t42 + t56 * t130 - t68 * t88;
t94 = t151 + t88 * t9 + (-t22 * t88 - t23 * t86) * qJD(5);
t13 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t89 = qJD(1) ^ 2;
t93 = -t89 * mrSges(5,1) + t86 * t13 + (t138 * t88 - t70 * t86) * qJD(5);
t92 = t31 * mrSges(6,1) + t5 * mrSges(7,1) + t74 * Ifges(7,3) + t61 * Ifges(7,5) + t60 * Ifges(7,6) + t115 - (Ifges(6,4) * t86 + t88 * Ifges(6,2)) * qJD(1) / 0.2e1 - t6 * mrSges(7,2);
t73 = Ifges(7,3) * t117;
t66 = t111 * qJD(1);
t55 = -pkin(4) - t114;
t54 = (mrSges(6,1) * t86 + mrSges(6,2) * t88) * t124;
t15 = t56 * t140 + t42 * t85;
t14 = -t56 * t141 + t42 * t87;
t12 = t22 * t87 + t66 * t85;
t11 = -t22 * t85 + t66 * t87;
t8 = t38 * Ifges(7,1) + t39 * Ifges(7,4) + Ifges(7,5) * t117;
t7 = t38 * Ifges(7,4) + t39 * Ifges(7,2) + Ifges(7,6) * t117;
t4 = t110 * t87 + t85 * t95;
t3 = t110 * t85 - t87 * t95;
t19 = [t58 * mrSges(5,1) - t59 * mrSges(5,2) + t14 * t24 + t15 * t25 + t3 * t40 + t4 * t41 + t55 * t54 - t67 * t62 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + 0.2e1 * qJD(3) * mrSges(4,3) + 0.2e1 * qJD(2) * t113) * qJD(1) + m(5) * (t58 * t114 + t59 * t137 + t36 * t67 + t37 * t68) + m(4) * (qJD(2) * t79 - qJD(3) * t72 + (qJ(2) * qJD(2) - qJD(3) * t84) * qJD(1)) + m(6) * (-t31 * t67 - t55 * t58) + m(7) * (t1 * t15 + t14 * t2 + t3 * t6 + t4 * t5) + (t68 * t70 + t9 * mrSges(6,3) + t58 * mrSges(6,1) - t73 / 0.2e1 + m(6) * (t23 * t68 + t56 * t9) + (0.3e1 / 0.2e1 * t78 + t116 + (-m(6) * t22 - t181) * t56 - t175) * qJD(5) + t172) * t88 + (-t58 * mrSges(6,2) + t8 * t161 + t7 * t162 + t56 * t13 + t103 * t167 + t105 * t168 + t138 * t68 + (mrSges(6,3) + t106) * t10 + (-t1 * t85 - t2 * t87) * mrSges(7,3) + m(6) * (-t22 * t68 + t152) + m(7) * (t16 * t68 + t152) + (t101 * t163 + t102 * t166 + t104 * t165 + t16 * t107 - t87 * t20 / 0.2e1 + t21 * t162 + (t5 * t85 - t6 * t87) * mrSges(7,3)) * qJD(6) + (-t56 * t70 + t115 + (-m(6) * t56 - mrSges(6,3)) * t23 + ((-0.3e1 / 0.2e1 * Ifges(6,4) + t154 / 0.2e1 - t153 / 0.2e1) * t86 + (-0.3e1 / 0.2e1 * Ifges(6,2) - Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(6,1)) * t88) * qJD(1) + t92) * qJD(5)) * t86; -t97 * t24 + t51 * t25 + t81 * t54 + t176 * t41 + t178 * t40 + (mrSges(5,2) * t81 - t113) * t89 + t93 * t82 + m(6) * (t82 * t94 - t146) + m(5) * (t59 * t82 - t146) + ((-t79 - qJD(3)) * m(4) + t170 * t82 - t174 * t81) * qJD(1) + (t1 * t51 + t176 * t5 + t178 * t6 - t2 * t97 + t82 * t96) * m(7); -t89 * mrSges(4,3) + t48 * t24 + t50 * t25 + (-t89 * mrSges(5,2) - t54) * t82 + t177 * t41 + t179 * t40 + t93 * t81 + m(6) * (t81 * t94 + t145) + m(5) * (t59 * t81 + t145) + ((t72 + qJD(2)) * m(4) + t170 * t81 + t174 * t82) * qJD(1) + (t1 * t50 + t177 * t5 + t179 * t6 + t2 * t48 + t81 * t96) * m(7); (-t13 + (t40 * t87 - t41 * t85 - t121 + t70) * qJD(5) + m(6) * (-t10 + t132) + m(7) * (t6 * t129 - t5 * t131 - t10)) * t88 + (t99 * qJD(6) + (-t118 + t138) * qJD(5) + m(6) * (t9 - t133) + m(7) * (qJD(5) * t16 - t5 * t126 - t6 * t127 + t109) + t100) * t86; t102 * t167 + t104 * t168 + t7 * t161 + t85 * t8 / 0.2e1 - t22 * t70 - t12 * t40 - t11 * t41 - t9 * mrSges(6,2) - pkin(5) * t13 - t138 * t23 + (-mrSges(6,1) - t107) * t10 + t109 * mrSges(7,3) + t90 * qJD(6) + ((t23 * mrSges(6,3) + Ifges(6,4) * t180 + t101 * t182 + t115 - t92) * t86 + (t183 + t116 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t134 + t175) * t88) * qJD(1) + (-pkin(5) * t10 - t11 * t5 - t12 * t6 - t16 * t23) * m(7) + (t100 + m(7) * t109 + (-m(7) * t108 + t99) * qJD(6)) * pkin(8); t73 - t16 * (mrSges(7,1) * t61 + mrSges(7,2) * t60) + (Ifges(7,1) * t60 - t157) * t165 + t20 * t164 + (Ifges(7,5) * t60 - Ifges(7,6) * t61) * t163 - t5 * t40 + t6 * t41 + (t5 * t60 + t6 * t61) * mrSges(7,3) + (-Ifges(7,2) * t61 + t21 + t57) * t166 - t172;];
tauc  = t19(:);
