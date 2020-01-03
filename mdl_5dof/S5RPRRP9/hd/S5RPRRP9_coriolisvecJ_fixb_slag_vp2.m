% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:36
% EndTime: 2019-12-31 18:48:44
% DurationCPUTime: 3.38s
% Computational Cost: add. (3311->298), mult. (9162->392), div. (0->0), fcn. (6781->6), ass. (0->139)
t189 = mrSges(5,1) + mrSges(6,1);
t182 = Ifges(5,1) + Ifges(6,1);
t181 = Ifges(6,4) + Ifges(5,5);
t110 = qJD(3) + qJD(4);
t113 = sin(qJ(4));
t111 = sin(pkin(8));
t114 = sin(qJ(3));
t112 = cos(pkin(8));
t115 = cos(qJ(3));
t141 = t112 * t115;
t128 = t111 * t114 - t141;
t124 = t128 * qJD(1);
t122 = t113 * t124;
t157 = cos(qJ(4));
t97 = t111 * t115 + t112 * t114;
t89 = t97 * qJD(1);
t119 = t157 * t89 - t122;
t80 = t157 * t124;
t65 = t113 * t89 + t80;
t152 = Ifges(6,5) * t65;
t61 = Ifges(5,4) * t65;
t179 = t181 * t110 + t182 * t119 + t152 - t61;
t188 = t179 / 0.2e1;
t180 = -Ifges(5,6) + Ifges(6,6);
t137 = qJD(4) * t113;
t146 = pkin(6) + qJ(2);
t101 = t146 * t111;
t98 = qJD(1) * t101;
t102 = t146 * t112;
t99 = qJD(1) * t102;
t71 = -t114 * t98 + t115 * t99;
t50 = -pkin(7) * t124 + t71;
t134 = t157 * t50;
t70 = -t114 * t99 - t115 * t98;
t49 = -pkin(7) * t89 + t70;
t15 = t113 * t49 + t134;
t187 = -pkin(3) * t137 + t15;
t186 = (m(3) * qJ(2) + mrSges(3,3)) * (t111 ^ 2 + t112 ^ 2);
t48 = qJD(3) * pkin(3) + t49;
t13 = t113 * t48 + t134;
t133 = -t112 * pkin(2) - pkin(1);
t78 = pkin(3) * t128 + t133;
t72 = qJD(1) * t78 + qJD(2);
t14 = t65 * pkin(4) - qJ(5) * t119 + t72;
t60 = Ifges(6,5) * t119;
t21 = t110 * Ifges(6,6) + t65 * Ifges(6,3) + t60;
t153 = Ifges(5,4) * t119;
t22 = -t65 * Ifges(5,2) + t110 * Ifges(5,6) + t153;
t185 = -t13 * mrSges(5,3) + t21 / 0.2e1 - t22 / 0.2e1 + t14 * mrSges(6,1) + t72 * mrSges(5,1);
t33 = pkin(4) * t119 + qJ(5) * t65;
t142 = t113 * t50;
t12 = t157 * t48 - t142;
t175 = qJD(5) - t12;
t10 = -t110 * pkin(4) + t175;
t184 = t72 * mrSges(5,2) + t10 * mrSges(6,2) - t12 * mrSges(5,3) - t14 * mrSges(6,3);
t147 = -Ifges(5,4) + Ifges(6,5);
t149 = t65 * mrSges(5,3);
t51 = -mrSges(5,2) * t110 - t149;
t54 = -mrSges(6,2) * t65 + mrSges(6,3) * t110;
t145 = t51 + t54;
t155 = mrSges(5,3) * t119;
t178 = -mrSges(6,2) * t119 + t189 * t110 - t155;
t11 = t110 * qJ(5) + t13;
t176 = t11 * t119;
t174 = t180 * t119;
t170 = -t65 / 0.2e1;
t169 = t65 / 0.2e1;
t167 = -t119 / 0.2e1;
t166 = t119 / 0.2e1;
t90 = t128 * qJD(3);
t165 = -t90 / 0.2e1;
t91 = t97 * qJD(3);
t164 = -t91 / 0.2e1;
t163 = pkin(3) * t89;
t73 = -t115 * t101 - t102 * t114;
t58 = -pkin(7) * t97 + t73;
t74 = -t114 * t101 + t115 * t102;
t59 = -pkin(7) * t128 + t74;
t127 = -t113 * t59 + t157 * t58;
t125 = t97 * qJD(2);
t46 = -qJD(1) * t125 - qJD(3) * t71;
t83 = qJD(1) * t90;
t117 = pkin(7) * t83 + t46;
t105 = qJD(2) * t141;
t138 = qJD(3) * t115;
t139 = qJD(2) * t111;
t45 = qJD(1) * t105 - t98 * t138 + (-qJD(1) * t139 - qJD(3) * t99) * t114;
t84 = qJD(1) * t91;
t41 = -pkin(7) * t84 + t45;
t5 = t13 * qJD(4) + t113 * t41 - t117 * t157;
t162 = t127 * t5;
t69 = -t113 * t128 + t157 * t97;
t161 = t5 * t69;
t160 = t84 * pkin(3);
t159 = -t110 / 0.2e1;
t156 = mrSges(4,3) * t89;
t154 = Ifges(4,4) * t89;
t151 = pkin(3) * t113;
t31 = -t157 * t83 + (-qJD(4) * t89 - t84) * t113 - qJD(4) * t80;
t150 = t31 * mrSges(6,2);
t148 = mrSges(6,2) + mrSges(5,3);
t136 = t157 * pkin(3);
t132 = t84 * mrSges(4,1) - t83 * mrSges(4,2);
t131 = qJD(4) * t157;
t130 = pkin(3) * t131;
t20 = t113 * t58 + t157 * t59;
t126 = -t113 * t97 - t128 * t157;
t4 = t113 * t117 + t48 * t131 - t137 * t50 + t157 * t41;
t123 = mrSges(4,3) * t124;
t2 = qJD(5) * t110 + t4;
t32 = -qJD(4) * t122 - t113 * t83 + t131 * t89 + t157 * t84;
t120 = -t4 * mrSges(5,2) + t2 * mrSges(6,3) + t180 * t32 + t181 * t31 - t189 * t5;
t55 = t105 - t101 * t138 + (-qJD(3) * t102 - t139) * t114;
t56 = -qJD(3) * t74 - t125;
t118 = pkin(7) * t90 + t56;
t107 = -t136 - pkin(4);
t106 = qJ(5) + t151;
t104 = t130 + qJD(5);
t100 = qJD(1) * t133 + qJD(2);
t85 = Ifges(4,4) * t124;
t77 = qJD(3) * mrSges(4,1) - t156;
t76 = -qJD(3) * mrSges(4,2) - t123;
t63 = t89 * Ifges(4,1) + Ifges(4,5) * qJD(3) - t85;
t62 = -Ifges(4,2) * t124 + Ifges(4,6) * qJD(3) + t154;
t43 = -pkin(7) * t91 + t55;
t40 = qJD(4) * t69 - t113 * t90 + t157 * t91;
t39 = qJD(4) * t126 - t113 * t91 - t157 * t90;
t35 = mrSges(5,1) * t65 + mrSges(5,2) * t119;
t34 = mrSges(6,1) * t65 - mrSges(6,3) * t119;
t26 = t32 * mrSges(6,1);
t25 = t31 * mrSges(5,2);
t18 = -pkin(4) * t126 - qJ(5) * t69 + t78;
t17 = t163 + t33;
t16 = t157 * t49 - t142;
t9 = pkin(3) * t91 + pkin(4) * t40 - qJ(5) * t39 - qJD(5) * t69;
t8 = qJD(4) * t20 + t113 * t43 - t118 * t157;
t7 = qJD(4) * t127 + t113 * t118 + t157 * t43;
t6 = pkin(4) * t32 - qJ(5) * t31 - qJD(5) * t119 + t160;
t1 = [t55 * t76 + t56 * t77 + t9 * t34 + m(4) * (t45 * t74 + t46 * t73 + t55 * t71 + t56 * t70) + t128 * Ifges(4,2) * t84 + m(6) * (t10 * t8 + t11 * t7 + t14 * t9 + t18 * t6 + t2 * t20 - t162) + t100 * (mrSges(4,1) * t91 - mrSges(4,2) * t90) + qJD(3) * (-Ifges(4,5) * t90 - Ifges(4,6) * t91) / 0.2e1 + (t180 * t40 + t181 * t39) * t110 / 0.2e1 + (t147 * t40 + t182 * t39) * t166 + t39 * t188 + (-t128 * (-Ifges(4,4) * t90 - Ifges(4,2) * t91) / 0.2e1 + 0.2e1 * t186 * qJD(2)) * qJD(1) + t78 * t25 + t145 * t7 + (-t11 * mrSges(6,2) - Ifges(5,2) * t170 + Ifges(6,3) * t169 + t185) * t40 + (Ifges(5,4) * t170 + Ifges(6,5) * t169 + t184) * t39 + t133 * t132 + t18 * t26 + t6 * (-mrSges(6,1) * t126 - mrSges(6,3) * t69) + (mrSges(5,1) * t78 + t147 * t69 - (Ifges(5,2) + Ifges(6,3)) * t126 - t148 * t20) * t32 + (t126 * t2 + t161) * mrSges(6,2) + (t126 * t4 + t161) * mrSges(5,3) + (-mrSges(6,3) * t18 - t126 * t147 - t127 * t148 + t182 * t69) * t31 + (t165 * t89 - t83 * t97) * Ifges(4,1) + m(5) * (-t12 * t8 + t13 * t7 - t162 + t20 * t4 + (t72 * t91 + t78 * t84) * pkin(3)) + (t84 * (-mrSges(5,1) * t126 + mrSges(5,2) * t69) + t91 * t35) * pkin(3) + (-t128 * t45 - t46 * t97 + t70 * t90 - t71 * t91 + t73 * t83 - t74 * t84) * mrSges(4,3) + (t128 * t83 + t89 * t164 - t84 * t97) * Ifges(4,4) + t62 * t164 + t63 * t165 - t178 * t8; t25 + t32 * mrSges(5,1) - t31 * mrSges(6,3) + t26 - m(4) * (-t124 * t71 - t70 * t89) + t76 * t124 + t89 * t77 + t132 + t145 * t65 + t178 * t119 + (-t10 * t119 + t11 * t65 + t6) * m(6) + (t119 * t12 + t13 * t65 + t160) * m(5) - t186 * qJD(1) ^ 2; t104 * t54 + t89 * t62 / 0.2e1 - Ifges(4,5) * t83 - Ifges(4,6) * t84 - t45 * mrSges(4,2) + t46 * mrSges(4,1) - t17 * t34 - t145 * t16 + t120 + t51 * t130 - t89 * (-Ifges(4,1) * t124 - t154) / 0.2e1 - t100 * (t89 * mrSges(4,1) - mrSges(4,2) * t124) - qJD(3) * (-Ifges(4,5) * t124 - Ifges(4,6) * t89) / 0.2e1 + (-Ifges(5,4) * t169 - Ifges(6,5) * t170 - t181 * t159 - t182 * t167 + t184 + t188) * t65 + ((-t157 * t5 + t113 * t4 + (-t113 * t12 + t13 * t157) * qJD(4)) * pkin(3) + t12 * t15 - t13 * t16 - t163 * t72) * m(5) - t35 * t163 + (-t136 * t31 - t151 * t32) * mrSges(5,3) + t187 * t178 + (t106 * t2 + t107 * t5 - t14 * t17 + (t104 - t16) * t11 - t187 * t10) * m(6) + (-Ifges(5,2) * t169 + Ifges(6,3) * t170 + t147 * t167 - t185) * t119 + (-t106 * t32 + t176) * mrSges(6,2) + (-Ifges(4,2) * t89 + t63 - t85) * t124 / 0.2e1 + (-t76 - t123) * t70 + (t77 + t156) * t71 + t174 * t159 + t107 * t150; -t72 * (mrSges(5,1) * t119 - mrSges(5,2) * t65) + (Ifges(6,3) * t119 - t152) * t170 - t14 * (mrSges(6,1) * t119 + mrSges(6,3) * t65) + t22 * t166 + qJD(5) * t54 - t33 * t34 + t120 + (-pkin(4) * t31 - qJ(5) * t32 + t10 * t65 + t176) * mrSges(6,2) + (t178 + t155) * t13 + (-t145 - t149) * t12 + (-t181 * t65 + t174) * t159 + (-pkin(4) * t5 + qJ(5) * t2 - t10 * t13 + t11 * t175 - t14 * t33) * m(6) + (-Ifges(5,2) * t119 + t179 - t61) * t169 + (-t182 * t65 - t153 + t21 + t60) * t167; t150 - t110 * t54 + t119 * t34 + 0.2e1 * (t5 / 0.2e1 + t11 * t159 + t14 * t166) * m(6);];
tauc = t1(:);
