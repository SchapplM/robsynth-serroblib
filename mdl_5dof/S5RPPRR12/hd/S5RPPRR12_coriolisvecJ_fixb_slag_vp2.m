% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:07:02
% DurationCPUTime: 2.09s
% Computational Cost: add. (2326->247), mult. (5317->340), div. (0->0), fcn. (3468->6), ass. (0->121)
t172 = -Ifges(5,1) / 0.2e1;
t87 = sin(qJ(5));
t88 = cos(qJ(5));
t101 = Ifges(6,5) * t88 - Ifges(6,6) * t87;
t139 = Ifges(6,4) * t88;
t103 = -Ifges(6,2) * t87 + t139;
t140 = Ifges(6,4) * t87;
t105 = Ifges(6,1) * t88 - t140;
t106 = mrSges(6,1) * t87 + mrSges(6,2) * t88;
t144 = sin(qJ(4));
t145 = cos(qJ(4));
t86 = -pkin(1) - qJ(3);
t70 = t86 * qJD(1) + qJD(2);
t119 = -pkin(6) * qJD(1) + t70;
t84 = sin(pkin(8));
t53 = t119 * t84;
t85 = cos(pkin(8));
t54 = t119 * t85;
t35 = t144 * t54 + t145 * t53;
t29 = qJD(4) * pkin(7) + t35;
t63 = t144 * t85 + t145 * t84;
t57 = t63 * qJD(1);
t62 = t144 * t84 - t145 * t85;
t58 = t62 * qJD(1);
t83 = qJD(1) * qJ(2);
t79 = qJD(3) + t83;
t80 = t84 * pkin(3);
t67 = qJD(1) * t80 + t79;
t30 = pkin(4) * t57 + pkin(7) * t58 + t67;
t8 = -t29 * t87 + t30 * t88;
t9 = t29 * t88 + t30 * t87;
t110 = t8 * t88 + t9 * t87;
t46 = qJD(4) * t87 - t58 * t88;
t141 = Ifges(6,4) * t46;
t45 = qJD(4) * t88 + t58 * t87;
t56 = qJD(5) + t57;
t15 = Ifges(6,2) * t45 + Ifges(6,6) * t56 + t141;
t153 = t88 / 0.2e1;
t155 = t56 / 0.2e1;
t156 = t46 / 0.2e1;
t157 = t45 / 0.2e1;
t44 = Ifges(6,4) * t45;
t16 = Ifges(6,1) * t46 + Ifges(6,5) * t56 + t44;
t34 = -t144 * t53 + t145 * t54;
t28 = -qJD(4) * pkin(4) - t34;
t171 = -t110 * mrSges(6,3) + t101 * t155 + t103 * t157 + t105 * t156 + t28 * t106 - t87 * t15 / 0.2e1 + t16 * t153;
t93 = t63 * qJD(3);
t19 = -qJD(1) * t93 + t34 * qJD(4);
t125 = qJD(1) * qJD(2);
t51 = qJD(4) * t57;
t52 = qJD(4) * t58;
t33 = -pkin(4) * t52 + pkin(7) * t51 + t125;
t1 = qJD(5) * t8 + t19 * t88 + t33 * t87;
t2 = -qJD(5) * t9 - t19 * t87 + t33 * t88;
t112 = t1 * t88 - t2 * t87;
t31 = -mrSges(6,2) * t56 + mrSges(6,3) * t45;
t32 = mrSges(6,1) * t56 - mrSges(6,3) * t46;
t159 = -m(6) * t110 - t87 * t31 - t88 * t32;
t24 = t45 * qJD(5) - t51 * t88;
t17 = -mrSges(6,1) * t52 - mrSges(6,3) * t24;
t25 = -t46 * qJD(5) + t51 * t87;
t18 = mrSges(6,2) * t52 + mrSges(6,3) * t25;
t170 = m(6) * t112 + t159 * qJD(5) - t87 * t17 + t88 * t18;
t169 = Ifges(5,2) / 0.2e1;
t55 = Ifges(5,4) * t57;
t160 = t58 * t172 - t55 / 0.2e1;
t168 = t67 * mrSges(5,2) + Ifges(5,5) * qJD(4) + t160 + t171;
t166 = t57 * t169;
t129 = t84 ^ 2 + t85 ^ 2;
t163 = mrSges(4,3) * t129;
t162 = t62 * qJD(3);
t109 = -t8 * t87 + t9 * t88;
t161 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t24 + Ifges(6,6) * t25;
t154 = t87 / 0.2e1;
t146 = -pkin(6) + t86;
t143 = m(3) * qJ(2);
t20 = -t162 * qJD(1) + t35 * qJD(4);
t65 = t146 * t84;
t66 = t146 * t85;
t42 = t144 * t65 - t145 * t66;
t138 = t20 * t42;
t137 = t20 * t62;
t131 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t45 + mrSges(6,2) * t46 - t58 * mrSges(5,3);
t108 = mrSges(4,1) * t84 + mrSges(4,2) * t85;
t130 = mrSges(5,1) * t57 - mrSges(5,2) * t58 + qJD(1) * t108;
t128 = m(4) * qJD(3);
t76 = qJ(2) + t80;
t122 = -t52 * mrSges(5,1) - t51 * mrSges(5,2);
t121 = qJD(4) * t145;
t120 = qJD(4) * t144;
t118 = qJD(1) * t129;
t111 = t1 * t87 + t2 * t88;
t107 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t104 = Ifges(6,1) * t87 + t139;
t102 = Ifges(6,2) * t88 + t140;
t100 = Ifges(6,5) * t87 + Ifges(6,6) * t88;
t98 = t88 * t31 - t87 * t32;
t59 = -t85 * t120 - t84 * t121;
t60 = -t84 * t120 + t85 * t121;
t96 = t34 * t59 + t35 * t60;
t40 = pkin(4) * t63 + pkin(7) * t62 + t76;
t43 = t144 * t66 + t145 * t65;
t12 = t40 * t88 - t43 * t87;
t13 = t40 * t87 + t43 * t88;
t47 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t57;
t95 = t47 + t98;
t92 = t67 * mrSges(5,1) + t8 * mrSges(6,1) - t9 * mrSges(6,2) + Ifges(5,4) * t58 + t46 * Ifges(6,5) - Ifges(5,6) * qJD(4) + t45 * Ifges(6,6) + t56 * Ifges(6,3) + t166;
t89 = qJD(1) ^ 2;
t50 = Ifges(6,3) * t52;
t41 = -pkin(4) * t58 + pkin(7) * t57;
t38 = pkin(4) * t60 - pkin(7) * t59 + qJD(2);
t27 = t43 * qJD(4) - t162;
t26 = -t42 * qJD(4) - t93;
t11 = t34 * t88 + t41 * t87;
t10 = -t34 * t87 + t41 * t88;
t7 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t6 = t24 * Ifges(6,1) + t25 * Ifges(6,4) - t52 * Ifges(6,5);
t5 = t24 * Ifges(6,4) + t25 * Ifges(6,2) - t52 * Ifges(6,6);
t4 = -t13 * qJD(5) - t26 * t87 + t38 * t88;
t3 = t12 * qJD(5) + t26 * t88 + t38 * t87;
t14 = [t76 * t122 + t26 * t47 + t42 * t7 + t3 * t31 + t4 * t32 + t12 * t17 + t13 * t18 + (t166 + t92) * t60 + (t160 + t168) * t59 + t131 * t27 + 0.2e1 * qJD(3) * t163 * qJD(1) + m(5) * (t19 * t43 + t26 * t35 - t27 * t34 + t138) + m(6) * (t1 * t13 + t12 * t2 + t27 * t28 + t3 * t9 + t4 * t8 + t138) + (-t86 * t118 - t129 * t70) * t128 + (-t42 * t51 + t43 * t52 - t96) * mrSges(5,3) + (-t50 / 0.2e1 + Ifges(5,4) * t51 + mrSges(5,1) * t125 - t19 * mrSges(5,3) - (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t52 + t161) * t63 + (-t24 * t105 / 0.2e1 - t25 * t103 / 0.2e1 + Ifges(5,1) * t51 - t88 * t6 / 0.2e1 + t5 * t154 - mrSges(5,2) * t125 + (-mrSges(5,3) - t106) * t20 + t111 * mrSges(6,3) + (t109 * mrSges(6,3) + t100 * t155 + t102 * t157 + t104 * t156 + t28 * t107 + t15 * t153 + t16 * t154) * qJD(5) - (-t101 / 0.2e1 + Ifges(5,4)) * t52) * t62 + (t130 + ((2 * mrSges(3,3)) + t108 + 0.2e1 * t143) * qJD(1) + m(4) * (t79 + t83) + m(5) * (qJD(1) * t76 + t67)) * qJD(2); (-mrSges(3,3) - t143) * t89 + (-t51 * mrSges(5,3) + t7) * t62 - t131 * t59 + t95 * t60 + m(5) * (t96 + t137) + m(6) * (t109 * t60 - t28 * t59 + t137) + (m(5) * t19 + t52 * mrSges(5,3) + t170) * t63 + (-m(4) * t79 - m(5) * t67 - t129 * t128 - t130 + t159) * qJD(1); t88 * t17 + t87 * t18 + t131 * t58 + t98 * qJD(5) + (m(4) + m(5)) * t125 - t89 * t163 + t95 * t57 - m(5) * (t34 * t58 - t35 * t57) + m(4) * t70 * t118 + t122 + (t56 * t109 + t28 * t58 + t111) * m(6); t24 * t104 / 0.2e1 + t25 * t102 / 0.2e1 + t5 * t153 + t6 * t154 - Ifges(5,5) * t51 - t34 * t47 - t11 * t31 - t10 * t32 - t19 * mrSges(5,2) - pkin(4) * t7 - t131 * t35 + (-mrSges(5,1) + t107) * t20 + t112 * mrSges(6,3) - (t35 * mrSges(5,3) - t92) * t58 - (t55 / 0.2e1 + t34 * mrSges(5,3) - (t172 + t169) * t58 - t168) * t57 + t171 * qJD(5) - (t100 / 0.2e1 - Ifges(5,6)) * t52 + (-pkin(4) * t20 - t10 * t8 - t11 * t9 - t28 * t35) * m(6) + t170 * pkin(7); -t50 - t28 * (mrSges(6,1) * t46 + mrSges(6,2) * t45) - t46 * (Ifges(6,1) * t45 - t141) / 0.2e1 + t15 * t156 - t56 * (Ifges(6,5) * t45 - Ifges(6,6) * t46) / 0.2e1 - t8 * t31 + t9 * t32 + (t45 * t8 + t46 * t9) * mrSges(6,3) - (-Ifges(6,2) * t46 + t16 + t44) * t45 / 0.2e1 + t161;];
tauc = t14(:);
