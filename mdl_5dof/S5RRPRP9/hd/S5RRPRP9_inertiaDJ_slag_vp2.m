% Calculate time derivative of joint inertia matrix for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:43
% EndTime: 2019-12-31 20:05:48
% DurationCPUTime: 2.02s
% Computational Cost: add. (1650->293), mult. (4133->435), div. (0->0), fcn. (3407->6), ass. (0->121)
t163 = Ifges(6,4) + Ifges(5,5);
t161 = Ifges(6,6) - Ifges(5,6);
t162 = -Ifges(6,2) - Ifges(5,3);
t115 = cos(pkin(8));
t118 = cos(qJ(4));
t114 = sin(pkin(8));
t116 = sin(qJ(4));
t144 = t114 * t116;
t92 = -t118 * t115 + t144;
t84 = t92 * qJD(4);
t93 = t114 * t118 + t115 * t116;
t85 = t93 * qJD(4);
t160 = t161 * t85 - t163 * t84;
t159 = m(6) * qJ(5) + mrSges(6,3);
t119 = cos(qJ(2));
t117 = sin(qJ(2));
t141 = t115 * t117;
t97 = -pkin(2) * t119 - t117 * qJ(3) - pkin(1);
t90 = t115 * t97;
t52 = -pkin(7) * t141 + t90 + (-pkin(6) * t114 - pkin(3)) * t119;
t143 = t114 * t117;
t140 = t115 * t119;
t71 = pkin(6) * t140 + t114 * t97;
t61 = -pkin(7) * t143 + t71;
t150 = t116 * t52 + t118 * t61;
t138 = qJD(2) * t117;
t131 = pkin(6) * t138;
t83 = -t117 * qJD(3) + (pkin(2) * t117 - qJ(3) * t119) * qJD(2);
t59 = t114 * t131 + t115 * t83;
t30 = (pkin(3) * t117 - pkin(7) * t140) * qJD(2) + t59;
t142 = t114 * t119;
t72 = t114 * t83;
t42 = t72 + (-pkin(6) * t141 - pkin(7) * t142) * qJD(2);
t4 = -qJD(4) * t150 - t116 * t42 + t118 * t30;
t158 = 2 * m(4);
t157 = 2 * m(5);
t156 = 2 * m(6);
t155 = -0.2e1 * pkin(1);
t154 = 0.2e1 * pkin(6);
t153 = t115 / 0.2e1;
t151 = pkin(7) + qJ(3);
t147 = mrSges(4,2) * t115;
t146 = Ifges(4,4) * t114;
t145 = Ifges(4,4) * t115;
t137 = qJD(2) * t119;
t130 = t114 * t137;
t76 = mrSges(4,1) * t130 + t137 * t147;
t110 = pkin(6) * t137;
t88 = pkin(3) * t130 + t110;
t96 = pkin(3) * t143 + t117 * pkin(6);
t136 = qJD(4) * t116;
t135 = qJD(4) * t117;
t134 = qJD(4) * t118;
t133 = t116 * qJD(3);
t132 = t118 * qJD(3);
t107 = -pkin(3) * t115 - pkin(2);
t38 = -t135 * t93 - t137 * t92;
t39 = t134 * t141 - t135 * t144 + t137 * t93;
t11 = t39 * mrSges(5,1) + t38 * mrSges(5,2);
t46 = t85 * mrSges(5,1) - t84 * mrSges(5,2);
t10 = t39 * mrSges(6,1) - t38 * mrSges(6,3);
t45 = t85 * mrSges(6,1) + t84 * mrSges(6,3);
t127 = qJD(4) * t151;
t99 = t151 * t115;
t24 = t115 * t132 - t99 * t136 + (-t118 * t127 - t133) * t114;
t25 = t115 * t133 + t99 * t134 + (-t116 * t127 + t132) * t114;
t128 = t151 * t114;
t62 = t116 * t99 + t118 * t128;
t63 = -t116 * t128 + t118 * t99;
t129 = t63 * t24 + t25 * t62;
t21 = -mrSges(6,1) * t138 + t38 * mrSges(6,2);
t126 = t115 * Ifges(4,1) - t146;
t125 = -t114 * Ifges(4,2) + t145;
t124 = -Ifges(4,5) * t115 + Ifges(4,6) * t114;
t15 = -t116 * t61 + t118 * t52;
t121 = t162 * t138 - t161 * t39 - t163 * t38;
t3 = t116 * t30 + t118 * t42 + t52 * t134 - t136 * t61;
t95 = -mrSges(4,1) * t119 - mrSges(4,3) * t141;
t94 = mrSges(4,2) * t119 - mrSges(4,3) * t143;
t87 = (mrSges(4,1) * t117 - mrSges(4,3) * t140) * qJD(2);
t86 = (-mrSges(4,2) * t117 - mrSges(4,3) * t142) * qJD(2);
t75 = t92 * t117;
t74 = t93 * t117;
t70 = -pkin(6) * t142 + t90;
t69 = (t117 * Ifges(4,5) + t119 * t126) * qJD(2);
t68 = (t117 * Ifges(4,6) + t119 * t125) * qJD(2);
t67 = mrSges(6,1) * t119 - t75 * mrSges(6,2);
t66 = -mrSges(5,1) * t119 + t75 * mrSges(5,3);
t65 = mrSges(5,2) * t119 - t74 * mrSges(5,3);
t64 = -t74 * mrSges(6,2) - mrSges(6,3) * t119;
t60 = -t115 * t131 + t72;
t58 = Ifges(5,1) * t93 - Ifges(5,4) * t92;
t57 = Ifges(6,1) * t93 + Ifges(6,5) * t92;
t56 = Ifges(5,4) * t93 - Ifges(5,2) * t92;
t55 = Ifges(6,5) * t93 + Ifges(6,3) * t92;
t54 = mrSges(6,1) * t92 - mrSges(6,3) * t93;
t51 = pkin(4) * t92 - qJ(5) * t93 + t107;
t50 = -Ifges(5,1) * t84 - Ifges(5,4) * t85;
t49 = -Ifges(6,1) * t84 + Ifges(6,5) * t85;
t48 = -Ifges(5,4) * t84 - Ifges(5,2) * t85;
t47 = -Ifges(6,5) * t84 + Ifges(6,3) * t85;
t41 = mrSges(6,1) * t74 + mrSges(6,3) * t75;
t29 = -Ifges(5,1) * t75 - Ifges(5,4) * t74 - Ifges(5,5) * t119;
t28 = -Ifges(6,1) * t75 - Ifges(6,4) * t119 + Ifges(6,5) * t74;
t27 = -Ifges(5,4) * t75 - Ifges(5,2) * t74 - Ifges(5,6) * t119;
t26 = -Ifges(6,5) * t75 - Ifges(6,6) * t119 + Ifges(6,3) * t74;
t22 = -mrSges(5,2) * t138 - t39 * mrSges(5,3);
t20 = mrSges(5,1) * t138 - t38 * mrSges(5,3);
t19 = -t39 * mrSges(6,2) + mrSges(6,3) * t138;
t18 = pkin(4) * t85 + qJ(5) * t84 - qJD(5) * t93;
t17 = pkin(4) * t74 + qJ(5) * t75 + t96;
t14 = pkin(4) * t119 - t15;
t12 = -qJ(5) * t119 + t150;
t9 = Ifges(5,1) * t38 - Ifges(5,4) * t39 + Ifges(5,5) * t138;
t8 = Ifges(6,1) * t38 + Ifges(6,4) * t138 + Ifges(6,5) * t39;
t7 = Ifges(5,4) * t38 - Ifges(5,2) * t39 + Ifges(5,6) * t138;
t6 = Ifges(6,5) * t38 + Ifges(6,6) * t138 + Ifges(6,3) * t39;
t5 = t39 * pkin(4) - t38 * qJ(5) + qJD(5) * t75 + t88;
t2 = -pkin(4) * t138 - t4;
t1 = qJ(5) * t138 - qJD(5) * t119 + t3;
t13 = [0.2e1 * t60 * t94 + 0.2e1 * t59 * t95 + 0.2e1 * t96 * t11 + 0.2e1 * t71 * t86 + 0.2e1 * t70 * t87 + t74 * t6 - t74 * t7 - t75 * t8 - t75 * t9 + 0.2e1 * t1 * t64 + 0.2e1 * t3 * t65 + 0.2e1 * t4 * t66 + 0.2e1 * t2 * t67 + 0.2e1 * t5 * t41 + 0.2e1 * t12 * t19 + 0.2e1 * t15 * t20 + 0.2e1 * t14 * t21 + 0.2e1 * t17 * t10 + ((mrSges(3,2) * t155 + 0.2e1 * (Ifges(3,4) + t124) * t119) * t119 + (mrSges(3,1) * t155 + (-0.2e1 * Ifges(3,4) - t124) * t117 - t163 * t75 + t161 * t74 + (-(2 * Ifges(4,3)) + pkin(6) ^ 2 * t158 - (2 * Ifges(3,2)) + (2 * Ifges(3,1)) + (mrSges(4,1) * t114 + t147) * t154 + t115 * t126 - t114 * t125 + t162) * t119) * t117) * qJD(2) + 0.2e1 * t150 * t22 + (t15 * t4 + t150 * t3 + t88 * t96) * t157 + 0.2e1 * t88 * (mrSges(5,1) * t74 - mrSges(5,2) * t75) + (t1 * t12 + t14 * t2 + t17 * t5) * t156 + (t70 * t59 + t71 * t60) * t158 + t121 * t119 + (t28 + t29) * t38 + (-t27 + t26) * t39 + (-t114 * t68 + t115 * t69 + t76 * t154) * t117; m(6) * (t1 * t63 + t24 * t12 + t25 * t14 + t18 * t17 + t2 * t62 + t51 * t5) - (t49 / 0.2e1 + t50 / 0.2e1) * t75 + t96 * t46 + t107 * t11 - pkin(2) * t76 + t5 * t54 + t18 * t41 + t17 * t45 + t51 * t10 + (t47 / 0.2e1 - t48 / 0.2e1) * t74 - t160 * t119 / 0.2e1 + (t55 / 0.2e1 - t56 / 0.2e1) * t39 + (t57 / 0.2e1 + t58 / 0.2e1) * t38 + m(5) * (t107 * t88 - t25 * t15 + t150 * t24 + t3 * t63 - t4 * t62) + (-t1 * t92 - t12 * t85 - t14 * t84 + t2 * t93) * mrSges(6,2) + (t15 * t84 - t150 * t85 - t3 * t92 - t4 * t93) * mrSges(5,3) + m(4) * ((-t114 * t70 + t115 * t71) * qJD(3) + (-t114 * t59 + t115 * t60) * qJ(3)) + (t19 + t22) * t63 - (t28 / 0.2e1 + t29 / 0.2e1) * t84 + (-t27 / 0.2e1 + t26 / 0.2e1) * t85 + (t88 * mrSges(5,1) + t6 / 0.2e1 - t7 / 0.2e1) * t92 + (t88 * mrSges(5,2) + t8 / 0.2e1 + t9 / 0.2e1) * t93 + (t69 / 0.2e1 - qJ(3) * t87 - qJD(3) * t95 - t59 * mrSges(4,3)) * t114 + (t68 / 0.2e1 + qJ(3) * t86 + qJD(3) * t94 + t60 * mrSges(4,3)) * t115 + (t64 + t65) * t24 + (-t66 + t67) * t25 + ((Ifges(4,5) * t114 / 0.2e1 + Ifges(4,6) * t153 - Ifges(3,6) + pkin(6) * mrSges(3,2) + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t93 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t92) * t117 + (Ifges(3,5) + (Ifges(4,1) * t114 + t145) * t153 - t114 * (Ifges(4,2) * t115 + t146) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t115 + mrSges(4,2) * t114 - mrSges(3,1)) * pkin(6)) * t119) * qJD(2) + (t21 - t20) * t62; 0.2e1 * t107 * t46 + 0.2e1 * t18 * t54 + 0.2e1 * t51 * t45 + (t49 + t50) * t93 + (t47 - t48) * t92 + (t55 - t56) * t85 - (t57 + t58) * t84 + (t51 * t18 + t129) * t156 + t129 * t157 + (qJ(3) * t158 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t114 ^ 2 + t115 ^ 2) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * (-t24 * t92 + t25 * t93 - t62 * t84 - t63 * t85); m(4) * t110 + m(5) * t88 + m(6) * t5 + t10 + t11 + t76; m(6) * t18 + t45 + t46; 0; -pkin(4) * t21 + m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t12) + qJD(5) * t64 + qJ(5) * t19 + t1 * mrSges(6,3) - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t2 * mrSges(6,1) - t121; m(6) * qJD(5) * t63 + (pkin(4) * t84 - qJ(5) * t85 - qJD(5) * t92) * mrSges(6,2) + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t25 + (-mrSges(5,2) + t159) * t24 + t160; 0; 0.2e1 * t159 * qJD(5); m(6) * t2 + t21; m(6) * t25 - t84 * mrSges(6,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
