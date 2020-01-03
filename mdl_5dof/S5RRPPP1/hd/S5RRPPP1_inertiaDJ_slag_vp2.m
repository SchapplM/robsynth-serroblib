% Calculate time derivative of joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:23:27
% DurationCPUTime: 2.00s
% Computational Cost: add. (1100->265), mult. (3201->395), div. (0->0), fcn. (2627->6), ass. (0->126)
t149 = Ifges(4,1) + Ifges(6,3);
t153 = Ifges(4,4) - Ifges(6,6);
t152 = Ifges(6,4) + Ifges(5,5);
t151 = Ifges(4,5) + Ifges(6,5);
t150 = Ifges(6,2) + Ifges(5,3);
t147 = -Ifges(5,6) + Ifges(6,6);
t100 = sin(qJ(2));
t121 = qJD(2) * t100;
t97 = sin(pkin(5));
t112 = t97 * t121;
t101 = cos(qJ(2));
t98 = cos(pkin(8));
t127 = t100 * t98;
t96 = sin(pkin(8));
t99 = cos(pkin(5));
t64 = (t101 * t96 + t127 * t99) * qJD(2);
t133 = t96 * t99;
t65 = qJD(2) * t101 * t98 - t121 * t133;
t146 = t152 * t112 + t147 * t65 + t150 * t64;
t145 = t151 * t112 + t149 * t65 - t153 * t64;
t128 = qJ(3) * t97;
t79 = -pkin(2) * t101 - t100 * t128 - pkin(1);
t110 = qJ(3) * t99 + pkin(7);
t86 = t110 * t100;
t144 = (t79 * t97 - t86 * t99) * t98;
t134 = t96 * t97;
t120 = qJD(3) * t100;
t126 = t101 * t97;
t53 = -t97 * t120 + (pkin(2) * t100 - qJ(3) * t126) * qJD(2);
t105 = qJD(2) * t110;
t125 = t101 * t99;
t54 = qJD(3) * t125 - t100 * t105;
t55 = -t101 * t105 - t120 * t99;
t10 = t55 * t133 + t53 * t134 + t98 * t54;
t5 = -(qJ(4) * t121 - qJD(4) * t101) * t97 - t10;
t143 = 2 * m(4);
t142 = 2 * m(5);
t141 = 2 * m(6);
t17 = t99 * t53 - t55 * t97;
t140 = 0.2e1 * t17;
t135 = t53 * t98;
t132 = t97 * t98;
t131 = t98 * t99;
t130 = pkin(3) + qJ(5);
t87 = t110 * t101;
t75 = t96 * t87;
t129 = pkin(3) * t126 + t75;
t74 = pkin(2) * t133 + t98 * t128;
t124 = qJD(3) * t96;
t123 = qJD(3) * t98;
t122 = qJD(4) * t96;
t119 = 0.2e1 * t96;
t118 = 0.2e1 * t101;
t117 = -Ifges(4,4) + t147;
t116 = Ifges(5,4) - t151;
t115 = -Ifges(4,6) + t152;
t16 = -t86 * t133 + t79 * t134 + t98 * t87;
t114 = t97 * t124;
t113 = -pkin(2) * t98 - pkin(3);
t37 = t65 * mrSges(5,1) + mrSges(5,2) * t112;
t36 = -t64 * mrSges(6,1) + mrSges(6,2) * t112;
t27 = -t65 * mrSges(6,2) + t64 * mrSges(6,3);
t111 = -qJ(4) * t96 - pkin(2);
t32 = t99 * t79 + t86 * t97;
t44 = t96 * t54;
t109 = -t131 * t55 + t44;
t56 = -qJ(4) * t99 - t74;
t108 = Ifges(4,5) / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t107 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1 - Ifges(4,6) / 0.2e1;
t34 = t65 * mrSges(6,1) - mrSges(6,3) * t112;
t70 = t125 * t96 + t127;
t104 = -qJ(4) * t70 + t32;
t13 = qJ(4) * t126 - t16;
t102 = -qJ(4) * t65 - qJD(4) * t70 + t17;
t92 = t96 * t128;
t85 = mrSges(5,1) * t134 + mrSges(5,2) * t99;
t84 = mrSges(6,1) * t132 + mrSges(6,2) * t99;
t83 = -mrSges(5,1) * t132 - mrSges(5,3) * t99;
t82 = mrSges(6,1) * t134 - mrSges(6,3) * t99;
t81 = -mrSges(4,2) * t99 + mrSges(4,3) * t132;
t80 = mrSges(4,1) * t99 - mrSges(4,3) * t134;
t78 = qJD(4) * t99 + t123 * t97;
t77 = -qJD(5) * t99 + t114;
t73 = (-mrSges(6,2) * t96 - mrSges(6,3) * t98) * t97;
t72 = (mrSges(5,2) * t98 - mrSges(5,3) * t96) * t97;
t71 = pkin(2) * t131 - t92;
t69 = t100 * t96 - t125 * t98;
t68 = (-qJD(5) * t98 - t122) * t97;
t62 = (-pkin(3) * t98 + t111) * t97;
t61 = t113 * t99 + t92;
t60 = t65 * mrSges(5,3);
t58 = t65 * mrSges(4,2);
t52 = -mrSges(4,1) * t126 - t70 * mrSges(4,3);
t51 = mrSges(4,2) * t126 - t69 * mrSges(4,3);
t50 = t70 * mrSges(5,1) - mrSges(5,2) * t126;
t49 = -t69 * mrSges(6,1) - mrSges(6,2) * t126;
t48 = t69 * mrSges(5,1) + mrSges(5,3) * t126;
t47 = t70 * mrSges(6,1) + mrSges(6,3) * t126;
t43 = (-t130 * t98 + t111) * t97;
t42 = pkin(4) * t132 - t56;
t39 = mrSges(4,1) * t112 - mrSges(4,3) * t65;
t38 = -mrSges(4,2) * t112 - mrSges(4,3) * t64;
t35 = mrSges(5,1) * t64 - mrSges(5,3) * t112;
t33 = pkin(4) * t134 + t92 + (-qJ(5) + t113) * t99;
t31 = -mrSges(5,2) * t69 - mrSges(5,3) * t70;
t30 = -mrSges(6,2) * t70 + mrSges(6,3) * t69;
t29 = -t64 * mrSges(5,2) - t60;
t28 = t64 * mrSges(4,1) + t58;
t25 = t65 * Ifges(4,4) - t64 * Ifges(4,2) + Ifges(4,6) * t112;
t24 = Ifges(4,5) * t65 - Ifges(4,6) * t64 + Ifges(4,3) * t112;
t23 = Ifges(5,1) * t112 - Ifges(5,4) * t65 + Ifges(5,5) * t64;
t22 = Ifges(6,1) * t112 + Ifges(6,4) * t64 + Ifges(6,5) * t65;
t21 = Ifges(5,4) * t112 - t65 * Ifges(5,2) + t64 * Ifges(5,6);
t15 = -t75 + t144;
t14 = t129 - t144;
t12 = pkin(3) * t69 + t104;
t11 = -t69 * pkin(4) - t13;
t9 = -t44 + (t53 * t97 + t55 * t99) * t98;
t8 = t130 * t69 + t104;
t7 = t86 * t131 + t70 * pkin(4) + (qJ(5) * t101 - t79 * t98) * t97 + t129;
t6 = (-pkin(3) * t121 - t135) * t97 + t109;
t4 = pkin(3) * t64 + t102;
t3 = -t64 * pkin(4) - t5;
t2 = t65 * pkin(4) + (qJD(5) * t101 - t121 * t130 - t135) * t97 + t109;
t1 = qJD(5) * t69 + t130 * t64 + t102;
t18 = [((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t101) * t118 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t100 + (Ifges(3,1) - Ifges(3,2)) * t118 + (-t116 * t70 + t115 * t69 + (-Ifges(5,1) - Ifges(6,1) - Ifges(4,3)) * t126) * t97) * t100) * qJD(2) + (mrSges(4,2) * t140 + t145 - t21) * t70 + (mrSges(4,1) * t140 + t146 - t25) * t69 + ((Ifges(5,2) + t149) * t70 + t117 * t69) * t65 + (t117 * t70 + (Ifges(4,2) + t150) * t69) * t64 + 0.2e1 * t6 * t50 + 0.2e1 * t10 * t51 + 0.2e1 * t9 * t52 + 0.2e1 * t32 * t28 + 0.2e1 * t7 * t34 + 0.2e1 * t13 * t35 + 0.2e1 * t11 * t36 + 0.2e1 * t14 * t37 + 0.2e1 * t16 * t38 + 0.2e1 * t15 * t39 + 0.2e1 * t2 * t47 + 0.2e1 * t5 * t48 + 0.2e1 * t3 * t49 + 0.2e1 * t8 * t27 + 0.2e1 * t12 * t29 + 0.2e1 * t1 * t30 + 0.2e1 * t4 * t31 + (-t115 * t64 + t116 * t65 - t22 - t23 - t24) * t126 + (t1 * t8 + t11 * t3 + t2 * t7) * t141 + (t12 * t4 + t13 * t5 + t14 * t6) * t142 + (t10 * t16 + t15 * t9 + t17 * t32) * t143; m(6) * (t43 * t1 + t11 * t78 + t33 * t2 + t42 * t3 + t68 * t8 + t7 * t77) + (Ifges(3,5) * t101 - Ifges(3,6) * t100 + (-mrSges(3,1) * t101 + mrSges(3,2) * t100) * pkin(7)) * qJD(2) + (t17 * (-mrSges(4,1) * t98 + mrSges(4,2) * t96) - t64 * (Ifges(4,4) * t96 + Ifges(4,2) * t98) / 0.2e1 - t65 * (-Ifges(5,2) * t96 - Ifges(5,6) * t98) / 0.2e1 + t98 * t25 / 0.2e1 - t96 * t21 / 0.2e1 - pkin(2) * t28 - t31 * t122 + (t98 * t51 + (t50 - t52) * t96) * qJD(3) + m(4) * (-pkin(2) * t17 + t123 * t16 - t124 * t15) + m(5) * (-t12 * t122 + t124 * t14) + ((-t107 * t98 + t108 * t96) * t97 + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t99) * t121 + (t147 * t96 - t150 * t98) * t64 / 0.2e1 + (t149 * t96 + t153 * t98) * t65 / 0.2e1 + t145 * t96 / 0.2e1 - t146 * t98 / 0.2e1) * t97 + (t49 - t48) * t78 + t71 * t39 + t4 * t72 + t1 * t73 + t74 * t38 + t77 * t47 + t9 * t80 + t10 * t81 + t2 * t82 + t5 * t83 + t3 * t84 + t6 * t85 + t68 * t30 + t56 * t35 + t61 * t37 + t62 * t29 + t33 * t34 + t42 * t36 + t43 * t27 + m(4) * (t10 * t74 + t71 * t9) + m(5) * (-t13 * t78 + t4 * t62 + t5 * t56 + t6 * t61) + (t22 / 0.2e1 + t23 / 0.2e1 + t24 / 0.2e1 + t108 * t65 + t107 * t64) * t99; 0.2e1 * t68 * t73 + 0.2e1 * t77 * t82 + (t33 * t77 + t43 * t68) * t141 + ((-m(5) * t62 - t72) * qJD(4) * t119 + (0.2e1 * t98 * t81 + (-t80 + t85) * t119 + (-t71 * t96 + t74 * t98) * t143 + t61 * t96 * t142) * qJD(3)) * t97 + (-0.2e1 * m(5) * t56 + t42 * t141 - 0.2e1 * t83 + 0.2e1 * t84) * t78; t58 - t60 + (-mrSges(5,2) + mrSges(4,1)) * t64 + m(4) * t17 + m(5) * t4 + m(6) * t1 + t27; -m(5) * t122 * t97 + m(6) * t68; 0; m(5) * t6 + m(6) * t2 + t34 + t37; m(5) * t114 + m(6) * t77; 0; 0; m(6) * t3 + t36; m(6) * t78; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;
