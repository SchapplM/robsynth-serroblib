% Calculate joint inertia matrix for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:27
% EndTime: 2019-12-31 22:14:31
% DurationCPUTime: 1.25s
% Computational Cost: add. (1351->306), mult. (3103->422), div. (0->0), fcn. (3080->8), ass. (0->113)
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t151 = t110 ^ 2 + t113 ^ 2;
t150 = 2 * pkin(8);
t108 = sin(pkin(5));
t115 = cos(qJ(2));
t134 = t108 * t115;
t109 = cos(pkin(5));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t112 = sin(qJ(2));
t135 = t108 * t112;
t59 = t109 * t111 + t114 * t135;
t33 = t110 * t59 + t113 * t134;
t34 = -t110 * t134 + t113 * t59;
t58 = -t109 * t114 + t111 * t135;
t7 = Ifges(5,5) * t34 - Ifges(5,6) * t33 + Ifges(5,3) * t58;
t8 = Ifges(6,4) * t34 + Ifges(6,2) * t58 + Ifges(6,6) * t33;
t149 = t7 + t8;
t148 = pkin(1) * t115;
t147 = pkin(8) * t114;
t84 = pkin(7) * t135;
t60 = t109 * t148 - t84;
t146 = t60 * mrSges(3,1);
t61 = t109 * t112 * pkin(1) + pkin(7) * t134;
t145 = t61 * mrSges(3,2);
t144 = -Ifges(6,2) - Ifges(5,3);
t44 = t84 + (-pkin(2) - t148) * t109;
t16 = pkin(3) * t58 - pkin(9) * t59 + t44;
t45 = pkin(8) * t109 + t61;
t46 = (-pkin(2) * t115 - pkin(8) * t112 - pkin(1)) * t108;
t24 = t111 * t46 + t114 * t45;
t18 = -pkin(9) * t134 + t24;
t4 = t110 * t16 + t113 * t18;
t143 = -Ifges(4,5) * t59 + Ifges(4,6) * t58;
t132 = t111 * t113;
t67 = t114 * mrSges(6,1) + mrSges(6,2) * t132;
t133 = t110 * t111;
t142 = Ifges(6,4) * t132 + Ifges(6,6) * t133;
t70 = -pkin(3) * t114 - pkin(9) * t111 - pkin(2);
t43 = t110 * t70 + t113 * t147;
t75 = Ifges(5,5) * t110 + Ifges(5,6) * t113;
t141 = Ifges(4,5) * t111 + Ifges(4,6) * t114;
t140 = Ifges(5,4) * t110;
t139 = Ifges(5,4) * t113;
t138 = Ifges(6,5) * t110;
t137 = Ifges(6,5) * t113;
t136 = t113 * t70;
t131 = t151 * pkin(9) ^ 2;
t6 = Ifges(6,5) * t34 + Ifges(6,6) * t58 + Ifges(6,3) * t33;
t9 = Ifges(5,4) * t34 - Ifges(5,2) * t33 + Ifges(5,6) * t58;
t130 = t6 / 0.2e1 - t9 / 0.2e1;
t129 = Ifges(3,5) * t135 + Ifges(3,6) * t134 + Ifges(3,3) * t109;
t10 = Ifges(6,1) * t34 + Ifges(6,4) * t58 + Ifges(6,5) * t33;
t11 = Ifges(5,1) * t34 - Ifges(5,4) * t33 + Ifges(5,5) * t58;
t128 = t10 / 0.2e1 + t11 / 0.2e1;
t47 = -Ifges(6,6) * t114 + (Ifges(6,3) * t110 + t137) * t111;
t50 = -Ifges(5,6) * t114 + (-Ifges(5,2) * t110 + t139) * t111;
t127 = t47 / 0.2e1 - t50 / 0.2e1;
t51 = -Ifges(6,4) * t114 + (Ifges(6,1) * t113 + t138) * t111;
t52 = -Ifges(5,5) * t114 + (Ifges(5,1) * t113 - t140) * t111;
t126 = t51 / 0.2e1 + t52 / 0.2e1;
t74 = -Ifges(6,3) * t113 + t138;
t77 = Ifges(5,2) * t113 + t140;
t125 = t74 / 0.2e1 - t77 / 0.2e1;
t76 = Ifges(6,4) * t110 - Ifges(6,6) * t113;
t124 = t75 / 0.2e1 + t76 / 0.2e1;
t79 = Ifges(6,1) * t110 - t137;
t80 = Ifges(5,1) * t110 + t139;
t123 = t79 / 0.2e1 + t80 / 0.2e1;
t22 = -t58 * mrSges(6,1) + t34 * mrSges(6,2);
t23 = -t111 * t45 + t114 * t46;
t17 = pkin(3) * t134 - t23;
t121 = Ifges(5,5) * t132 - Ifges(5,6) * t133;
t120 = t110 * mrSges(5,1) + t113 * mrSges(5,2);
t119 = t110 * mrSges(6,1) - t113 * mrSges(6,3);
t118 = -pkin(4) * t110 + qJ(5) * t113;
t3 = -t110 * t18 + t113 * t16;
t117 = pkin(8) ^ 2;
t107 = t114 ^ 2;
t105 = t111 ^ 2;
t102 = t105 * t117;
t81 = Ifges(4,1) * t111 + Ifges(4,4) * t114;
t78 = Ifges(4,4) * t111 + Ifges(4,2) * t114;
t73 = -mrSges(4,1) * t114 + mrSges(4,2) * t111;
t72 = -mrSges(5,1) * t113 + mrSges(5,2) * t110;
t71 = -mrSges(6,1) * t113 - mrSges(6,3) * t110;
t69 = -pkin(4) * t113 - qJ(5) * t110 - pkin(3);
t68 = -mrSges(6,2) * t133 - mrSges(6,3) * t114;
t66 = -mrSges(5,1) * t114 - mrSges(5,3) * t132;
t65 = mrSges(5,2) * t114 - mrSges(5,3) * t133;
t63 = t120 * t111;
t62 = t119 * t111;
t53 = (pkin(8) - t118) * t111;
t49 = -Ifges(6,2) * t114 + t142;
t48 = -Ifges(5,3) * t114 + t121;
t42 = -t110 * t147 + t136;
t38 = -t136 + (pkin(8) * t110 + pkin(4)) * t114;
t37 = -qJ(5) * t114 + t43;
t36 = -mrSges(4,1) * t134 - mrSges(4,3) * t59;
t35 = mrSges(4,2) * t134 - mrSges(4,3) * t58;
t27 = mrSges(4,1) * t58 + mrSges(4,2) * t59;
t26 = Ifges(4,1) * t59 - Ifges(4,4) * t58 - Ifges(4,5) * t134;
t25 = Ifges(4,4) * t59 - Ifges(4,2) * t58 - Ifges(4,6) * t134;
t21 = mrSges(5,1) * t58 - mrSges(5,3) * t34;
t20 = -mrSges(5,2) * t58 - mrSges(5,3) * t33;
t19 = -mrSges(6,2) * t33 + mrSges(6,3) * t58;
t13 = mrSges(5,1) * t33 + mrSges(5,2) * t34;
t12 = mrSges(6,1) * t33 - mrSges(6,3) * t34;
t5 = pkin(4) * t33 - qJ(5) * t34 + t17;
t2 = -pkin(4) * t58 - t3;
t1 = qJ(5) * t58 + t4;
t14 = [0.2e1 * t1 * t19 + 0.2e1 * t5 * t12 + 0.2e1 * t17 * t13 + 0.2e1 * t2 * t22 + 0.2e1 * t4 * t20 + 0.2e1 * t3 * t21 + 0.2e1 * t23 * t36 + 0.2e1 * t24 * t35 + t59 * t26 + 0.2e1 * t44 * t27 + Ifges(2,3) + (t10 + t11) * t34 + (t6 - t9) * t33 + (-t25 + t149) * t58 + (t129 - 0.2e1 * t145 + 0.2e1 * t146) * t109 + ((-0.2e1 * t60 * mrSges(3,3) + Ifges(3,5) * t109 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t112) * t108) * t112 + (0.2e1 * t61 * mrSges(3,3) + Ifges(3,6) * t109 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t112 + (Ifges(3,2) + Ifges(4,3)) * t115) * t108 + t143) * t115) * t108 + m(3) * (pkin(1) ^ 2 * t108 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2 + t44 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t17 ^ 2 + t3 ^ 2 + t4 ^ 2); -t141 * t134 / 0.2e1 + t129 + t126 * t34 + t127 * t33 + (t26 / 0.2e1 - t23 * mrSges(4,3) + t128 * t113 + t130 * t110 + (-t36 + t13) * pkin(8)) * t111 + m(4) * (-pkin(2) * t44 + (-t23 * t111 + t24 * t114) * pkin(8)) + t4 * t65 + t3 * t66 + t2 * t67 + t1 * t68 + t44 * t73 + t59 * t81 / 0.2e1 + t5 * t62 + t17 * t63 + t42 * t21 + t43 * t20 + t53 * t12 + t37 * t19 + t38 * t22 - pkin(2) * t27 + t146 - t145 + (-t78 / 0.2e1 + t48 / 0.2e1 + t49 / 0.2e1) * t58 + m(5) * (pkin(8) * t111 * t17 + t3 * t42 + t4 * t43) + m(6) * (t1 * t37 + t2 * t38 + t5 * t53) + (t25 / 0.2e1 - t7 / 0.2e1 - t8 / 0.2e1 + pkin(8) * t35 + t24 * mrSges(4,3)) * t114; -0.2e1 * pkin(2) * t73 + 0.2e1 * t37 * t68 + 0.2e1 * t38 * t67 + 0.2e1 * t42 * t66 + 0.2e1 * t43 * t65 + 0.2e1 * t53 * t62 + Ifges(3,3) + (t105 + t107) * mrSges(4,3) * t150 + (-t49 - t48 + t78) * t114 + m(5) * (t42 ^ 2 + t43 ^ 2 + t102) + m(6) * (t37 ^ 2 + t38 ^ 2 + t53 ^ 2) + m(4) * (pkin(2) ^ 2 + t107 * t117 + t102) + (t63 * t150 + t81 + (t51 + t52) * t113 + (t47 - t50) * t110) * t111; -Ifges(4,3) * t134 + t23 * mrSges(4,1) - t24 * mrSges(4,2) - pkin(3) * t13 + t69 * t12 + t17 * t72 + t5 * t71 + t124 * t58 + t123 * t34 + t125 * t33 + (t1 * mrSges(6,2) + t4 * mrSges(5,3) + (t19 + t20) * pkin(9) - t130) * t113 + (t2 * mrSges(6,2) - t3 * mrSges(5,3) + (-t21 + t22) * pkin(9) + t128) * t110 + m(6) * (t5 * t69 + (t1 * t113 + t110 * t2) * pkin(9)) + m(5) * (-pkin(3) * t17 + (-t110 * t3 + t113 * t4) * pkin(9)) - t143; -pkin(3) * t63 + t69 * t62 + (m(6) * t69 + t71) * t53 + (-pkin(8) * mrSges(4,2) - t124) * t114 + (t37 * mrSges(6,2) + t43 * mrSges(5,3) - t127) * t113 + (t38 * mrSges(6,2) - t42 * mrSges(5,3) + t126) * t110 + ((t65 + t68) * t113 + (-t66 + t67) * t110 + m(6) * (t110 * t38 + t113 * t37) + m(5) * (-t110 * t42 + t113 * t43)) * pkin(9) + (t123 * t113 + t125 * t110 + (-m(5) * pkin(3) - mrSges(4,1) + t72) * pkin(8)) * t111 + t141; -0.2e1 * pkin(3) * t72 + 0.2e1 * t69 * t71 + Ifges(4,3) + (-t74 + t77) * t113 + (t79 + t80) * t110 + m(6) * (t69 ^ 2 + t131) + m(5) * (pkin(3) ^ 2 + t131) + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * pkin(9) * t151; m(6) * (-pkin(4) * t2 + qJ(5) * t1) + t1 * mrSges(6,3) + qJ(5) * t19 - t4 * mrSges(5,2) + t3 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t22 + t149; -t43 * mrSges(5,2) - pkin(4) * t67 + m(6) * (-pkin(4) * t38 + qJ(5) * t37) + qJ(5) * t68 + t37 * mrSges(6,3) - t38 * mrSges(6,1) + t42 * mrSges(5,1) + t144 * t114 + t121 + t142; t118 * mrSges(6,2) + (m(6) * t118 - t119 - t120) * pkin(9) + t76 + t75; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) - t144; m(6) * t2 + t22; m(6) * t38 + t67; (m(6) * pkin(9) + mrSges(6,2)) * t110; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
