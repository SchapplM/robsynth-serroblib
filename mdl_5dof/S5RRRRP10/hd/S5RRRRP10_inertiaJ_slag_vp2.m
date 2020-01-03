% Calculate joint inertia matrix for
% S5RRRRP10
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:11
% EndTime: 2019-12-31 22:08:14
% DurationCPUTime: 1.05s
% Computational Cost: add. (1336->303), mult. (3072->418), div. (0->0), fcn. (3077->8), ass. (0->115)
t151 = 2 * pkin(8);
t150 = -2 * mrSges(6,3);
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t111 = sin(pkin(5));
t118 = cos(qJ(2));
t134 = t111 * t118;
t112 = cos(pkin(5));
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t115 = sin(qJ(2));
t135 = t111 * t115;
t59 = t112 * t114 + t117 * t135;
t35 = -t113 * t59 - t116 * t134;
t36 = -t113 * t134 + t116 * t59;
t58 = -t112 * t117 + t114 * t135;
t6 = Ifges(6,5) * t36 + Ifges(6,6) * t35 + Ifges(6,3) * t58;
t7 = Ifges(5,5) * t36 + Ifges(5,6) * t35 + Ifges(5,3) * t58;
t149 = t6 + t7;
t74 = -mrSges(5,1) * t116 + mrSges(5,2) * t113;
t148 = -m(5) * pkin(3) + t74;
t146 = pkin(1) * t118;
t145 = pkin(8) * t117;
t87 = pkin(7) * t135;
t60 = t112 * t146 - t87;
t144 = t60 * mrSges(3,1);
t61 = t112 * t115 * pkin(1) + pkin(7) * t134;
t143 = t61 * mrSges(3,2);
t142 = -Ifges(5,3) - Ifges(6,3);
t141 = -qJ(5) - pkin(9);
t45 = t87 + (-pkin(2) - t146) * t112;
t17 = pkin(3) * t58 - pkin(9) * t59 + t45;
t46 = pkin(8) * t112 + t61;
t47 = (-pkin(2) * t118 - pkin(8) * t115 - pkin(1)) * t111;
t25 = t114 * t47 + t117 * t46;
t19 = -pkin(9) * t134 + t25;
t4 = t113 * t17 + t116 * t19;
t140 = -Ifges(4,5) * t59 + Ifges(4,6) * t58;
t132 = t114 * t116;
t133 = t113 * t114;
t62 = mrSges(6,1) * t133 + mrSges(6,2) * t132;
t71 = -pkin(3) * t117 - pkin(9) * t114 - pkin(2);
t44 = t113 * t71 + t116 * t145;
t139 = Ifges(5,4) * t113;
t138 = Ifges(5,4) * t116;
t137 = Ifges(6,4) * t113;
t136 = Ifges(6,4) * t116;
t77 = Ifges(6,5) * t113 + Ifges(6,6) * t116;
t78 = Ifges(5,5) * t113 + Ifges(5,6) * t116;
t131 = Ifges(4,5) * t114 + Ifges(4,6) * t117;
t130 = t113 ^ 2 + t116 ^ 2;
t8 = Ifges(6,4) * t36 + Ifges(6,2) * t35 + Ifges(6,6) * t58;
t9 = Ifges(5,4) * t36 + Ifges(5,2) * t35 + Ifges(5,6) * t58;
t129 = -t8 / 0.2e1 - t9 / 0.2e1;
t128 = Ifges(3,5) * t135 + Ifges(3,6) * t134 + Ifges(3,3) * t112;
t10 = Ifges(6,1) * t36 + Ifges(6,4) * t35 + Ifges(6,5) * t58;
t11 = Ifges(5,1) * t36 + Ifges(5,4) * t35 + Ifges(5,5) * t58;
t127 = t10 / 0.2e1 + t11 / 0.2e1;
t50 = -Ifges(6,6) * t117 + (-Ifges(6,2) * t113 + t136) * t114;
t51 = -Ifges(5,6) * t117 + (-Ifges(5,2) * t113 + t138) * t114;
t126 = t50 / 0.2e1 + t51 / 0.2e1;
t52 = -Ifges(6,5) * t117 + (Ifges(6,1) * t116 - t137) * t114;
t53 = -Ifges(5,5) * t117 + (Ifges(5,1) * t116 - t139) * t114;
t125 = t52 / 0.2e1 + t53 / 0.2e1;
t124 = t77 / 0.2e1 + t78 / 0.2e1;
t79 = Ifges(6,2) * t116 + t137;
t80 = Ifges(5,2) * t116 + t139;
t123 = t79 / 0.2e1 + t80 / 0.2e1;
t82 = Ifges(6,1) * t113 + t136;
t83 = Ifges(5,1) * t113 + t138;
t122 = t82 / 0.2e1 + t83 / 0.2e1;
t12 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t73 = -t116 * mrSges(6,1) + t113 * mrSges(6,2);
t3 = -t113 * t19 + t116 * t17;
t24 = -t114 * t46 + t117 * t47;
t18 = pkin(3) * t134 - t24;
t121 = mrSges(5,1) * t113 + mrSges(5,2) * t116;
t120 = pkin(8) ^ 2;
t110 = t117 ^ 2;
t108 = t114 ^ 2;
t106 = t108 * t120;
t97 = -pkin(4) * t116 - pkin(3);
t94 = Ifges(5,5) * t132;
t93 = Ifges(6,5) * t132;
t84 = Ifges(4,1) * t114 + Ifges(4,4) * t117;
t81 = Ifges(4,4) * t114 + Ifges(4,2) * t117;
t76 = t141 * t116;
t75 = -mrSges(4,1) * t117 + mrSges(4,2) * t114;
t72 = t141 * t113;
t70 = (pkin(4) * t113 + pkin(8)) * t114;
t69 = -mrSges(5,1) * t117 - mrSges(5,3) * t132;
t68 = -mrSges(6,1) * t117 - mrSges(6,3) * t132;
t67 = mrSges(5,2) * t117 - mrSges(5,3) * t133;
t66 = mrSges(6,2) * t117 - mrSges(6,3) * t133;
t65 = t116 * t71;
t63 = t121 * t114;
t49 = -Ifges(5,6) * t133 - Ifges(5,3) * t117 + t94;
t48 = -Ifges(6,6) * t133 - Ifges(6,3) * t117 + t93;
t43 = -t113 * t145 + t65;
t39 = -mrSges(4,1) * t134 - mrSges(4,3) * t59;
t38 = mrSges(4,2) * t134 - mrSges(4,3) * t58;
t37 = -qJ(5) * t133 + t44;
t29 = -qJ(5) * t132 + t65 + (-pkin(8) * t113 - pkin(4)) * t117;
t28 = mrSges(4,1) * t58 + mrSges(4,2) * t59;
t27 = Ifges(4,1) * t59 - Ifges(4,4) * t58 - Ifges(4,5) * t134;
t26 = Ifges(4,4) * t59 - Ifges(4,2) * t58 - Ifges(4,6) * t134;
t23 = mrSges(5,1) * t58 - mrSges(5,3) * t36;
t22 = mrSges(6,1) * t58 - mrSges(6,3) * t36;
t21 = -mrSges(5,2) * t58 + mrSges(5,3) * t35;
t20 = -mrSges(6,2) * t58 + mrSges(6,3) * t35;
t13 = -mrSges(5,1) * t35 + mrSges(5,2) * t36;
t5 = -pkin(4) * t35 + t18;
t2 = qJ(5) * t35 + t4;
t1 = pkin(4) * t58 - qJ(5) * t36 + t3;
t14 = [0.2e1 * t1 * t22 + 0.2e1 * t5 * t12 + 0.2e1 * t18 * t13 + 0.2e1 * t2 * t20 + 0.2e1 * t4 * t21 + 0.2e1 * t3 * t23 + 0.2e1 * t24 * t39 + 0.2e1 * t25 * t38 + t59 * t27 + 0.2e1 * t45 * t28 + Ifges(2,3) + (t10 + t11) * t36 + (t8 + t9) * t35 + (-t26 + t149) * t58 + (t128 - 0.2e1 * t143 + 0.2e1 * t144) * t112 + ((-0.2e1 * t60 * mrSges(3,3) + Ifges(3,5) * t112 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t115) * t111) * t115 + (0.2e1 * t61 * mrSges(3,3) + Ifges(3,6) * t112 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t115 + (Ifges(3,2) + Ifges(4,3)) * t118) * t111 + t140) * t118) * t111 + m(3) * (pkin(1) ^ 2 * t111 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2 + t45 ^ 2) + m(5) * (t18 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); t128 - t131 * t134 / 0.2e1 + t125 * t36 + t126 * t35 + (t27 / 0.2e1 - t24 * mrSges(4,3) + t127 * t116 + t129 * t113 + (-t39 + t13) * pkin(8)) * t114 + t45 * t75 + t59 * t84 / 0.2e1 + t5 * t62 + t18 * t63 + t2 * t66 + t4 * t67 + t1 * t68 + t3 * t69 + t70 * t12 + t37 * t20 + t43 * t23 + t44 * t21 + t29 * t22 + m(4) * (-pkin(2) * t45 + (-t24 * t114 + t25 * t117) * pkin(8)) - pkin(2) * t28 - t143 + t144 + (-t81 / 0.2e1 + t48 / 0.2e1 + t49 / 0.2e1) * t58 + m(5) * (pkin(8) * t114 * t18 + t3 * t43 + t4 * t44) + m(6) * (t1 * t29 + t2 * t37 + t5 * t70) + (t26 / 0.2e1 - t6 / 0.2e1 - t7 / 0.2e1 + pkin(8) * t38 + t25 * mrSges(4,3)) * t117; -0.2e1 * pkin(2) * t75 + 0.2e1 * t29 * t68 + 0.2e1 * t37 * t66 + 0.2e1 * t43 * t69 + 0.2e1 * t44 * t67 + 0.2e1 * t70 * t62 + Ifges(3,3) + (t108 + t110) * mrSges(4,3) * t151 + (-t49 - t48 + t81) * t117 + m(5) * (t43 ^ 2 + t44 ^ 2 + t106) + m(6) * (t29 ^ 2 + t37 ^ 2 + t70 ^ 2) + m(4) * (pkin(2) ^ 2 + t110 * t120 + t106) + (t63 * t151 + t84 + (t52 + t53) * t116 + (-t50 - t51) * t113) * t114; -Ifges(4,3) * t134 + t24 * mrSges(4,1) - t25 * mrSges(4,2) - pkin(3) * t13 + t97 * t12 - t76 * t20 + t72 * t22 + t5 * t73 + t124 * t58 + t122 * t36 + t123 * t35 + m(6) * (t1 * t72 - t2 * t76 + t5 * t97) + (t2 * mrSges(6,3) + t4 * mrSges(5,3) + (m(5) * t4 + t21) * pkin(9) - t129) * t116 + (-t1 * mrSges(6,3) - t3 * mrSges(5,3) + (-m(5) * t3 - t23) * pkin(9) + t127) * t113 - t140 + t148 * t18; -t76 * t66 + t97 * t62 - pkin(3) * t63 + t72 * t68 + t70 * t73 + m(6) * (t29 * t72 - t37 * t76 + t70 * t97) - t124 * t117 + (-t117 * mrSges(4,2) + (-mrSges(4,1) + t148) * t114) * pkin(8) + (t44 * mrSges(5,3) + t37 * mrSges(6,3) + t122 * t114 + (m(5) * t44 + t67) * pkin(9) + t126) * t116 + (-t43 * mrSges(5,3) - t29 * mrSges(6,3) - t123 * t114 + (-m(5) * t43 - t69) * pkin(9) + t125) * t113 + t131; -0.2e1 * pkin(3) * t74 + 0.2e1 * t97 * t73 + Ifges(4,3) + 0.2e1 * t130 * pkin(9) * mrSges(5,3) + m(6) * (t72 ^ 2 + t76 ^ 2 + t97 ^ 2) + m(5) * (pkin(9) ^ 2 * t130 + pkin(3) ^ 2) + (t76 * t150 + t79 + t80) * t116 + (t72 * t150 + t82 + t83) * t113; mrSges(5,1) * t3 + mrSges(6,1) * t1 - mrSges(5,2) * t4 - mrSges(6,2) * t2 + (m(6) * t1 + t22) * pkin(4) + t149; mrSges(5,1) * t43 + mrSges(6,1) * t29 - mrSges(5,2) * t44 - mrSges(6,2) * t37 + t93 + t94 + t142 * t117 + (-Ifges(5,6) - Ifges(6,6)) * t133 + (m(6) * t29 + t68) * pkin(4); mrSges(6,1) * t72 + mrSges(6,2) * t76 - t121 * pkin(9) + (m(6) * t72 - t113 * mrSges(6,3)) * pkin(4) + t78 + t77; (m(6) * pkin(4) + 0.2e1 * mrSges(6,1)) * pkin(4) - t142; m(6) * t5 + t12; m(6) * t70 + t62; m(6) * t97 + t73; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
