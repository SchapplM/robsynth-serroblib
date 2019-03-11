% Calculate joint inertia matrix for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:01
% EndTime: 2019-03-09 12:22:03
% DurationCPUTime: 1.14s
% Computational Cost: add. (2270->309), mult. (4585->418), div. (0->0), fcn. (4846->8), ass. (0->110)
t150 = 2 * pkin(7);
t120 = sin(qJ(2));
t116 = sin(pkin(10));
t117 = cos(pkin(10));
t119 = sin(qJ(4));
t122 = cos(qJ(4));
t93 = t116 * t122 + t117 * t119;
t81 = t93 * t120;
t92 = -t116 * t119 + t117 * t122;
t82 = t92 * t120;
t149 = -Ifges(5,5) * t82 + Ifges(5,6) * t81;
t148 = 2 * mrSges(7,3);
t147 = -t116 / 0.2e1;
t146 = t117 / 0.2e1;
t123 = cos(qJ(2));
t145 = pkin(7) * t123;
t111 = t120 * pkin(7);
t143 = Ifges(7,2) + Ifges(6,3);
t142 = pkin(8) + qJ(3);
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t138 = t117 * t120;
t97 = -pkin(2) * t123 - qJ(3) * t120 - pkin(1);
t89 = t117 * t97;
t61 = -pkin(8) * t138 + t89 + (-pkin(7) * t116 - pkin(3)) * t123;
t139 = t116 * t120;
t73 = t116 * t97 + t117 * t145;
t67 = -pkin(8) * t139 + t73;
t31 = -t119 * t67 + t122 * t61;
t15 = -pkin(4) * t123 - pkin(9) * t82 + t31;
t32 = t119 * t61 + t122 * t67;
t24 = -pkin(9) * t81 + t32;
t6 = t118 * t15 + t121 * t24;
t133 = t142 * t117;
t134 = t142 * t116;
t69 = -t119 * t134 + t122 * t133;
t141 = Ifges(4,4) * t116;
t140 = Ifges(4,4) * t117;
t48 = -t118 * t81 + t121 * t82;
t36 = t123 * mrSges(7,1) + t48 * mrSges(7,2);
t84 = mrSges(4,1) * t139 + mrSges(4,2) * t138;
t96 = pkin(3) * t139 + t111;
t137 = t116 ^ 2 + t117 ^ 2;
t68 = -t119 * t133 - t122 * t134;
t127 = -t93 * pkin(9) + t68;
t49 = pkin(9) * t92 + t69;
t21 = t118 * t49 - t121 * t127;
t23 = t118 * t127 + t121 * t49;
t136 = t21 ^ 2 + t23 ^ 2;
t135 = Ifges(5,3) + t143;
t106 = -pkin(3) * t117 - pkin(2);
t50 = t81 * mrSges(5,1) + t82 * mrSges(5,2);
t64 = -t92 * mrSges(5,1) + t93 * mrSges(5,2);
t47 = t118 * t82 + t121 * t81;
t17 = t47 * mrSges(6,1) + t48 * mrSges(6,2);
t59 = t118 * t93 - t121 * t92;
t60 = t118 * t92 + t121 * t93;
t26 = t59 * mrSges(6,1) + t60 * mrSges(6,2);
t16 = t47 * mrSges(7,1) - t48 * mrSges(7,3);
t25 = t59 * mrSges(7,1) - t60 * mrSges(7,3);
t132 = (-Ifges(7,4) - Ifges(6,5)) * t48 + (Ifges(6,6) - Ifges(7,6)) * t47;
t98 = -t117 * mrSges(4,1) + t116 * mrSges(4,2);
t63 = pkin(4) * t81 + t96;
t72 = -t116 * t145 + t89;
t131 = -t72 * t116 + t73 * t117;
t5 = -t118 * t24 + t121 * t15;
t75 = -pkin(4) * t92 + t106;
t130 = (mrSges(6,1) * t121 - mrSges(6,2) * t118) * pkin(4);
t2 = -qJ(6) * t123 + t6;
t3 = pkin(5) * t123 - t5;
t129 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) - t132;
t55 = Ifges(7,6) * t59;
t56 = Ifges(6,6) * t59;
t57 = Ifges(6,5) * t60;
t58 = Ifges(7,4) * t60;
t128 = t55 - t56 + t57 + t58 + (-mrSges(6,2) + mrSges(7,3)) * t23 + (-mrSges(6,1) - mrSges(7,1)) * t21;
t125 = pkin(7) ^ 2;
t115 = t123 ^ 2;
t114 = t120 ^ 2;
t110 = t114 * t125;
t107 = -pkin(4) * t121 - pkin(5);
t105 = pkin(4) * t118 + qJ(6);
t100 = Ifges(4,1) * t116 + t140;
t99 = Ifges(4,2) * t117 + t141;
t95 = -mrSges(4,1) * t123 - mrSges(4,3) * t138;
t94 = mrSges(4,2) * t123 - mrSges(4,3) * t139;
t87 = Ifges(5,5) * t93;
t86 = Ifges(5,6) * t92;
t80 = -t123 * Ifges(4,5) + (Ifges(4,1) * t117 - t141) * t120;
t79 = -t123 * Ifges(4,6) + (-Ifges(4,2) * t116 + t140) * t120;
t71 = -mrSges(5,1) * t123 - mrSges(5,3) * t82;
t70 = mrSges(5,2) * t123 - mrSges(5,3) * t81;
t66 = Ifges(5,1) * t93 + Ifges(5,4) * t92;
t65 = Ifges(5,4) * t93 + Ifges(5,2) * t92;
t46 = Ifges(5,1) * t82 - Ifges(5,4) * t81 - Ifges(5,5) * t123;
t45 = Ifges(5,4) * t82 - Ifges(5,2) * t81 - Ifges(5,6) * t123;
t35 = -mrSges(6,1) * t123 - mrSges(6,3) * t48;
t34 = mrSges(6,2) * t123 - mrSges(6,3) * t47;
t33 = -mrSges(7,2) * t47 - mrSges(7,3) * t123;
t30 = Ifges(6,1) * t60 - Ifges(6,4) * t59;
t29 = Ifges(7,1) * t60 + Ifges(7,5) * t59;
t28 = Ifges(6,4) * t60 - Ifges(6,2) * t59;
t27 = Ifges(7,5) * t60 + Ifges(7,3) * t59;
t20 = pkin(5) * t59 - qJ(6) * t60 + t75;
t11 = Ifges(6,1) * t48 - Ifges(6,4) * t47 - Ifges(6,5) * t123;
t10 = Ifges(7,1) * t48 - Ifges(7,4) * t123 + Ifges(7,5) * t47;
t9 = Ifges(6,4) * t48 - Ifges(6,2) * t47 - Ifges(6,6) * t123;
t8 = Ifges(7,5) * t48 - Ifges(7,6) * t123 + Ifges(7,3) * t47;
t7 = pkin(5) * t47 - qJ(6) * t48 + t63;
t1 = [0.2e1 * t7 * t16 + 0.2e1 * t63 * t17 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t36 + 0.2e1 * t31 * t71 + 0.2e1 * t32 * t70 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 - t81 * t45 + t82 * t46 + 0.2e1 * t96 * t50 + 0.2e1 * t72 * t95 + 0.2e1 * t73 * t94 + Ifges(2,3) + (t10 + t11) * t48 + (t8 - t9) * t47 + (t114 + t115) * mrSges(3,3) * t150 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t120 - t116 * t79 + t117 * t80 + t84 * t150) * t120 + m(3) * (pkin(1) ^ 2 + t115 * t125 + t110) + m(4) * (t72 ^ 2 + t73 ^ 2 + t110) + m(5) * (t31 ^ 2 + t32 ^ 2 + t96 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t63 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2) + t135) * t123 + (-Ifges(4,5) * t117 + Ifges(4,6) * t116 + (2 * Ifges(3,4))) * t120 + t132 + t149) * t123; (t10 / 0.2e1 + t11 / 0.2e1 + t3 * mrSges(7,2) - t5 * mrSges(6,3)) * t60 + (t8 / 0.2e1 - t9 / 0.2e1 - t2 * mrSges(7,2) - t6 * mrSges(6,3)) * t59 + (t29 / 0.2e1 + t30 / 0.2e1) * t48 + (t27 / 0.2e1 - t28 / 0.2e1) * t47 + (-t31 * t93 + t32 * t92) * mrSges(5,3) + (t100 * t146 + t99 * t147 + Ifges(3,5) + (t98 - mrSges(3,1)) * pkin(7)) * t120 + (-pkin(7) * mrSges(3,2) + Ifges(4,5) * t147 - Ifges(4,6) * t117 / 0.2e1 - t87 / 0.2e1 - t86 / 0.2e1 - t57 / 0.2e1 + t56 / 0.2e1 - t58 / 0.2e1 - t55 / 0.2e1 + Ifges(3,6)) * t123 + m(4) * (-pkin(2) * t111 + t131 * qJ(3)) + t131 * mrSges(4,3) + t116 * t80 / 0.2e1 + t106 * t50 + t92 * t45 / 0.2e1 + t93 * t46 / 0.2e1 + t96 * t64 - t81 * t65 / 0.2e1 + t82 * t66 / 0.2e1 - pkin(2) * t84 + t63 * t26 + t69 * t70 + t68 * t71 + t75 * t17 + t7 * t25 + t20 * t16 + (t33 + t34) * t23 + (-t35 + t36) * t21 + (-t116 * t95 + t117 * t94) * qJ(3) + m(5) * (t106 * t96 + t31 * t68 + t32 * t69) + m(6) * (-t21 * t5 + t23 * t6 + t63 * t75) + m(7) * (t2 * t23 + t20 * t7 + t21 * t3) + t79 * t146; -0.2e1 * pkin(2) * t98 + t116 * t100 + 0.2e1 * t106 * t64 + t117 * t99 + 0.2e1 * t20 * t25 + 0.2e1 * t75 * t26 + t92 * t65 + t93 * t66 + Ifges(3,3) + (t29 + t30) * t60 + (t27 - t28) * t59 + m(4) * (t137 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t106 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t75 ^ 2 + t136) + m(7) * (t20 ^ 2 + t136) + 0.2e1 * (-t68 * t93 + t69 * t92) * mrSges(5,3) + 0.2e1 * t137 * qJ(3) * mrSges(4,3) + 0.2e1 * (t21 * t60 - t23 * t59) * (mrSges(7,2) + mrSges(6,3)); m(4) * t111 + m(5) * t96 + m(6) * t63 + m(7) * t7 + t16 + t17 + t50 + t84; -m(4) * pkin(2) + m(5) * t106 + m(6) * t75 + m(7) * t20 + t25 + t26 + t64 + t98; m(4) + m(5) + m(6) + m(7); m(7) * (t105 * t2 + t107 * t3) - t135 * t123 + t105 * t33 + t107 * t36 + t31 * mrSges(5,1) - t32 * mrSges(5,2) + t129 + (t121 * t35 + t118 * t34 + m(6) * (t118 * t6 + t121 * t5)) * pkin(4) - t149; m(7) * (t105 * t23 + t107 * t21) + t87 + t86 - t69 * mrSges(5,2) + t68 * mrSges(5,1) + (-t105 * t59 + t107 * t60) * mrSges(7,2) + (m(6) * (t118 * t23 - t121 * t21) + (-t118 * t59 - t121 * t60) * mrSges(6,3)) * pkin(4) + t128; 0; -0.2e1 * t107 * mrSges(7,1) + t105 * t148 + 0.2e1 * t130 + m(7) * (t105 ^ 2 + t107 ^ 2) + m(6) * (t118 ^ 2 + t121 ^ 2) * pkin(4) ^ 2 + t135; -pkin(5) * t36 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) + qJ(6) * t33 - t143 * t123 + t129; m(7) * (-pkin(5) * t21 + qJ(6) * t23) + (-pkin(5) * t60 - qJ(6) * t59) * mrSges(7,2) + t128; 0; m(7) * (-pkin(5) * t107 + qJ(6) * t105) + t130 + (qJ(6) + t105) * mrSges(7,3) + (pkin(5) - t107) * mrSges(7,1) + t143; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t148 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t143; m(7) * t3 + t36; m(7) * t21 + t60 * mrSges(7,2); 0; m(7) * t107 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
