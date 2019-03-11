% Calculate joint inertia matrix for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:18
% EndTime: 2019-03-09 17:14:21
% DurationCPUTime: 1.21s
% Computational Cost: add. (1329->324), mult. (2503->423), div. (0->0), fcn. (2283->6), ass. (0->116)
t135 = Ifges(6,3) + Ifges(7,3);
t109 = sin(qJ(3));
t112 = cos(qJ(3));
t149 = t109 ^ 2 + t112 ^ 2;
t148 = 2 * pkin(7);
t137 = mrSges(6,2) + mrSges(7,2);
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t114 = -pkin(3) - pkin(4);
t71 = qJ(4) * t111 + t108 * t114;
t147 = t137 * t71;
t146 = t137 * t108;
t110 = sin(qJ(2));
t126 = t110 * t112;
t127 = t109 * t110;
t145 = -Ifges(5,6) * t127 + (-Ifges(5,4) - Ifges(4,5)) * t126;
t144 = m(7) * pkin(5);
t143 = pkin(8) - pkin(9);
t142 = pkin(3) * t109;
t113 = cos(qJ(2));
t141 = pkin(7) * t113;
t119 = t108 * t109 + t111 * t112;
t140 = t119 * t71;
t70 = -qJ(4) * t108 + t111 * t114;
t139 = t70 * mrSges(6,1);
t138 = -mrSges(6,1) - mrSges(7,1);
t136 = Ifges(5,2) + Ifges(4,3);
t103 = t113 * pkin(3);
t73 = -pkin(2) * t113 - pkin(8) * t110 - pkin(1);
t88 = t109 * t141;
t16 = pkin(4) * t113 + t103 + t88 + (-pkin(9) * t110 - t73) * t112;
t36 = t109 * t73 + t112 * t141;
t29 = -qJ(4) * t113 + t36;
t18 = pkin(9) * t127 + t29;
t4 = t108 * t16 + t111 * t18;
t60 = -t108 * t112 + t109 * t111;
t47 = t60 * t110;
t31 = -mrSges(7,2) * t113 + mrSges(7,3) * t47;
t32 = -mrSges(6,2) * t113 + mrSges(6,3) * t47;
t134 = t31 + t32;
t80 = t143 * t109;
t81 = t143 * t112;
t28 = t108 * t80 + t111 * t81;
t66 = t113 * mrSges(5,1) + mrSges(5,2) * t126;
t133 = Ifges(4,4) * t109;
t132 = Ifges(4,4) * t112;
t131 = Ifges(5,5) * t109;
t130 = Ifges(5,5) * t112;
t129 = Ifges(5,6) * t113;
t128 = t149 * pkin(8) ^ 2;
t72 = -t112 * pkin(3) - t109 * qJ(4) - pkin(2);
t48 = t119 * t110;
t12 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t19 = mrSges(7,1) * t119 + t60 * mrSges(7,2);
t3 = -t108 * t18 + t111 * t16;
t27 = -t108 * t81 + t111 * t80;
t35 = t112 * t73 - t88;
t51 = t112 * pkin(4) - t72;
t124 = t109 * mrSges(4,1) + t112 * mrSges(4,2);
t123 = t109 * mrSges(5,1) - t112 * mrSges(5,3);
t122 = qJ(4) * t112 - t142;
t120 = (Ifges(6,5) + Ifges(7,5)) * t48 + (Ifges(6,6) + Ifges(7,6)) * t47 + t135 * t113;
t83 = qJ(4) * t126;
t26 = t83 + (t114 * t109 - pkin(7)) * t110;
t10 = -qJ(6) * t60 + t27;
t11 = -qJ(6) * t119 + t28;
t54 = Ifges(7,6) * t119;
t55 = Ifges(6,6) * t119;
t56 = Ifges(7,5) * t60;
t57 = Ifges(6,5) * t60;
t118 = t27 * mrSges(6,1) + t10 * mrSges(7,1) - t28 * mrSges(6,2) - t11 * mrSges(7,2) - t54 - t55 + t56 + t57;
t1 = pkin(5) * t113 - qJ(6) * t48 + t3;
t2 = qJ(6) * t47 + t4;
t117 = t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + t120;
t116 = pkin(7) ^ 2;
t107 = t113 ^ 2;
t105 = t110 ^ 2;
t99 = t105 * t116;
t97 = Ifges(5,4) * t109;
t96 = Ifges(4,5) * t109;
t95 = Ifges(4,6) * t112;
t79 = Ifges(4,1) * t109 + t132;
t78 = Ifges(5,1) * t109 - t130;
t77 = Ifges(4,2) * t112 + t133;
t76 = -Ifges(5,3) * t112 + t131;
t75 = -mrSges(4,1) * t112 + mrSges(4,2) * t109;
t74 = -mrSges(5,1) * t112 - mrSges(5,3) * t109;
t69 = t71 ^ 2;
t68 = -pkin(5) + t70;
t67 = -mrSges(5,2) * t127 - mrSges(5,3) * t113;
t65 = -mrSges(4,1) * t113 - mrSges(4,3) * t126;
t64 = mrSges(4,2) * t113 - mrSges(4,3) * t127;
t52 = t108 * t71;
t50 = t124 * t110;
t49 = t123 * t110;
t46 = -t83 + (pkin(7) + t142) * t110;
t45 = -Ifges(4,5) * t113 + (Ifges(4,1) * t112 - t133) * t110;
t44 = -Ifges(5,4) * t113 + (Ifges(5,1) * t112 + t131) * t110;
t43 = -Ifges(4,6) * t113 + (-Ifges(4,2) * t109 + t132) * t110;
t42 = -t129 + (Ifges(5,3) * t109 + t130) * t110;
t34 = mrSges(6,1) * t113 - mrSges(6,3) * t48;
t33 = mrSges(7,1) * t113 - mrSges(7,3) * t48;
t30 = t103 - t35;
t25 = pkin(5) * t119 + t51;
t24 = Ifges(6,1) * t60 - Ifges(6,4) * t119;
t23 = Ifges(7,1) * t60 - Ifges(7,4) * t119;
t22 = Ifges(6,4) * t60 - Ifges(6,2) * t119;
t21 = Ifges(7,4) * t60 - Ifges(7,2) * t119;
t20 = mrSges(6,1) * t119 + mrSges(6,2) * t60;
t13 = -mrSges(6,1) * t47 + mrSges(6,2) * t48;
t9 = Ifges(6,1) * t48 + Ifges(6,4) * t47 + Ifges(6,5) * t113;
t8 = Ifges(7,1) * t48 + Ifges(7,4) * t47 + Ifges(7,5) * t113;
t7 = Ifges(6,4) * t48 + Ifges(6,2) * t47 + Ifges(6,6) * t113;
t6 = Ifges(7,4) * t48 + Ifges(7,2) * t47 + Ifges(7,6) * t113;
t5 = -pkin(5) * t47 + t26;
t14 = [0.2e1 * t1 * t33 + 0.2e1 * t5 * t12 + 0.2e1 * t26 * t13 + 0.2e1 * t2 * t31 + 0.2e1 * t29 * t67 + 0.2e1 * t3 * t34 + 0.2e1 * t30 * t66 + 0.2e1 * t4 * t32 + 0.2e1 * t35 * t65 + 0.2e1 * t36 * t64 + 0.2e1 * t46 * t49 + Ifges(2,3) + (t8 + t9) * t48 + (t6 + t7) * t47 + (t105 + t107) * mrSges(3,3) * t148 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t110 + t50 * t148 + (t44 + t45) * t112 + (t42 - t43) * t109) * t110 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t136) * t113 + (Ifges(4,6) * t109 + (2 * Ifges(3,4))) * t110 + t120 + t145) * t113 + m(3) * (pkin(1) ^ 2 + t107 * t116 + t99) + m(4) * (t35 ^ 2 + t36 ^ 2 + t99) + m(5) * (t29 ^ 2 + t30 ^ 2 + t46 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t26 ^ 2 + t3 ^ 2 + t4 ^ 2); Ifges(3,5) * t110 - pkin(2) * t50 + t10 * t33 + t11 * t31 + t25 * t12 + t51 * t13 + t5 * t19 + t26 * t20 + t27 * t34 + t28 * t32 + t46 * t74 + t72 * t49 + (t23 / 0.2e1 + t24 / 0.2e1) * t48 + (t21 / 0.2e1 + t22 / 0.2e1) * t47 + (-t97 / 0.2e1 - t96 / 0.2e1 - t95 / 0.2e1 + t56 / 0.2e1 - t54 / 0.2e1 + t57 / 0.2e1 - t55 / 0.2e1 + Ifges(3,6)) * t113 + (-t113 * mrSges(3,2) + (-mrSges(3,1) + t75) * t110) * pkin(7) + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + t8 / 0.2e1 + t9 / 0.2e1) * t60 - (t2 * mrSges(7,3) + t4 * mrSges(6,3) + t6 / 0.2e1 + t7 / 0.2e1) * t119 + (t129 / 0.2e1 - t42 / 0.2e1 + t43 / 0.2e1 + t29 * mrSges(5,2) + t36 * mrSges(4,3) + (t78 / 0.2e1 + t79 / 0.2e1) * t110 + (t64 + t67) * pkin(8)) * t112 + (t44 / 0.2e1 + t45 / 0.2e1 + t30 * mrSges(5,2) - t35 * mrSges(4,3) + (t76 / 0.2e1 - t77 / 0.2e1) * t110 + (-t65 + t66) * pkin(8)) * t109 + m(4) * (-pkin(2) * pkin(7) * t110 + (-t35 * t109 + t36 * t112) * pkin(8)) + m(5) * (t46 * t72 + (t30 * t109 + t29 * t112) * pkin(8)) + m(6) * (t26 * t51 + t27 * t3 + t28 * t4) + m(7) * (t1 * t10 + t11 * t2 + t25 * t5); -0.2e1 * pkin(2) * t75 + 0.2e1 * t25 * t19 + 0.2e1 * t51 * t20 + 0.2e1 * t72 * t74 + Ifges(3,3) + (-t76 + t77) * t112 + (t79 + t78) * t109 + (-0.2e1 * mrSges(6,3) * t27 - 0.2e1 * mrSges(7,3) * t10 + t23 + t24) * t60 - (0.2e1 * mrSges(6,3) * t28 + 0.2e1 * mrSges(7,3) * t11 + t21 + t22) * t119 + m(5) * (t72 ^ 2 + t128) + m(4) * (pkin(2) ^ 2 + t128) + m(6) * (t27 ^ 2 + t28 ^ 2 + t51 ^ 2) + m(7) * (t10 ^ 2 + t11 ^ 2 + t25 ^ 2) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(8) * t149; -t117 + t134 * t71 + m(7) * (t1 * t68 + t2 * t71) + m(6) * (t3 * t70 + t4 * t71) + m(5) * (-pkin(3) * t30 + qJ(4) * t29) - t136 * t113 - Ifges(4,6) * t127 - pkin(3) * t66 + qJ(4) * t67 + t68 * t33 + t70 * t34 + t29 * mrSges(5,3) - t30 * mrSges(5,1) + t35 * mrSges(4,1) - t36 * mrSges(4,2) - t145; -t118 + (m(5) * t122 - t123 - t124) * pkin(8) + m(7) * (t10 * t68 + t11 * t71) + m(6) * (t27 * t70 + t28 * t71) + (-t60 * t68 - t140) * mrSges(7,3) + (-t60 * t70 - t140) * mrSges(6,3) + t122 * mrSges(5,2) - Ifges(5,6) * t112 + t95 + t96 + t97; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t139 - 0.2e1 * t68 * mrSges(7,1) + 0.2e1 * qJ(4) * mrSges(5,3) + 0.2e1 * t147 + m(6) * (t70 ^ 2 + t69) + m(7) * (t68 ^ 2 + t69) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t135 + t136; (t33 + t34) * t111 + t134 * t108 + m(7) * (t1 * t111 + t108 * t2) + m(6) * (t108 * t4 + t111 * t3) + m(5) * t30 + t66; (m(5) * pkin(8) + mrSges(5,2)) * t109 + m(7) * (t10 * t111 + t108 * t11) + m(6) * (t108 * t28 + t111 * t27) + (mrSges(7,3) + mrSges(6,3)) * (-t108 * t119 - t111 * t60); -m(5) * pkin(3) - mrSges(5,1) + t138 * t111 + t146 + m(6) * (t111 * t70 + t52) + m(7) * (t111 * t68 + t52); m(5) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t108 ^ 2 + t111 ^ 2); (m(7) * t1 + t33) * pkin(5) + t117; (m(7) * t10 - t60 * mrSges(7,3)) * pkin(5) + t118; t68 * t144 + t139 - t147 + (-pkin(5) + t68) * mrSges(7,1) - t135; -t146 + (-t138 + t144) * t111; (0.2e1 * mrSges(7,1) + t144) * pkin(5) + t135; m(7) * t5 + t12; m(7) * t25 + t19; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
