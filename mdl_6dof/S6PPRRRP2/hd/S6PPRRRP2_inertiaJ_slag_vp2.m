% Calculate joint inertia matrix for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:06
% EndTime: 2019-03-08 18:56:08
% DurationCPUTime: 0.79s
% Computational Cost: add. (833->228), mult. (2186->319), div. (0->0), fcn. (2311->12), ass. (0->100)
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t137 = t80 ^ 2 + t83 ^ 2;
t75 = sin(pkin(7));
t85 = cos(qJ(3));
t113 = t75 * t85;
t82 = sin(qJ(3));
t114 = t75 * t82;
t78 = cos(pkin(7));
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t35 = t84 * t114 + t81 * t78;
t18 = t83 * t113 + t80 * t35;
t20 = -t80 * t113 + t83 * t35;
t136 = t18 * t80 + t20 * t83;
t77 = cos(pkin(12));
t112 = t77 * t78;
t74 = sin(pkin(12));
t76 = sin(pkin(6));
t79 = cos(pkin(6));
t16 = t79 * t114 + (t82 * t112 + t74 * t85) * t76;
t31 = -t76 * t77 * t75 + t79 * t78;
t11 = t16 * t84 + t31 * t81;
t14 = -t79 * t113 + (-t112 * t85 + t74 * t82) * t76;
t4 = t11 * t80 - t14 * t83;
t6 = t11 * t83 + t14 * t80;
t135 = t4 * t80 + t6 * t83;
t134 = 2 * pkin(9);
t133 = m(6) + m(7);
t132 = mrSges(6,3) + mrSges(7,2);
t131 = -m(7) * pkin(5) - mrSges(7,1);
t130 = -mrSges(6,1) + t131;
t129 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t9 = t16 * t81 - t31 * t84;
t8 = t9 ^ 2;
t128 = t14 ^ 2;
t33 = t81 * t114 - t84 * t78;
t32 = t33 ^ 2;
t127 = pkin(9) * t84;
t7 = t33 * t9;
t124 = t9 * t81;
t123 = Ifges(6,4) * t80;
t122 = Ifges(6,4) * t83;
t121 = Ifges(7,5) * t80;
t120 = Ifges(7,5) * t83;
t119 = Ifges(6,6) * t84;
t118 = Ifges(7,6) * t84;
t115 = t33 * t81;
t111 = t80 * t81;
t110 = t81 * t83;
t45 = -t84 * pkin(4) - t81 * pkin(10) - pkin(3);
t109 = t83 * t45;
t108 = t84 * mrSges(5,3);
t48 = -t84 * mrSges(5,1) + t81 * mrSges(5,2);
t107 = mrSges(4,1) - t48;
t106 = Ifges(7,2) + Ifges(6,3);
t40 = t84 * mrSges(6,2) - mrSges(6,3) * t111;
t43 = -mrSges(7,2) * t111 - t84 * mrSges(7,3);
t105 = t40 + t43;
t41 = -t84 * mrSges(6,1) - mrSges(6,3) * t110;
t42 = t84 * mrSges(7,1) + mrSges(7,2) * t110;
t104 = -t41 + t42;
t47 = -t83 * mrSges(6,1) + t80 * mrSges(6,2);
t103 = t47 - mrSges(5,1);
t24 = t83 * t127 + t80 * t45;
t102 = t137 * pkin(10) ^ 2;
t46 = -t83 * mrSges(7,1) - t80 * mrSges(7,3);
t101 = t46 + t103;
t100 = -Ifges(7,6) * t111 + (-Ifges(7,4) - Ifges(6,5)) * t110;
t97 = t135 * pkin(10);
t94 = t136 * pkin(10);
t90 = t80 * mrSges(7,1) - t83 * mrSges(7,3);
t36 = t90 * t81;
t91 = t80 * mrSges(6,1) + t83 * mrSges(6,2);
t37 = t91 * t81;
t93 = t81 * mrSges(5,3) + t36 + t37;
t89 = -pkin(5) * t80 + qJ(6) * t83;
t87 = pkin(9) ^ 2;
t73 = t84 ^ 2;
t71 = t81 ^ 2;
t68 = t75 ^ 2;
t66 = t71 * t87;
t64 = Ifges(7,4) * t80;
t63 = Ifges(6,5) * t80;
t62 = Ifges(6,6) * t83;
t59 = t68 * t85 ^ 2;
t52 = Ifges(6,1) * t80 + t122;
t51 = Ifges(7,1) * t80 - t120;
t50 = Ifges(6,2) * t83 + t123;
t49 = -Ifges(7,3) * t83 + t121;
t44 = -t83 * pkin(5) - t80 * qJ(6) - pkin(4);
t30 = (pkin(9) - t89) * t81;
t29 = -Ifges(6,5) * t84 + (Ifges(6,1) * t83 - t123) * t81;
t28 = -Ifges(7,4) * t84 + (Ifges(7,1) * t83 + t121) * t81;
t27 = -t119 + (-Ifges(6,2) * t80 + t122) * t81;
t26 = -t118 + (Ifges(7,3) * t80 + t120) * t81;
t23 = -t80 * t127 + t109;
t22 = -t109 + (pkin(9) * t80 + pkin(5)) * t84;
t21 = -t84 * qJ(6) + t24;
t1 = [m(2) + m(5) * (t11 ^ 2 + t128 + t8) + m(4) * (t16 ^ 2 + t31 ^ 2 + t128) + m(3) * (t79 ^ 2 + (t74 ^ 2 + t77 ^ 2) * t76 ^ 2) + t133 * (t4 ^ 2 + t6 ^ 2 + t8); m(3) * t79 + m(5) * (t35 * t11 - t14 * t113 + t7) + m(4) * (t78 * t31 + (-t14 * t85 + t16 * t82) * t75) + t133 * (t18 * t4 + t20 * t6 + t7); m(3) + m(5) * (t35 ^ 2 + t32 + t59) + m(4) * (t68 * t82 ^ 2 + t78 ^ 2 + t59) + t133 * (t18 ^ 2 + t20 ^ 2 + t32); t11 * t108 - t16 * mrSges(4,2) + t105 * t6 + t104 * t4 - t107 * t14 + t93 * t9 + m(7) * (t21 * t6 + t22 * t4 + t30 * t9) + m(6) * (pkin(9) * t124 - t23 * t4 + t24 * t6) + m(5) * (-pkin(3) * t14 + (t11 * t84 + t124) * pkin(9)); t35 * t108 + t105 * t20 + t104 * t18 + (-t82 * mrSges(4,2) + t107 * t85) * t75 + t93 * t33 + m(6) * (pkin(9) * t115 - t23 * t18 + t24 * t20) + m(7) * (t22 * t18 + t21 * t20 + t30 * t33) + m(5) * (pkin(3) * t113 + (t35 * t84 + t115) * pkin(9)); -0.2e1 * pkin(3) * t48 + 0.2e1 * t21 * t43 + 0.2e1 * t22 * t42 + 0.2e1 * t23 * t41 + 0.2e1 * t24 * t40 + 0.2e1 * t30 * t36 + Ifges(4,3) + (t71 + t73) * mrSges(5,3) * t134 + m(7) * (t21 ^ 2 + t22 ^ 2 + t30 ^ 2) + m(6) * (t23 ^ 2 + t24 ^ 2 + t66) + m(5) * (pkin(3) ^ 2 + t73 * t87 + t66) + ((Ifges(5,2) + t106) * t84 + t100) * t84 + (Ifges(5,1) * t81 + 0.2e1 * Ifges(5,4) * t84 + t37 * t134 + (t28 + t29) * t83 + (t26 - t27 + t119) * t80) * t81; -t11 * mrSges(5,2) + t101 * t9 + m(7) * (t44 * t9 + t97) + m(6) * (-pkin(4) * t9 + t97) + t132 * t135; -t35 * mrSges(5,2) + t101 * t33 + m(6) * (-pkin(4) * t33 + t94) + m(7) * (t44 * t33 + t94) + t132 * t136; -pkin(4) * t37 + t44 * t36 + (m(7) * t44 + t46) * t30 + (-pkin(9) * mrSges(5,2) - t64 / 0.2e1 - t63 / 0.2e1 - t62 / 0.2e1 + Ifges(5,6)) * t84 + (t24 * mrSges(6,3) + t21 * mrSges(7,2) + t118 / 0.2e1 - t26 / 0.2e1 + t27 / 0.2e1) * t83 + (t22 * mrSges(7,2) - t23 * mrSges(6,3) + t28 / 0.2e1 + t29 / 0.2e1) * t80 + (t105 * t83 + t104 * t80 + m(6) * (-t23 * t80 + t24 * t83) + m(7) * (t21 * t83 + t22 * t80)) * pkin(10) + (Ifges(5,5) + (t51 / 0.2e1 + t52 / 0.2e1) * t83 + (t49 / 0.2e1 - t50 / 0.2e1) * t80 + (-m(6) * pkin(4) + t103) * pkin(9)) * t81; -0.2e1 * pkin(4) * t47 + 0.2e1 * t44 * t46 + Ifges(5,3) + (-t49 + t50) * t83 + (t51 + t52) * t80 + m(7) * (t44 ^ 2 + t102) + m(6) * (pkin(4) ^ 2 + t102) + 0.2e1 * t132 * pkin(10) * t137; t129 * t6 + t130 * t4; t129 * t20 + t130 * t18; -Ifges(6,6) * t111 + t21 * mrSges(7,3) + qJ(6) * t43 + m(7) * (-pkin(5) * t22 + qJ(6) * t21) - t24 * mrSges(6,2) + t23 * mrSges(6,1) - t22 * mrSges(7,1) - pkin(5) * t42 - t106 * t84 - t100; -Ifges(7,6) * t83 + t62 + t63 + t64 + t89 * mrSges(7,2) + (m(7) * t89 - t90 - t91) * pkin(10); 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t106; m(7) * t4; m(7) * t18; m(7) * t22 + t42; (m(7) * pkin(10) + mrSges(7,2)) * t80; t131; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
