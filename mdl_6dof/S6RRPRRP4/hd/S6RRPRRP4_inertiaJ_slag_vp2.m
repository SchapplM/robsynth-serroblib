% Calculate joint inertia matrix for
% S6RRPRRP4
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:37
% EndTime: 2019-03-09 11:52:40
% DurationCPUTime: 1.29s
% Computational Cost: add. (2087->281), mult. (3922->396), div. (0->0), fcn. (4244->8), ass. (0->95)
t145 = Ifges(7,4) + Ifges(6,5);
t144 = -Ifges(6,6) + Ifges(7,6);
t133 = Ifges(7,2) + Ifges(6,3);
t143 = m(7) + m(6);
t104 = cos(qJ(4));
t102 = sin(qJ(2));
t105 = cos(qJ(2));
t98 = sin(pkin(10));
t99 = cos(pkin(10));
t73 = t102 * t99 + t105 * t98;
t121 = t104 * t73;
t72 = t102 * t98 - t99 * t105;
t142 = Ifges(5,5) * t121 + Ifges(5,3) * t72;
t100 = sin(qJ(5));
t101 = sin(qJ(4));
t103 = cos(qJ(5));
t76 = t100 * t101 - t103 * t104;
t78 = t100 * t104 + t101 * t103;
t141 = t144 * t76 + t145 * t78;
t132 = -qJ(3) - pkin(7);
t116 = t132 * t102;
t81 = t132 * t105;
t55 = -t99 * t116 - t81 * t98;
t140 = t55 ^ 2;
t139 = 2 * mrSges(7,3);
t138 = 0.2e1 * t55;
t91 = -pkin(2) * t105 - pkin(1);
t137 = 0.2e1 * t91;
t86 = pkin(2) * t98 + pkin(8);
t135 = pkin(9) + t86;
t122 = t101 * t73;
t45 = pkin(3) * t72 - pkin(8) * t73 + t91;
t57 = t98 * t116 - t99 * t81;
t19 = t101 * t45 + t104 * t57;
t15 = -pkin(9) * t122 + t19;
t18 = -t101 * t57 + t104 * t45;
t9 = pkin(4) * t72 - pkin(9) * t121 + t18;
t6 = t100 * t9 + t103 * t15;
t33 = t78 * t73;
t20 = -mrSges(6,2) * t72 - mrSges(6,3) * t33;
t23 = -mrSges(7,2) * t33 + mrSges(7,3) * t72;
t131 = t20 + t23;
t34 = t76 * t73;
t21 = mrSges(6,1) * t72 + mrSges(6,3) * t34;
t22 = -t72 * mrSges(7,1) - t34 * mrSges(7,2);
t130 = -t21 + t22;
t127 = Ifges(5,5) * t101 + Ifges(5,6) * t104;
t126 = t101 ^ 2 + t104 ^ 2;
t125 = t102 ^ 2 + t105 ^ 2;
t124 = Ifges(5,4) * t101;
t123 = Ifges(5,4) * t104;
t118 = t135 * t101;
t65 = t135 * t104;
t42 = t100 * t65 + t103 * t118;
t44 = -t100 * t118 + t103 * t65;
t120 = t42 ^ 2 + t44 ^ 2;
t88 = -pkin(2) * t99 - pkin(3);
t49 = t76 * mrSges(7,1) - t78 * mrSges(7,3);
t50 = t76 * mrSges(6,1) + t78 * mrSges(6,2);
t27 = pkin(4) * t122 + t55;
t115 = -pkin(5) * t76 + qJ(6) * t78;
t5 = -t100 * t15 + t103 * t9;
t80 = -t104 * mrSges(5,1) + t101 * mrSges(5,2);
t114 = mrSges(5,1) * t101 + mrSges(5,2) * t104;
t113 = t133 * t72 + t144 * t33 - t145 * t34;
t79 = -pkin(4) * t104 + t88;
t112 = (mrSges(6,1) * t103 - mrSges(6,2) * t100) * pkin(4);
t111 = -t49 - t50;
t110 = t141 + (-mrSges(6,2) + mrSges(7,3)) * t44 + (-mrSges(6,1) - mrSges(7,1)) * t42;
t2 = qJ(6) * t72 + t6;
t3 = -pkin(5) * t72 - t5;
t109 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) + t113;
t90 = -pkin(4) * t103 - pkin(5);
t87 = pkin(4) * t100 + qJ(6);
t83 = Ifges(5,1) * t101 + t123;
t82 = Ifges(5,2) * t104 + t124;
t61 = t73 * mrSges(4,2);
t54 = Ifges(6,1) * t78 - Ifges(6,4) * t76;
t53 = Ifges(7,1) * t78 + Ifges(7,5) * t76;
t52 = Ifges(6,4) * t78 - Ifges(6,2) * t76;
t51 = Ifges(7,5) * t78 + Ifges(7,3) * t76;
t47 = mrSges(5,1) * t72 - mrSges(5,3) * t121;
t46 = -mrSges(5,2) * t72 - mrSges(5,3) * t122;
t41 = t114 * t73;
t35 = -t115 + t79;
t25 = Ifges(5,5) * t72 + (Ifges(5,1) * t104 - t124) * t73;
t24 = Ifges(5,6) * t72 + (-Ifges(5,2) * t101 + t123) * t73;
t17 = mrSges(6,1) * t33 - mrSges(6,2) * t34;
t16 = mrSges(7,1) * t33 + mrSges(7,3) * t34;
t13 = -Ifges(6,1) * t34 - Ifges(6,4) * t33 + Ifges(6,5) * t72;
t12 = -Ifges(7,1) * t34 + Ifges(7,4) * t72 + Ifges(7,5) * t33;
t11 = -Ifges(6,4) * t34 - Ifges(6,2) * t33 + Ifges(6,6) * t72;
t10 = -Ifges(7,5) * t34 + Ifges(7,6) * t72 + Ifges(7,3) * t33;
t7 = pkin(5) * t33 + qJ(6) * t34 + t27;
t1 = [-0.2e1 * pkin(1) * (-t105 * mrSges(3,1) + t102 * mrSges(3,2)) + t102 * (Ifges(3,1) * t102 + Ifges(3,4) * t105) + t105 * (Ifges(3,4) * t102 + Ifges(3,2) * t105) + t61 * t137 + t41 * t138 + 0.2e1 * t19 * t46 + 0.2e1 * t18 * t47 + 0.2e1 * t27 * t17 + 0.2e1 * t6 * t20 + 0.2e1 * t5 * t21 + 0.2e1 * t3 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t7 * t16 + Ifges(2,3) - (t12 + t13) * t34 + (t10 - t11) * t33 + 0.2e1 * t125 * pkin(7) * mrSges(3,3) + (mrSges(4,3) * t138 + Ifges(4,1) * t73 - t101 * t24 + t104 * t25) * t73 + (mrSges(4,1) * t137 - 0.2e1 * t57 * mrSges(4,3) + Ifges(4,2) * t72 + (-Ifges(5,6) * t101 - (2 * Ifges(4,4))) * t73 + t113 + t142) * t72 + m(3) * (t125 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t57 ^ 2 + t91 ^ 2 + t140) + m(5) * (t18 ^ 2 + t19 ^ 2 + t140) + m(6) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2); (t80 - mrSges(4,1)) * t55 + m(5) * (t55 * t88 + (-t18 * t101 + t19 * t104) * t86) + (t24 / 0.2e1 + t73 * t83 / 0.2e1 + t86 * t46 + t19 * mrSges(5,3)) * t104 + (t25 / 0.2e1 - t73 * t82 / 0.2e1 - t86 * t47 - t18 * mrSges(5,3)) * t101 + m(6) * (t27 * t79 - t42 * t5 + t44 * t6) + m(7) * (t2 * t44 + t3 * t42 + t35 * t7) + (t3 * mrSges(7,2) - t5 * mrSges(6,3) + t12 / 0.2e1 + t13 / 0.2e1) * t78 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) - t11 / 0.2e1 + t10 / 0.2e1) * t76 + (m(4) * (-t55 * t99 + t57 * t98) + (-t72 * t98 - t73 * t99) * mrSges(4,3)) * pkin(2) - (t53 / 0.2e1 + t54 / 0.2e1) * t34 + (t51 / 0.2e1 - t52 / 0.2e1) * t33 + Ifges(3,6) * t105 + Ifges(3,5) * t102 + t88 * t41 + (-t102 * mrSges(3,1) - t105 * mrSges(3,2)) * pkin(7) + t79 * t17 - Ifges(4,6) * t72 + Ifges(4,5) * t73 + t7 * t49 + t27 * t50 - t57 * mrSges(4,2) + t35 * t16 + t130 * t42 + t131 * t44 + (t127 + t141) * t72 / 0.2e1; t101 * t83 + t104 * t82 + 0.2e1 * t35 * t49 + 0.2e1 * t79 * t50 + 0.2e1 * t88 * t80 + Ifges(3,3) + Ifges(4,3) + (t53 + t54) * t78 + (t51 - t52) * t76 + m(7) * (t35 ^ 2 + t120) + m(6) * (t79 ^ 2 + t120) + m(5) * (t126 * t86 ^ 2 + t88 ^ 2) + m(4) * (t98 ^ 2 + t99 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t99 - mrSges(4,2) * t98) * pkin(2) + 0.2e1 * t126 * t86 * mrSges(5,3) + 0.2e1 * (t42 * t78 - t44 * t76) * (mrSges(7,2) + mrSges(6,3)); t72 * mrSges(4,1) + t101 * t46 + t104 * t47 + t61 + t131 * t78 + t130 * t76 + m(7) * (t2 * t78 + t3 * t76) + m(6) * (-t5 * t76 + t6 * t78) + m(5) * (t101 * t19 + t104 * t18) + m(4) * t91; t143 * (t42 * t76 + t78 * t44); m(5) * t126 + m(4) + t143 * (t76 ^ 2 + t78 ^ 2); (m(6) * (t100 * t6 + t103 * t5) + t100 * t20 + t103 * t21) * pkin(4) + m(7) * (t2 * t87 + t3 * t90) + t87 * t23 + t90 * t22 - Ifges(5,6) * t122 + t109 + t18 * mrSges(5,1) - t19 * mrSges(5,2) + t142; m(7) * (t42 * t90 + t44 * t87) - t114 * t86 + (-t76 * t87 + t78 * t90) * mrSges(7,2) + (m(6) * (t100 * t44 - t103 * t42) + (-t100 * t76 - t103 * t78) * mrSges(6,3)) * pkin(4) + t110 + t127; m(7) * (t76 * t90 + t78 * t87) + m(6) * (t100 * t78 - t103 * t76) * pkin(4) + t111 - t80; -0.2e1 * t90 * mrSges(7,1) + t87 * t139 + Ifges(5,3) + 0.2e1 * t112 + m(7) * (t87 ^ 2 + t90 ^ 2) + m(6) * (t100 ^ 2 + t103 ^ 2) * pkin(4) ^ 2 + t133; m(7) * (-pkin(5) * t3 + qJ(6) * t2) - pkin(5) * t22 + qJ(6) * t23 + t109; m(7) * (-pkin(5) * t42 + qJ(6) * t44) + (-pkin(5) * t78 - qJ(6) * t76) * mrSges(7,2) + t110; m(7) * t115 + t111; m(7) * (-pkin(5) * t90 + qJ(6) * t87) + t112 + (qJ(6) + t87) * mrSges(7,3) + (pkin(5) - t90) * mrSges(7,1) + t133; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t139 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t133; m(7) * t3 + t22; m(7) * t42 + t78 * mrSges(7,2); m(7) * t76; m(7) * t90 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
