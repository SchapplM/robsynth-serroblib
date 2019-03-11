% Calculate joint inertia matrix for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:41
% EndTime: 2019-03-09 06:18:44
% DurationCPUTime: 1.19s
% Computational Cost: add. (1978->268), mult. (3781->373), div. (0->0), fcn. (4122->8), ass. (0->95)
t144 = Ifges(7,4) + Ifges(6,5);
t143 = -Ifges(6,6) + Ifges(7,6);
t131 = Ifges(7,2) + Ifges(6,3);
t142 = m(6) + m(7);
t103 = cos(qJ(4));
t101 = sin(qJ(3));
t133 = cos(qJ(3));
t97 = sin(pkin(10));
t98 = cos(pkin(10));
t72 = t101 * t98 + t133 * t97;
t119 = t103 * t72;
t71 = t101 * t97 - t133 * t98;
t141 = Ifges(5,5) * t119 + Ifges(5,3) * t71;
t100 = sin(qJ(4));
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t74 = t100 * t99 - t102 * t103;
t76 = t100 * t102 + t103 * t99;
t140 = t143 * t74 + t144 * t76;
t130 = pkin(7) + qJ(2);
t78 = t130 * t97;
t79 = t130 * t98;
t49 = t101 * t79 + t133 * t78;
t139 = t49 ^ 2;
t94 = t98 ^ 2;
t138 = 2 * mrSges(7,3);
t137 = 0.2e1 * t49;
t86 = -pkin(2) * t98 - pkin(1);
t136 = 0.2e1 * t86;
t134 = -pkin(9) - pkin(8);
t120 = t100 * t72;
t38 = pkin(3) * t71 - pkin(8) * t72 + t86;
t51 = -t101 * t78 + t133 * t79;
t19 = t100 * t38 + t103 * t51;
t15 = -pkin(9) * t120 + t19;
t18 = -t100 * t51 + t103 * t38;
t9 = pkin(4) * t71 - pkin(9) * t119 + t18;
t6 = t102 * t15 + t99 * t9;
t32 = t76 * t72;
t20 = -mrSges(6,2) * t71 - mrSges(6,3) * t32;
t23 = -mrSges(7,2) * t32 + mrSges(7,3) * t71;
t129 = t20 + t23;
t33 = t74 * t72;
t21 = mrSges(6,1) * t71 + mrSges(6,3) * t33;
t22 = -t71 * mrSges(7,1) - t33 * mrSges(7,2);
t128 = -t21 + t22;
t125 = Ifges(5,5) * t100 + Ifges(5,6) * t103;
t124 = t97 ^ 2 + t94;
t123 = t100 ^ 2 + t103 ^ 2;
t122 = Ifges(5,4) * t100;
t121 = Ifges(5,4) * t103;
t116 = t134 * t100;
t83 = t134 * t103;
t55 = -t102 * t116 - t83 * t99;
t57 = -t102 * t83 + t99 * t116;
t118 = t55 ^ 2 + t57 ^ 2;
t89 = -pkin(4) * t103 - pkin(3);
t115 = -t98 * mrSges(3,1) + t97 * mrSges(3,2);
t43 = t74 * mrSges(7,1) - t76 * mrSges(7,3);
t44 = t74 * mrSges(6,1) + t76 * mrSges(6,2);
t26 = pkin(4) * t120 + t49;
t113 = -pkin(5) * t74 + qJ(6) * t76;
t5 = t102 * t9 - t15 * t99;
t80 = -t103 * mrSges(5,1) + t100 * mrSges(5,2);
t112 = mrSges(5,1) * t100 + mrSges(5,2) * t103;
t111 = t131 * t71 + t143 * t32 - t144 * t33;
t110 = (mrSges(6,1) * t102 - mrSges(6,2) * t99) * pkin(4);
t109 = -t43 - t44;
t108 = t140 + (-mrSges(6,2) + mrSges(7,3)) * t57 + (-mrSges(6,1) - mrSges(7,1)) * t55;
t2 = qJ(6) * t71 + t6;
t3 = -pkin(5) * t71 - t5;
t107 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) + t111;
t88 = -pkin(4) * t102 - pkin(5);
t85 = pkin(4) * t99 + qJ(6);
t82 = Ifges(5,1) * t100 + t121;
t81 = Ifges(5,2) * t103 + t122;
t60 = t72 * mrSges(4,2);
t48 = Ifges(6,1) * t76 - Ifges(6,4) * t74;
t47 = Ifges(7,1) * t76 + Ifges(7,5) * t74;
t46 = Ifges(6,4) * t76 - Ifges(6,2) * t74;
t45 = Ifges(7,5) * t76 + Ifges(7,3) * t74;
t41 = mrSges(5,1) * t71 - mrSges(5,3) * t119;
t40 = -mrSges(5,2) * t71 - mrSges(5,3) * t120;
t39 = -t113 + t89;
t37 = t112 * t72;
t25 = Ifges(5,5) * t71 + (Ifges(5,1) * t103 - t122) * t72;
t24 = Ifges(5,6) * t71 + (-Ifges(5,2) * t100 + t121) * t72;
t17 = mrSges(6,1) * t32 - mrSges(6,2) * t33;
t16 = mrSges(7,1) * t32 + mrSges(7,3) * t33;
t13 = -Ifges(6,1) * t33 - Ifges(6,4) * t32 + Ifges(6,5) * t71;
t12 = -Ifges(7,1) * t33 + Ifges(7,4) * t71 + Ifges(7,5) * t32;
t11 = -Ifges(6,4) * t33 - Ifges(6,2) * t32 + Ifges(6,6) * t71;
t10 = -Ifges(7,5) * t33 + Ifges(7,6) * t71 + Ifges(7,3) * t32;
t7 = pkin(5) * t32 + qJ(6) * t33 + t26;
t1 = [-0.2e1 * pkin(1) * t115 + Ifges(3,2) * t94 + t60 * t136 + 0.2e1 * t19 * t40 + 0.2e1 * t18 * t41 + t37 * t137 + 0.2e1 * t6 * t20 + 0.2e1 * t5 * t21 + 0.2e1 * t3 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t26 * t17 + 0.2e1 * t7 * t16 + Ifges(2,3) + (Ifges(3,1) * t97 + 0.2e1 * Ifges(3,4) * t98) * t97 - (t12 + t13) * t33 + (t10 - t11) * t32 + 0.2e1 * t124 * qJ(2) * mrSges(3,3) + (mrSges(4,3) * t137 + Ifges(4,1) * t72 - t100 * t24 + t103 * t25) * t72 + (mrSges(4,1) * t136 - 0.2e1 * t51 * mrSges(4,3) + Ifges(4,2) * t71 + (-Ifges(5,6) * t100 - (2 * Ifges(4,4))) * t72 + t111 + t141) * t71 + m(3) * (t124 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t51 ^ 2 + t86 ^ 2 + t139) + m(5) * (t18 ^ 2 + t19 ^ 2 + t139) + m(6) * (t26 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2); -m(3) * pkin(1) + t71 * mrSges(4,1) + t100 * t40 + t103 * t41 + t60 + t129 * t76 + t128 * t74 + m(7) * (t2 * t76 + t3 * t74) + m(6) * (-t5 * t74 + t6 * t76) + m(5) * (t100 * t19 + t103 * t18) + m(4) * t86 + t115; m(5) * t123 + m(3) + m(4) + t142 * (t74 ^ 2 + t76 ^ 2); t89 * t17 - Ifges(4,6) * t71 + Ifges(4,5) * t72 + t39 * t16 + t7 * t43 + t26 * t44 - t51 * mrSges(4,2) - pkin(3) * t37 + t129 * t57 + t128 * t55 + (t80 - mrSges(4,1)) * t49 - (t47 / 0.2e1 + t48 / 0.2e1) * t33 + (t45 / 0.2e1 - t46 / 0.2e1) * t32 + (t72 * t82 / 0.2e1 + pkin(8) * t40 + t19 * mrSges(5,3) + t24 / 0.2e1) * t103 + (-t72 * t81 / 0.2e1 - pkin(8) * t41 - t18 * mrSges(5,3) + t25 / 0.2e1) * t100 + m(5) * (-pkin(3) * t49 + (-t18 * t100 + t19 * t103) * pkin(8)) + m(6) * (t26 * t89 - t5 * t55 + t57 * t6) + m(7) * (t2 * t57 + t3 * t55 + t39 * t7) + (t3 * mrSges(7,2) - t5 * mrSges(6,3) + t12 / 0.2e1 + t13 / 0.2e1) * t76 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t10 / 0.2e1 - t11 / 0.2e1) * t74 + (t125 + t140) * t71 / 0.2e1; t142 * (t55 * t74 + t57 * t76); -0.2e1 * pkin(3) * t80 + t100 * t82 + t103 * t81 + 0.2e1 * t39 * t43 + 0.2e1 * t89 * t44 + Ifges(4,3) + m(6) * (t89 ^ 2 + t118) + m(7) * (t39 ^ 2 + t118) + m(5) * (t123 * pkin(8) ^ 2 + pkin(3) ^ 2) + (t47 + t48) * t76 + (t45 - t46) * t74 + 0.2e1 * t123 * pkin(8) * mrSges(5,3) + 0.2e1 * (t55 * t76 - t57 * t74) * (mrSges(7,2) + mrSges(6,3)); m(7) * (t2 * t85 + t3 * t88) + (m(6) * (t102 * t5 + t6 * t99) + t102 * t21 + t99 * t20) * pkin(4) - Ifges(5,6) * t120 + t88 * t22 + t85 * t23 + t18 * mrSges(5,1) - t19 * mrSges(5,2) + t107 + t141; m(7) * (t74 * t88 + t76 * t85) + m(6) * (-t102 * t74 + t76 * t99) * pkin(4) + t109 - t80; m(7) * (t55 * t88 + t57 * t85) - t112 * pkin(8) + (-t74 * t85 + t76 * t88) * mrSges(7,2) + (m(6) * (-t102 * t55 + t57 * t99) + (-t102 * t76 - t74 * t99) * mrSges(6,3)) * pkin(4) + t108 + t125; -0.2e1 * t88 * mrSges(7,1) + t85 * t138 + Ifges(5,3) + 0.2e1 * t110 + m(7) * (t85 ^ 2 + t88 ^ 2) + m(6) * (t102 ^ 2 + t99 ^ 2) * pkin(4) ^ 2 + t131; -pkin(5) * t22 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) + qJ(6) * t23 + t107; m(7) * t113 + t109; m(7) * (-pkin(5) * t55 + qJ(6) * t57) + (-pkin(5) * t76 - qJ(6) * t74) * mrSges(7,2) + t108; m(7) * (-pkin(5) * t88 + qJ(6) * t85) + t110 + (t85 + qJ(6)) * mrSges(7,3) + (-t88 + pkin(5)) * mrSges(7,1) + t131; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t138 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t131; m(7) * t3 + t22; m(7) * t74; m(7) * t55 + t76 * mrSges(7,2); m(7) * t88 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
