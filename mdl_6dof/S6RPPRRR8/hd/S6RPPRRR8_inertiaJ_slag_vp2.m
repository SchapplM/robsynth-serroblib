% Calculate joint inertia matrix for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:26
% EndTime: 2019-03-09 02:35:27
% DurationCPUTime: 0.75s
% Computational Cost: add. (1424->223), mult. (2569->323), div. (0->0), fcn. (2765->8), ass. (0->85)
t74 = sin(pkin(10));
t79 = sin(qJ(4));
t109 = cos(qJ(4));
t75 = cos(pkin(10));
t94 = t109 * t75;
t50 = t74 * t79 - t94;
t81 = cos(qJ(5));
t104 = t50 * t81;
t103 = t79 * t75;
t51 = t109 * t74 + t103;
t118 = -Ifges(6,5) * t104 + Ifges(6,3) * t51;
t76 = -pkin(1) - qJ(3);
t110 = -pkin(7) + t76;
t55 = t110 * t74;
t30 = -t110 * t94 + t55 * t79;
t117 = t30 ^ 2;
t116 = t50 ^ 2;
t70 = t75 ^ 2;
t115 = -2 * mrSges(5,3);
t113 = m(7) * pkin(5);
t111 = -pkin(9) - pkin(8);
t78 = sin(qJ(5));
t108 = Ifges(6,4) * t78;
t107 = Ifges(6,4) * t81;
t106 = t30 * t50;
t105 = t50 * t78;
t61 = t74 * pkin(3) + qJ(2);
t25 = pkin(4) * t51 + pkin(8) * t50 + t61;
t32 = t110 * t103 + t109 * t55;
t10 = t78 * t25 + t81 * t32;
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t54 = t77 * t81 + t78 * t80;
t88 = t77 * t78 - t80 * t81;
t102 = Ifges(7,5) * t54 - Ifges(7,6) * t88;
t56 = -t81 * mrSges(6,1) + t78 * mrSges(6,2);
t101 = t56 - mrSges(5,1);
t100 = t74 * mrSges(4,1) + t75 * mrSges(4,2);
t99 = Ifges(6,5) * t78 + Ifges(6,6) * t81;
t98 = t74 ^ 2 + t70;
t97 = t78 ^ 2 + t81 ^ 2;
t22 = t54 * t50;
t24 = t88 * t50;
t96 = Ifges(7,5) * t24 + Ifges(7,6) * t22 + Ifges(7,3) * t51;
t95 = m(4) * t98;
t93 = t98 * mrSges(4,3);
t92 = t97 * mrSges(6,3);
t21 = t54 * t51;
t23 = t88 * t51;
t91 = -t21 * mrSges(7,1) + t23 * mrSges(7,2);
t33 = mrSges(7,1) * t88 + t54 * mrSges(7,2);
t9 = t81 * t25 - t32 * t78;
t90 = t10 * t81 - t78 * t9;
t89 = -mrSges(6,1) * t78 - mrSges(6,2) * t81;
t59 = t111 * t78;
t60 = t111 * t81;
t37 = t59 * t80 + t60 * t77;
t38 = t59 * t77 - t60 * t80;
t87 = t37 * mrSges(7,1) - t38 * mrSges(7,2) + t102;
t4 = pkin(5) * t51 + pkin(9) * t104 + t9;
t7 = pkin(9) * t105 + t10;
t2 = t4 * t80 - t7 * t77;
t3 = t4 * t77 + t7 * t80;
t86 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t96;
t85 = (mrSges(7,1) * t80 - mrSges(7,2) * t77) * pkin(5);
t82 = qJ(2) ^ 2;
t63 = -pkin(5) * t81 - pkin(4);
t58 = Ifges(6,1) * t78 + t107;
t57 = Ifges(6,2) * t81 + t108;
t48 = t51 ^ 2;
t41 = t50 * mrSges(5,2);
t35 = Ifges(7,1) * t54 - Ifges(7,4) * t88;
t34 = Ifges(7,4) * t54 - Ifges(7,2) * t88;
t29 = mrSges(6,1) * t51 + mrSges(6,3) * t104;
t28 = -mrSges(6,2) * t51 + mrSges(6,3) * t105;
t26 = t89 * t50;
t15 = -pkin(5) * t105 + t30;
t14 = Ifges(6,5) * t51 + (-Ifges(6,1) * t81 + t108) * t50;
t13 = Ifges(6,6) * t51 + (Ifges(6,2) * t78 - t107) * t50;
t12 = mrSges(7,1) * t51 - mrSges(7,3) * t24;
t11 = -mrSges(7,2) * t51 + mrSges(7,3) * t22;
t8 = -mrSges(7,1) * t22 + mrSges(7,2) * t24;
t6 = Ifges(7,1) * t24 + Ifges(7,4) * t22 + Ifges(7,5) * t51;
t5 = Ifges(7,4) * t24 + Ifges(7,2) * t22 + Ifges(7,6) * t51;
t1 = [Ifges(4,1) * t70 - 0.2e1 * t61 * t41 + 0.2e1 * t30 * t26 + t22 * t5 + t24 * t6 + 0.2e1 * t10 * t28 + 0.2e1 * t9 * t29 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t15 * t8 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) + (-0.2e1 * Ifges(4,4) * t75 + Ifges(4,2) * t74) * t74 - 0.2e1 * t76 * t93 + (0.2e1 * mrSges(5,1) * t61 + Ifges(5,2) * t51 + t32 * t115 + t118 + t96) * t51 + (t30 * t115 + Ifges(5,1) * t50 + t78 * t13 - t81 * t14 + (Ifges(6,6) * t78 + (2 * Ifges(5,4))) * t51) * t50 + m(4) * (t98 * t76 ^ 2 + t82) + m(3) * ((pkin(1) ^ 2) + t82) + m(5) * (t32 ^ 2 + t61 ^ 2 + t117) + m(6) * (t10 ^ 2 + t9 ^ 2 + t117) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + 0.2e1 * (t100 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t116 * mrSges(5,3) - t23 * t11 - t21 * t12 + mrSges(3,2) + (t26 + t8) * t50 - t93 + (-mrSges(5,3) * t51 + t81 * t28 - t78 * t29) * t51 + m(7) * (t15 * t50 - t2 * t21 - t23 * t3) + m(6) * (t90 * t51 + t106) + m(5) * (t32 * t51 + t106) + t76 * t95; m(3) + m(7) * (t21 ^ 2 + t23 ^ 2 + t116) + m(6) * (t97 * t48 + t116) + m(5) * (t48 + t116) + t95; m(4) * qJ(2) + t51 * mrSges(5,1) + t54 * t11 - t88 * t12 + t78 * t28 + t81 * t29 - t41 + m(7) * (-t2 * t88 + t3 * t54) + m(6) * (t10 * t78 + t81 * t9) + m(5) * t61 + t100; m(7) * (t21 * t88 - t23 * t54); m(4) + m(5) + m(6) * t97 + m(7) * (t54 ^ 2 + t88 ^ 2); t63 * t8 - Ifges(5,6) * t51 - t88 * t5 / 0.2e1 + t54 * t6 / 0.2e1 - t32 * mrSges(5,2) + t15 * t33 + t22 * t34 / 0.2e1 + t24 * t35 / 0.2e1 + t37 * t12 + t38 * t11 - pkin(4) * t26 + t101 * t30 + (-t2 * t54 - t3 * t88) * mrSges(7,3) + (t13 / 0.2e1 + pkin(8) * t28 + t10 * mrSges(6,3)) * t81 + (t14 / 0.2e1 - pkin(8) * t29 - t9 * mrSges(6,3)) * t78 + m(6) * (-pkin(4) * t30 + t90 * pkin(8)) + m(7) * (t15 * t63 + t2 * t37 + t38 * t3) + (-Ifges(5,5) - t81 * t58 / 0.2e1 + t78 * t57 / 0.2e1) * t50 + (t99 + t102) * t51 / 0.2e1; (t21 * t54 + t23 * t88) * mrSges(7,3) + (-mrSges(5,2) + t92) * t51 + (t33 + t101) * t50 + m(7) * (-t21 * t37 - t23 * t38 + t50 * t63) + m(6) * (t97 * t51 * pkin(8) - pkin(4) * t50); m(7) * (-t37 * t88 + t38 * t54); -0.2e1 * pkin(4) * t56 + 0.2e1 * t63 * t33 - t88 * t34 + t54 * t35 + t81 * t57 + t78 * t58 + Ifges(5,3) + m(7) * (t37 ^ 2 + t38 ^ 2 + t63 ^ 2) + m(6) * (t97 * pkin(8) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t37 * t54 - t38 * t88) * mrSges(7,3) + 0.2e1 * pkin(8) * t92; Ifges(6,6) * t105 + t9 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (t2 * t80 + t3 * t77) + t77 * t11 + t80 * t12) * pkin(5) + t86 + t118; t89 * t51 + (-t21 * t80 - t23 * t77) * t113 + t91; (t54 * t77 - t80 * t88) * t113 - t56 - t33; t89 * pkin(8) + (m(7) * (t37 * t80 + t38 * t77) + (-t54 * t80 - t77 * t88) * mrSges(7,3)) * pkin(5) + t87 + t99; Ifges(6,3) + Ifges(7,3) + m(7) * (t77 ^ 2 + t80 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t85; t86; t91; -t33; t87; Ifges(7,3) + t85; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
