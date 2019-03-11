% Calculate joint inertia matrix for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:34
% EndTime: 2019-03-09 03:14:36
% DurationCPUTime: 0.98s
% Computational Cost: add. (1649->231), mult. (3204->316), div. (0->0), fcn. (3523->8), ass. (0->89)
t135 = Ifges(7,4) + Ifges(6,5);
t134 = Ifges(6,6) - Ifges(7,6);
t133 = Ifges(7,2) + Ifges(6,3);
t132 = m(7) + m(6);
t123 = cos(qJ(5));
t92 = sin(pkin(10));
t94 = cos(pkin(10));
t96 = sin(qJ(5));
t101 = t123 * t94 - t96 * t92;
t74 = t123 * t92 + t96 * t94;
t131 = t134 * t101 + t135 * t74;
t130 = -m(7) * pkin(5) - mrSges(7,1);
t124 = cos(qJ(3));
t117 = pkin(7) + qJ(2);
t93 = sin(pkin(9));
t77 = t117 * t93;
t95 = cos(pkin(9));
t79 = t117 * t95;
t97 = sin(qJ(3));
t52 = t124 * t77 + t79 * t97;
t129 = t52 ^ 2;
t91 = t95 ^ 2;
t128 = 0.2e1 * t52;
t85 = -pkin(2) * t95 - pkin(1);
t127 = 0.2e1 * t85;
t125 = t94 / 0.2e1;
t75 = t124 * t93 + t97 * t95;
t120 = t75 * t92;
t73 = -t124 * t95 + t93 * t97;
t39 = pkin(3) * t73 - qJ(4) * t75 + t85;
t55 = t124 * t79 - t97 * t77;
t17 = t92 * t39 + t94 * t55;
t13 = -pkin(8) * t120 + t17;
t119 = t75 * t94;
t16 = t94 * t39 - t55 * t92;
t7 = pkin(4) * t73 - pkin(8) * t119 + t16;
t4 = t123 * t13 + t96 * t7;
t122 = Ifges(5,4) * t92;
t121 = Ifges(5,4) * t94;
t116 = pkin(8) + qJ(4);
t33 = t74 * t75;
t18 = -mrSges(6,2) * t73 - mrSges(6,3) * t33;
t21 = -mrSges(7,2) * t33 + mrSges(7,3) * t73;
t115 = t18 + t21;
t34 = t101 * t75;
t19 = mrSges(6,1) * t73 - mrSges(6,3) * t34;
t20 = -t73 * mrSges(7,1) + t34 * mrSges(7,2);
t114 = -t19 + t20;
t37 = mrSges(5,1) * t120 + mrSges(5,2) * t119;
t111 = t92 ^ 2 + t94 ^ 2;
t110 = t93 ^ 2 + t91;
t107 = t116 * t92;
t78 = t116 * t94;
t50 = t107 * t123 + t78 * t96;
t54 = -t107 * t96 + t123 * t78;
t109 = t50 ^ 2 + t54 ^ 2;
t84 = -pkin(4) * t94 - pkin(3);
t106 = -t95 * mrSges(3,1) + t93 * mrSges(3,2);
t76 = -t94 * mrSges(5,1) + t92 * mrSges(5,2);
t15 = t33 * mrSges(6,1) + t34 * mrSges(6,2);
t44 = -mrSges(6,1) * t101 + t74 * mrSges(6,2);
t14 = t33 * mrSges(7,1) - t34 * mrSges(7,3);
t43 = -mrSges(7,1) * t101 - t74 * mrSges(7,3);
t24 = pkin(4) * t120 + t52;
t104 = pkin(5) * t101 + qJ(6) * t74;
t103 = -t16 * t92 + t17 * t94;
t102 = t133 * t73 - t134 * t33 + t135 * t34;
t3 = t123 * t7 - t96 * t13;
t100 = -t43 - t44;
t81 = Ifges(5,1) * t92 + t121;
t80 = Ifges(5,2) * t94 + t122;
t61 = t75 * mrSges(4,2);
t48 = Ifges(6,1) * t74 + Ifges(6,4) * t101;
t47 = Ifges(7,1) * t74 - Ifges(7,5) * t101;
t46 = Ifges(6,4) * t74 + Ifges(6,2) * t101;
t45 = Ifges(7,5) * t74 - Ifges(7,3) * t101;
t41 = mrSges(5,1) * t73 - mrSges(5,3) * t119;
t40 = -mrSges(5,2) * t73 - mrSges(5,3) * t120;
t38 = -t104 + t84;
t23 = t73 * Ifges(5,5) + (Ifges(5,1) * t94 - t122) * t75;
t22 = t73 * Ifges(5,6) + (-Ifges(5,2) * t92 + t121) * t75;
t11 = Ifges(6,1) * t34 - Ifges(6,4) * t33 + Ifges(6,5) * t73;
t10 = Ifges(7,1) * t34 + Ifges(7,4) * t73 + Ifges(7,5) * t33;
t9 = Ifges(6,4) * t34 - Ifges(6,2) * t33 + Ifges(6,6) * t73;
t8 = Ifges(7,5) * t34 + Ifges(7,6) * t73 + Ifges(7,3) * t33;
t5 = pkin(5) * t33 - qJ(6) * t34 + t24;
t2 = -t73 * pkin(5) - t3;
t1 = qJ(6) * t73 + t4;
t6 = [-0.2e1 * pkin(1) * t106 + Ifges(3,2) * t91 + t61 * t127 + t37 * t128 + 0.2e1 * t17 * t40 + 0.2e1 * t16 * t41 + 0.2e1 * t24 * t15 + 0.2e1 * t5 * t14 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t1 * t21 + Ifges(2,3) + (Ifges(3,1) * t93 + 0.2e1 * Ifges(3,4) * t95) * t93 + (t10 + t11) * t34 + (-t9 + t8) * t33 + 0.2e1 * t110 * qJ(2) * mrSges(3,3) + (mrSges(4,3) * t128 + Ifges(4,1) * t75 - t92 * t22 + t94 * t23) * t75 + (mrSges(4,1) * t127 - 0.2e1 * t55 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t73 + (Ifges(5,5) * t94 - Ifges(5,6) * t92 - (2 * Ifges(4,4))) * t75 + t102) * t73 + m(3) * (qJ(2) ^ 2 * t110 + pkin(1) ^ 2) + m(4) * (t55 ^ 2 + t85 ^ 2 + t129) + m(5) * (t16 ^ 2 + t17 ^ 2 + t129) + m(6) * (t24 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); -m(3) * pkin(1) + t73 * mrSges(4,1) + t92 * t40 + t94 * t41 + t61 + t115 * t74 - t114 * t101 + m(6) * (t101 * t3 + t4 * t74) + m(7) * (t1 * t74 - t101 * t2) + m(5) * (t16 * t94 + t17 * t92) + m(4) * t85 + t106; m(5) * t111 + m(3) + m(4) + t132 * (t101 ^ 2 + t74 ^ 2); (t48 / 0.2e1 + t47 / 0.2e1) * t34 + (t45 / 0.2e1 - t46 / 0.2e1) * t33 + t103 * mrSges(5,3) + m(5) * (-pkin(3) * t52 + qJ(4) * t103) + (t76 - mrSges(4,1)) * t52 + (t94 * t40 - t92 * t41) * qJ(4) + t92 * t23 / 0.2e1 + t84 * t15 - Ifges(4,6) * t73 - t55 * mrSges(4,2) - pkin(3) * t37 + t38 * t14 + t5 * t43 + t24 * t44 + t22 * t125 + (Ifges(4,5) - t92 * t80 / 0.2e1 + t81 * t125) * t75 + t114 * t50 + t115 * t54 + m(6) * (t24 * t84 - t3 * t50 + t4 * t54) + m(7) * (t1 * t54 + t2 * t50 + t38 * t5) - (-t9 / 0.2e1 + t8 / 0.2e1 - t4 * mrSges(6,3) - t1 * mrSges(7,2)) * t101 + (t10 / 0.2e1 + t11 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t74 + (Ifges(5,5) * t92 + Ifges(5,6) * t94 + t131) * t73 / 0.2e1; t132 * (-t101 * t50 + t54 * t74); -0.2e1 * pkin(3) * t76 + 0.2e1 * t38 * t43 + 0.2e1 * t84 * t44 + t94 * t80 + t92 * t81 + Ifges(4,3) + m(7) * (t38 ^ 2 + t109) + m(6) * (t84 ^ 2 + t109) + m(5) * (qJ(4) ^ 2 * t111 + pkin(3) ^ 2) + (t47 + t48) * t74 - (t45 - t46) * t101 + 0.2e1 * mrSges(5,3) * qJ(4) * t111 + 0.2e1 * (t101 * t54 + t50 * t74) * (mrSges(7,2) + mrSges(6,3)); m(5) * t52 + m(6) * t24 + m(7) * t5 + t14 + t15 + t37; 0; -m(5) * pkin(3) + m(6) * t84 + m(7) * t38 - t100 + t76; m(5) + t132; -pkin(5) * t20 + qJ(6) * t21 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) - t4 * mrSges(6,2) - t2 * mrSges(7,1) + t3 * mrSges(6,1) + t102; m(7) * t104 + t100; (-pkin(5) * t74 + qJ(6) * t101) * mrSges(7,2) + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t54 + (-mrSges(6,1) + t130) * t50 + t131; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t133; m(7) * t2 + t20; -m(7) * t101; m(7) * t50 + t74 * mrSges(7,2); 0; t130; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
