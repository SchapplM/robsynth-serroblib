% Calculate joint inertia matrix for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:05:43
% EndTime: 2018-11-23 16:05:43
% DurationCPUTime: 0.84s
% Computational Cost: add. (2548->254), mult. (4948->371), div. (0->0), fcn. (5723->10), ass. (0->97)
t101 = sin(pkin(11));
t103 = cos(pkin(11));
t106 = sin(qJ(5));
t109 = cos(qJ(5));
t83 = t101 * t109 + t103 * t106;
t102 = sin(pkin(10));
t104 = cos(pkin(10));
t107 = sin(qJ(3));
t131 = cos(qJ(3));
t84 = t102 * t131 + t107 * t104;
t45 = t83 * t84;
t81 = -t101 * t106 + t103 * t109;
t46 = t81 * t84;
t82 = t102 * t107 - t104 * t131;
t136 = Ifges(6,5) * t46 - Ifges(6,6) * t45 + Ifges(6,3) * t82;
t128 = pkin(7) + qJ(2);
t87 = t128 * t102;
t89 = t128 * t104;
t63 = t107 * t89 + t131 * t87;
t135 = t63 ^ 2;
t134 = 0.2e1 * t63;
t94 = -pkin(2) * t104 - pkin(1);
t133 = 0.2e1 * t94;
t100 = t104 ^ 2;
t130 = t82 * Ifges(5,5);
t129 = t82 * Ifges(5,6);
t127 = pkin(8) + qJ(4);
t120 = t103 * t84;
t50 = pkin(3) * t82 - qJ(4) * t84 + t94;
t66 = -t107 * t87 + t131 * t89;
t33 = -t101 * t66 + t103 * t50;
t13 = pkin(4) * t82 - pkin(8) * t120 + t33;
t121 = t101 * t84;
t34 = t101 * t50 + t103 * t66;
t22 = -pkin(8) * t121 + t34;
t7 = t106 * t13 + t109 * t22;
t105 = sin(qJ(6));
t108 = cos(qJ(6));
t54 = -t105 * t83 + t108 * t81;
t55 = t105 * t81 + t108 * t83;
t126 = Ifges(7,5) * t55 + Ifges(7,6) * t54;
t49 = mrSges(5,1) * t121 + mrSges(5,2) * t120;
t125 = Ifges(6,5) * t83 + Ifges(6,6) * t81;
t85 = t127 * t101;
t88 = t127 * t103;
t65 = -t106 * t85 + t109 * t88;
t124 = t101 ^ 2 + t103 ^ 2;
t123 = Ifges(5,4) * t101;
t122 = Ifges(5,4) * t103;
t119 = t102 ^ 2 + t100;
t26 = -t105 * t46 - t108 * t45;
t27 = -t105 * t45 + t108 * t46;
t118 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t82;
t93 = -pkin(4) * t103 - pkin(3);
t29 = t45 * mrSges(6,1) + t46 * mrSges(6,2);
t59 = -t81 * mrSges(6,1) + t83 * mrSges(6,2);
t10 = -t26 * mrSges(7,1) + t27 * mrSges(7,2);
t30 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t117 = -t104 * mrSges(3,1) + t102 * mrSges(3,2);
t86 = -t103 * mrSges(5,1) + t101 * mrSges(5,2);
t6 = -t106 * t22 + t109 * t13;
t62 = -t106 * t88 - t109 * t85;
t39 = pkin(4) * t121 + t63;
t40 = -pkin(9) * t83 + t62;
t41 = pkin(9) * t81 + t65;
t20 = -t105 * t41 + t108 * t40;
t21 = t105 * t40 + t108 * t41;
t116 = t20 * mrSges(7,1) - t21 * mrSges(7,2) + t126;
t4 = pkin(5) * t82 - pkin(9) * t46 + t6;
t5 = -pkin(9) * t45 + t7;
t2 = -t105 * t5 + t108 * t4;
t3 = t105 * t4 + t108 * t5;
t115 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t118;
t114 = (mrSges(7,1) * t108 - mrSges(7,2) * t105) * pkin(5);
t113 = -t30 - t59;
t91 = Ifges(5,1) * t101 + t122;
t90 = Ifges(5,2) * t103 + t123;
t72 = t84 * mrSges(4,2);
t67 = -pkin(5) * t81 + t93;
t61 = Ifges(6,1) * t83 + Ifges(6,4) * t81;
t60 = Ifges(6,4) * t83 + Ifges(6,2) * t81;
t57 = mrSges(5,1) * t82 - mrSges(5,3) * t120;
t56 = -mrSges(5,2) * t82 - mrSges(5,3) * t121;
t38 = t130 + (Ifges(5,1) * t103 - t123) * t84;
t37 = t129 + (-Ifges(5,2) * t101 + t122) * t84;
t36 = mrSges(6,1) * t82 - mrSges(6,3) * t46;
t35 = -mrSges(6,2) * t82 - mrSges(6,3) * t45;
t32 = Ifges(7,1) * t55 + Ifges(7,4) * t54;
t31 = Ifges(7,4) * t55 + Ifges(7,2) * t54;
t28 = pkin(5) * t45 + t39;
t17 = Ifges(6,1) * t46 - Ifges(6,4) * t45 + Ifges(6,5) * t82;
t16 = Ifges(6,4) * t46 - Ifges(6,2) * t45 + Ifges(6,6) * t82;
t15 = mrSges(7,1) * t82 - mrSges(7,3) * t27;
t14 = -mrSges(7,2) * t82 + mrSges(7,3) * t26;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t82;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t82;
t1 = [(Ifges(3,1) * t102 + 0.2e1 * Ifges(3,4) * t104) * t102 + 0.2e1 * t34 * t56 + 0.2e1 * t33 * t57 - t45 * t16 + t46 * t17 + (mrSges(4,1) * t133 - 0.2e1 * t66 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t82 + (Ifges(5,5) * t103 - Ifges(5,6) * t101 - (2 * Ifges(4,4))) * t84 + t118 + t136) * t82 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t36 + 0.2e1 * t39 * t29 + t26 * t8 + t27 * t9 + 0.2e1 * t28 * t10 + 0.2e1 * t3 * t14 + 0.2e1 * t2 * t15 + 0.2e1 * t119 * qJ(2) * mrSges(3,3) + t72 * t133 + t49 * t134 + m(7) * (t2 ^ 2 + t28 ^ 2 + t3 ^ 2) + m(6) * (t39 ^ 2 + t6 ^ 2 + t7 ^ 2) + Ifges(2,3) - 0.2e1 * pkin(1) * t117 + m(3) * (qJ(2) ^ 2 * t119 + pkin(1) ^ 2) + Ifges(3,2) * t100 + (mrSges(4,3) * t134 + Ifges(4,1) * t84 - t101 * t37 + t103 * t38) * t84 + m(4) * (t66 ^ 2 + t94 ^ 2 + t135) + m(5) * (t33 ^ 2 + t34 ^ 2 + t135); -m(3) * pkin(1) + t82 * mrSges(4,1) + t101 * t56 + t103 * t57 + t55 * t14 + t54 * t15 + t83 * t35 + t81 * t36 + t72 + m(7) * (t2 * t54 + t3 * t55) + m(6) * (t6 * t81 + t7 * t83) + m(5) * (t101 * t34 + t103 * t33) + m(4) * t94 + t117; m(3) + m(4) + m(5) * t124 + m(6) * (t81 ^ 2 + t83 ^ 2) + m(7) * (t54 ^ 2 + t55 ^ 2); (t86 - mrSges(4,1)) * t63 + m(6) * (t39 * t93 + t6 * t62 + t65 * t7) + m(7) * (t2 * t20 + t21 * t3 + t28 * t67) + t93 * t29 + t81 * t16 / 0.2e1 - Ifges(4,6) * t82 + t83 * t17 / 0.2e1 + Ifges(4,5) * t84 + t55 * t9 / 0.2e1 + t39 * t59 - t45 * t60 / 0.2e1 + t46 * t61 / 0.2e1 + t62 * t36 + t65 * t35 - t66 * mrSges(4,2) + t67 * t10 - pkin(3) * t49 + t54 * t8 / 0.2e1 + (t126 + t125) * t82 / 0.2e1 + t27 * t32 / 0.2e1 + t28 * t30 + t26 * t31 / 0.2e1 + t20 * t15 + t21 * t14 + m(5) * (-pkin(3) * t63 + (-t33 * t101 + t34 * t103) * qJ(4)) + (t84 * t91 / 0.2e1 + qJ(4) * t56 + t34 * mrSges(5,3) + t129 / 0.2e1 + t37 / 0.2e1) * t103 + (-t84 * t90 / 0.2e1 - qJ(4) * t57 - t33 * mrSges(5,3) + t38 / 0.2e1 + t130 / 0.2e1) * t101 + (-t6 * t83 + t7 * t81) * mrSges(6,3) + (-t2 * t55 + t3 * t54) * mrSges(7,3); m(6) * (t62 * t81 + t65 * t83) + m(7) * (t20 * t54 + t21 * t55); -0.2e1 * pkin(3) * t86 + t101 * t91 + t103 * t90 + 0.2e1 * t67 * t30 + t54 * t31 + t55 * t32 + 0.2e1 * t93 * t59 + t81 * t60 + t83 * t61 + Ifges(4,3) + m(7) * (t20 ^ 2 + t21 ^ 2 + t67 ^ 2) + m(6) * (t62 ^ 2 + t65 ^ 2 + t93 ^ 2) + m(5) * (qJ(4) ^ 2 * t124 + pkin(3) ^ 2) + 0.2e1 * (-t20 * t55 + t21 * t54) * mrSges(7,3) + 0.2e1 * (-t62 * t83 + t65 * t81) * mrSges(6,3) + 0.2e1 * t124 * qJ(4) * mrSges(5,3); m(5) * t63 + m(6) * t39 + m(7) * t28 + t10 + t29 + t49; 0; -m(5) * pkin(3) + m(6) * t93 + m(7) * t67 - t113 + t86; m(5) + m(6) + m(7); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (t105 * t14 + t108 * t15 + m(7) * (t105 * t3 + t108 * t2)) * pkin(5) + t115 + t136; m(7) * (t105 * t55 + t108 * t54) * pkin(5) + t113; t62 * mrSges(6,1) - t65 * mrSges(6,2) + (m(7) * (t105 * t21 + t108 * t20) + (t105 * t54 - t108 * t55) * mrSges(7,3)) * pkin(5) + t116 + t125; 0; Ifges(6,3) + Ifges(7,3) + m(7) * (t105 ^ 2 + t108 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t114; t115; -t30; t116; 0; Ifges(7,3) + t114; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
