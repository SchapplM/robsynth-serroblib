% Calculate joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:14:34
% EndTime: 2018-11-23 15:14:35
% DurationCPUTime: 0.77s
% Computational Cost: add. (1617->219), mult. (3340->323), div. (0->0), fcn. (3833->12), ass. (0->92)
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t60 = -mrSges(7,1) * t84 + mrSges(7,2) * t80;
t130 = -mrSges(6,1) + t60;
t129 = m(5) * pkin(3);
t120 = cos(qJ(5));
t76 = sin(pkin(12));
t78 = cos(pkin(12));
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t54 = -t76 * t82 + t78 * t85;
t55 = t76 * t85 + t78 * t82;
t81 = sin(qJ(5));
t39 = -t120 * t54 + t55 * t81;
t40 = t120 * t55 + t81 * t54;
t25 = t39 * mrSges(6,1) + t40 * mrSges(6,2);
t41 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t128 = -t25 - t41;
t109 = t80 * mrSges(7,3);
t23 = -mrSges(7,2) * t39 - t40 * t109;
t114 = t40 * t84;
t24 = mrSges(7,1) * t39 - mrSges(7,3) * t114;
t127 = t84 * t23 - t80 * t24;
t106 = -qJ(4) - pkin(8);
t100 = t106 * t85;
t99 = t106 * t82;
t43 = -t78 * t100 + t76 * t99;
t31 = pkin(9) * t54 + t43;
t42 = t76 * t100 + t78 * t99;
t91 = -t55 * pkin(9) + t42;
t15 = -t120 * t91 + t31 * t81;
t126 = t15 ^ 2;
t77 = sin(pkin(6));
t83 = sin(qJ(2));
t111 = t77 * t83;
t79 = cos(pkin(6));
t51 = -t82 * t111 + t79 * t85;
t52 = t85 * t111 + t79 * t82;
t32 = t51 * t78 - t52 * t76;
t33 = t51 * t76 + t52 * t78;
t19 = -t120 * t32 + t33 * t81;
t125 = t19 ^ 2;
t75 = t85 ^ 2;
t124 = 0.2e1 * t15;
t123 = -m(5) - m(6);
t122 = pkin(3) * t76;
t68 = -pkin(3) * t85 - pkin(2);
t44 = -pkin(4) * t54 + t68;
t14 = pkin(5) * t39 - pkin(10) * t40 + t44;
t17 = t120 * t31 + t81 * t91;
t3 = t14 * t80 + t17 * t84;
t121 = t3 * t84;
t119 = Ifges(7,4) * t80;
t118 = Ifges(7,4) * t84;
t86 = cos(qJ(2));
t110 = t77 * t86;
t21 = t120 * t33 + t81 * t32;
t11 = -t80 * t110 + t21 * t84;
t117 = t11 * t84;
t116 = t15 * t19;
t115 = t40 * t80;
t67 = pkin(3) * t78 + pkin(4);
t48 = t120 * t67 - t81 * t122;
t113 = t48 * mrSges(6,1);
t49 = t120 * t122 + t81 * t67;
t112 = t49 * mrSges(6,2);
t105 = Ifges(7,5) * t114 + Ifges(7,3) * t39;
t104 = Ifges(7,5) * t80 + Ifges(7,6) * t84;
t103 = t80 ^ 2 + t84 ^ 2;
t102 = t82 ^ 2 + t75;
t62 = Ifges(7,2) * t84 + t119;
t63 = Ifges(7,1) * t80 + t118;
t101 = t84 * t62 + t80 * t63 + Ifges(6,3);
t47 = pkin(10) + t49;
t98 = t103 * t47;
t2 = t14 * t84 - t17 * t80;
t97 = -t2 * t80 + t121;
t96 = mrSges(7,1) * t80 + mrSges(7,2) * t84;
t10 = -t84 * t110 - t21 * t80;
t95 = -t10 * t80 + t117;
t94 = -t51 * t82 + t52 * t85;
t93 = 0.2e1 * t103 * mrSges(7,3);
t92 = -t21 * mrSges(6,2) + mrSges(7,3) * t117 - t10 * t109 + t130 * t19;
t6 = Ifges(7,6) * t39 + (-Ifges(7,2) * t80 + t118) * t40;
t7 = Ifges(7,5) * t39 + (Ifges(7,1) * t84 - t119) * t40;
t90 = -t17 * mrSges(6,2) + mrSges(7,3) * t121 - t2 * t109 - t62 * t115 / 0.2e1 + t63 * t114 / 0.2e1 + Ifges(6,5) * t40 + t80 * t7 / 0.2e1 + t84 * t6 / 0.2e1 + (t104 / 0.2e1 - Ifges(6,6)) * t39 + t130 * t15;
t71 = t77 ^ 2;
t66 = t71 * t86 ^ 2;
t61 = -mrSges(4,1) * t85 + mrSges(4,2) * t82;
t46 = -pkin(5) - t48;
t22 = t96 * t40;
t1 = [m(2) + m(7) * (t10 ^ 2 + t11 ^ 2 + t125) + m(6) * (t21 ^ 2 + t125 + t66) + m(5) * (t32 ^ 2 + t33 ^ 2 + t66) + m(4) * (t51 ^ 2 + t52 ^ 2 + t66) + m(3) * (t71 * t83 ^ 2 + t79 ^ 2 + t66); t10 * t24 + t11 * t23 + t19 * t22 + (t19 * t40 - t21 * t39) * mrSges(6,3) + (-t32 * t55 + t33 * t54) * mrSges(5,3) + t94 * mrSges(4,3) + (-t83 * mrSges(3,2) + (mrSges(3,1) - t61 + t128) * t86) * t77 + m(7) * (t10 * t2 + t11 * t3 + t116) + m(6) * (-t44 * t110 + t17 * t21 + t116) + m(5) * (-t68 * t110 + t32 * t42 + t33 * t43) + m(4) * (pkin(2) * t110 + t94 * pkin(8)); Ifges(4,2) * t75 - 0.2e1 * pkin(2) * t61 + t22 * t124 + 0.2e1 * t2 * t24 + 0.2e1 * t3 * t23 + 0.2e1 * t44 * t25 + 0.2e1 * t68 * t41 + Ifges(3,3) + (Ifges(4,1) * t82 + 0.2e1 * Ifges(4,4) * t85) * t82 + (-0.2e1 * mrSges(5,3) * t42 + Ifges(5,1) * t55) * t55 + 0.2e1 * t102 * pkin(8) * mrSges(4,3) + (0.2e1 * mrSges(5,3) * t43 + 0.2e1 * Ifges(5,4) * t55 + Ifges(5,2) * t54) * t54 + (-0.2e1 * mrSges(6,3) * t17 + Ifges(6,2) * t39 + t105) * t39 + (mrSges(6,3) * t124 + Ifges(6,1) * t40 - t80 * t6 + t84 * t7 + (-Ifges(7,6) * t80 - (2 * Ifges(6,4))) * t39) * t40 + m(4) * (t102 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2 + t68 ^ 2) + m(6) * (t17 ^ 2 + t44 ^ 2 + t126) + m(7) * (t2 ^ 2 + t3 ^ 2 + t126); t51 * mrSges(4,1) + t32 * mrSges(5,1) - t52 * mrSges(4,2) - t33 * mrSges(5,2) + m(7) * (t19 * t46 + t95 * t47) + m(6) * (-t19 * t48 + t21 * t49) + (t32 * t78 + t33 * t76) * t129 + t92; (-t82 * mrSges(4,1) - t85 * mrSges(4,2)) * pkin(8) + (-t49 * t39 - t48 * t40) * mrSges(6,3) + Ifges(4,5) * t82 + Ifges(4,6) * t85 + (m(5) * (t42 * t78 + t43 * t76) + (t76 * t54 - t78 * t55) * mrSges(5,3)) * pkin(3) + Ifges(5,6) * t54 + Ifges(5,5) * t55 + t42 * mrSges(5,1) - t43 * mrSges(5,2) + t46 * t22 + m(7) * (t46 * t15 + t97 * t47) + t90 + m(6) * (-t15 * t48 + t17 * t49) + t127 * t47; 0.2e1 * t113 - 0.2e1 * t112 + 0.2e1 * t46 * t60 + Ifges(4,3) + Ifges(5,3) + t47 * t93 + m(7) * (t103 * t47 ^ 2 + t46 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2) + t101 + (0.2e1 * mrSges(5,1) * t78 - 0.2e1 * mrSges(5,2) * t76 + (t76 ^ 2 + t78 ^ 2) * t129) * pkin(3); m(7) * (t10 * t84 + t11 * t80) + t123 * t110; t80 * t23 + t84 * t24 + m(7) * (t2 * t84 + t3 * t80) + m(6) * t44 + m(5) * t68 - t128; 0; m(7) * t103 - t123; m(7) * (-pkin(5) * t19 + t95 * pkin(10)) + t92; t90 + (-m(7) * t15 - t22) * pkin(5) + (m(7) * t97 + t127) * pkin(10); m(7) * (-pkin(5) * t46 + pkin(10) * t98) - t112 + t113 + (-pkin(5) + t46) * t60 + (t103 * pkin(10) + t98) * mrSges(7,3) + t101; 0; -0.2e1 * pkin(5) * t60 + m(7) * (t103 * pkin(10) ^ 2 + pkin(5) ^ 2) + pkin(10) * t93 + t101; mrSges(7,1) * t10 - mrSges(7,2) * t11; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t115 + t105; -t96 * t47 + t104; -t60; -t96 * pkin(10) + t104; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
