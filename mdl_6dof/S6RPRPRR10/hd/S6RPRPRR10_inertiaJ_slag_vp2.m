% Calculate joint inertia matrix for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:08:55
% EndTime: 2018-11-23 16:08:56
% DurationCPUTime: 0.82s
% Computational Cost: add. (1586->267), mult. (3104->381), div. (0->0), fcn. (3197->8), ass. (0->107)
t105 = (-pkin(1) - pkin(7));
t135 = -2 * t105;
t102 = cos(qJ(6));
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t97 = sin(pkin(10));
t98 = cos(pkin(10));
t72 = -t100 * t97 + t103 * t98;
t73 = t100 * t98 + t103 * t97;
t99 = sin(qJ(6));
t39 = t102 * t72 - t73 * t99;
t40 = t102 * t73 + t72 * t99;
t14 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t43 = -t72 * mrSges(6,1) + t73 * mrSges(6,2);
t134 = -t14 - t43;
t104 = cos(qJ(3));
t59 = t73 * t104;
t61 = t72 * t104;
t31 = t59 * mrSges(6,1) + t61 * mrSges(6,2);
t121 = t104 * t98;
t122 = t104 * t97;
t62 = mrSges(5,1) * t122 + mrSges(5,2) * t121;
t26 = -t102 * t59 - t61 * t99;
t28 = t102 * t61 - t59 * t99;
t9 = -t26 * mrSges(7,1) + t28 * mrSges(7,2);
t133 = -t31 - t62 - t9;
t101 = sin(qJ(3));
t132 = Ifges(6,5) * t61 - Ifges(6,6) * t59 + Ifges(6,3) * t101;
t131 = 2 * qJ(2);
t130 = t97 / 0.2e1;
t129 = t98 / 0.2e1;
t128 = Ifges(5,4) * t97;
t127 = Ifges(5,4) * t98;
t126 = pkin(8) + qJ(4);
t76 = pkin(3) * t101 - qJ(4) * t104 + qJ(2);
t67 = t98 * t76;
t37 = -pkin(8) * t121 + t67 + (-t105 * t97 + pkin(4)) * t101;
t120 = t101 * t105;
t51 = t98 * t120 + t97 * t76;
t42 = -pkin(8) * t122 + t51;
t13 = t100 * t37 + t103 * t42;
t77 = t126 * t97;
t79 = t126 * t98;
t47 = -t100 * t77 + t103 * t79;
t78 = -t98 * mrSges(5,1) + t97 * mrSges(5,2);
t125 = -t78 + mrSges(4,1);
t124 = t97 ^ 2 + t98 ^ 2;
t94 = t101 ^ 2;
t95 = t104 ^ 2;
t123 = t95 + t94;
t119 = t104 * t105;
t118 = m(5) + m(6) + m(7);
t117 = Ifges(7,5) * t28 + Ifges(7,6) * t26 + Ifges(7,3) * t101;
t86 = -pkin(4) * t98 - pkin(3);
t116 = t123 * mrSges(4,3);
t58 = t73 * t101;
t60 = t72 * t101;
t25 = -t102 * t58 - t60 * t99;
t27 = t102 * t60 - t58 * t99;
t115 = t25 * mrSges(7,1) - t27 * mrSges(7,2);
t114 = qJ(4) * t124;
t12 = -t100 * t42 + t103 * t37;
t46 = -t100 * t79 - t103 * t77;
t68 = pkin(4) * t122 - t119;
t50 = -t120 * t97 + t67;
t113 = -t50 * t97 + t51 * t98;
t74 = -mrSges(5,2) * t101 - mrSges(5,3) * t122;
t75 = mrSges(5,1) * t101 - mrSges(5,3) * t121;
t112 = t98 * t74 - t97 * t75;
t29 = -pkin(9) * t73 + t46;
t30 = pkin(9) * t72 + t47;
t10 = t102 * t29 - t30 * t99;
t11 = t102 * t30 + t29 * t99;
t35 = Ifges(7,6) * t39;
t36 = Ifges(7,5) * t40;
t111 = t10 * mrSges(7,1) - t11 * mrSges(7,2) + t35 + t36;
t4 = pkin(5) * t101 - pkin(9) * t61 + t12;
t8 = -pkin(9) * t59 + t13;
t2 = t102 * t4 - t8 * t99;
t3 = t102 * t8 + t4 * t99;
t110 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t117;
t109 = (mrSges(7,1) * t102 - mrSges(7,2) * t99) * pkin(5);
t107 = qJ(2) ^ 2;
t96 = t105 ^ 2;
t88 = t95 * t105;
t87 = t95 * t96;
t81 = Ifges(5,1) * t97 + t127;
t80 = Ifges(5,2) * t98 + t128;
t65 = Ifges(6,5) * t73;
t64 = Ifges(6,6) * t72;
t57 = Ifges(5,5) * t101 + (Ifges(5,1) * t98 - t128) * t104;
t56 = Ifges(5,6) * t101 + (-Ifges(5,2) * t97 + t127) * t104;
t52 = -pkin(5) * t72 + t86;
t49 = mrSges(6,1) * t101 - mrSges(6,3) * t61;
t48 = -mrSges(6,2) * t101 - mrSges(6,3) * t59;
t45 = Ifges(6,1) * t73 + Ifges(6,4) * t72;
t44 = Ifges(6,4) * t73 + Ifges(6,2) * t72;
t41 = pkin(5) * t59 + t68;
t24 = Ifges(6,1) * t61 - Ifges(6,4) * t59 + Ifges(6,5) * t101;
t23 = Ifges(6,4) * t61 - Ifges(6,2) * t59 + Ifges(6,6) * t101;
t18 = mrSges(7,1) * t101 - mrSges(7,3) * t28;
t17 = -mrSges(7,2) * t101 + mrSges(7,3) * t26;
t16 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t15 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t6 = Ifges(7,1) * t28 + Ifges(7,4) * t26 + Ifges(7,5) * t101;
t5 = Ifges(7,4) * t28 + Ifges(7,2) * t26 + Ifges(7,6) * t101;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t131) + 0.2e1 * t12 * t49 + 0.2e1 * t13 * t48 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t59 * t23 + t61 * t24 + t26 * t5 + t28 * t6 + 0.2e1 * t68 * t31 + 0.2e1 * t41 * t9 + 0.2e1 * t50 * t75 + 0.2e1 * t51 * t74 + Ifges(3,1) + Ifges(2,3) + t116 * t135 + ((mrSges(4,2) * t131) + Ifges(4,1) * t104 + t62 * t135 - t97 * t56 + t98 * t57) * t104 + (mrSges(4,1) * t131 + (Ifges(5,3) + Ifges(4,2)) * t101 + (Ifges(5,5) * t98 - Ifges(5,6) * t97 - (2 * Ifges(4,4))) * t104 + t117 + t132) * t101 + m(4) * (t94 * t96 + t107 + t87) + (m(3) * (pkin(1) ^ 2 + t107)) + m(5) * (t50 ^ 2 + t51 ^ 2 + t87) + m(6) * (t12 ^ 2 + t13 ^ 2 + t68 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t41 ^ 2); -(m(3) * pkin(1)) + t27 * t17 + t25 * t18 + t60 * t48 - t58 * t49 + mrSges(3,2) + t112 * t101 - t116 + t133 * t104 + m(7) * (-t104 * t41 + t2 * t25 + t27 * t3) + m(6) * (-t104 * t68 - t12 * t58 + t13 * t60) + m(5) * (t101 * t113 + t88) + m(4) * (t105 * t94 + t88); m(3) + m(7) * (t25 ^ 2 + t27 ^ 2 + t95) + m(6) * (t58 ^ 2 + t60 ^ 2 + t95) + m(5) * (t124 * t94 + t95) + m(4) * t123; m(6) * (t12 * t46 + t13 * t47 + t68 * t86) + m(7) * (t10 * t2 + t11 * t3 + t41 * t52) + (-t2 * t40 + t3 * t39) * mrSges(7,3) + (-t12 * t73 + t13 * t72) * mrSges(6,3) + t56 * t129 + t57 * t130 + t86 * t31 + t68 * t43 + t72 * t23 / 0.2e1 + t73 * t24 / 0.2e1 - t59 * t44 / 0.2e1 + t61 * t45 / 0.2e1 - pkin(3) * t62 + t47 * t48 + t46 * t49 + t52 * t9 + t39 * t5 / 0.2e1 + t40 * t6 / 0.2e1 + t41 * t14 + t26 * t15 / 0.2e1 + t28 * t16 / 0.2e1 + t11 * t17 + t10 * t18 + (Ifges(4,5) - t97 * t80 / 0.2e1 + t81 * t129 + t125 * t105) * t104 + (Ifges(5,5) * t130 + Ifges(5,6) * t129 + t65 / 0.2e1 + t64 / 0.2e1 + t36 / 0.2e1 + t35 / 0.2e1 - Ifges(4,6) - (t105 * mrSges(4,2))) * t101 + t112 * qJ(4) + t113 * mrSges(5,3) + m(5) * (pkin(3) * t119 + qJ(4) * t113); (-t25 * t40 + t27 * t39) * mrSges(7,3) + (t58 * t73 + t60 * t72) * mrSges(6,3) + (mrSges(5,3) * t124 - mrSges(4,2)) * t101 + (t125 + t134) * t104 + m(7) * (t10 * t25 - t104 * t52 + t11 * t27) + m(6) * (-t104 * t86 - t46 * t58 + t47 * t60) + m(5) * (pkin(3) * t104 + t101 * t114); -0.2e1 * pkin(3) * t78 + 0.2e1 * t52 * t14 + t39 * t15 + t40 * t16 + 0.2e1 * t86 * t43 + t72 * t44 + t73 * t45 + t98 * t80 + t97 * t81 + Ifges(4,3) + m(7) * (t10 ^ 2 + t11 ^ 2 + t52 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t86 ^ 2) + m(5) * (qJ(4) ^ 2 * t124 + pkin(3) ^ 2) + 0.2e1 * (-t10 * t40 + t11 * t39) * mrSges(7,3) + 0.2e1 * (-t46 * t73 + t47 * t72) * mrSges(6,3) + 0.2e1 * mrSges(5,3) * t114; -m(5) * t119 + m(6) * t68 + m(7) * t41 - t133; -t118 * t104; -m(5) * pkin(3) + m(6) * t86 + m(7) * t52 - t134 + t78; t118; t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t102 * t2 + t3 * t99) + t99 * t17 + t102 * t18) * pkin(5) + t110 + t132; -t58 * mrSges(6,1) - t60 * mrSges(6,2) + m(7) * (t102 * t25 + t27 * t99) * pkin(5) + t115; t46 * mrSges(6,1) - t47 * mrSges(6,2) + t64 + t65 + (m(7) * (t10 * t102 + t11 * t99) + (-t102 * t40 + t39 * t99) * mrSges(7,3)) * pkin(5) + t111; 0; Ifges(6,3) + Ifges(7,3) + m(7) * (t102 ^ 2 + t99 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t109; t110; t115; t111; 0; Ifges(7,3) + t109; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
