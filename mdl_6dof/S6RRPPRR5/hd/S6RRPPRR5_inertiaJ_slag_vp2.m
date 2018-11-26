% Calculate joint inertia matrix for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2018-11-23 16:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:51:48
% EndTime: 2018-11-23 16:51:49
% DurationCPUTime: 1.11s
% Computational Cost: add. (1305->286), mult. (2712->387), div. (0->0), fcn. (2549->8), ass. (0->106)
t143 = Ifges(4,4) + Ifges(3,5);
t142 = Ifges(4,2) + Ifges(3,3) + Ifges(5,3);
t102 = cos(qJ(6));
t104 = cos(qJ(2));
t96 = sin(pkin(6));
t118 = t104 * t96;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t101 = sin(qJ(2));
t121 = t101 * t96;
t97 = cos(pkin(6));
t44 = -t100 * t97 + t103 * t121;
t99 = sin(qJ(6));
t26 = t102 * t118 - t44 * t99;
t43 = t100 * t121 + t103 * t97;
t12 = -mrSges(7,2) * t43 + mrSges(7,3) * t26;
t27 = t102 * t44 + t99 * t118;
t13 = mrSges(7,1) * t43 - mrSges(7,3) * t27;
t141 = t102 * t12 - t99 * t13;
t35 = -t96 * pkin(1) - pkin(2) * t118 - qJ(3) * t121;
t140 = -0.2e1 * t35;
t139 = t26 / 0.2e1;
t138 = t27 / 0.2e1;
t57 = Ifges(7,5) * t99 + Ifges(7,6) * t102;
t137 = t57 / 0.2e1;
t136 = t99 / 0.2e1;
t105 = -pkin(2) - pkin(3);
t135 = pkin(1) * t97;
t134 = -t102 / 0.2e1;
t133 = t102 / 0.2e1;
t132 = Ifges(7,4) * t99;
t45 = -pkin(8) * t121 + t104 * t135;
t131 = t45 * mrSges(3,1);
t46 = pkin(8) * t118 + t101 * t135;
t130 = t46 * mrSges(3,2);
t128 = -mrSges(5,1) - mrSges(4,1);
t29 = pkin(3) * t118 - t35;
t17 = (pkin(4) * t104 - pkin(9) * t101) * t96 + t29;
t34 = t97 * qJ(3) + t46;
t28 = -qJ(4) * t118 + t34;
t22 = -pkin(9) * t97 + t28;
t6 = t100 * t17 + t103 * t22;
t11 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t31 = mrSges(6,1) * t118 - mrSges(6,3) * t44;
t127 = -t31 + t11;
t126 = mrSges(5,1) * t118 + mrSges(5,2) * t121;
t125 = t102 ^ 2 + t99 ^ 2;
t93 = t100 ^ 2;
t95 = t103 ^ 2;
t124 = t93 + t95;
t123 = Ifges(7,4) * t102;
t122 = t100 * t99;
t98 = qJ(3) - pkin(9);
t119 = t103 * t98;
t117 = t100 * t102;
t91 = pkin(4) - t105;
t7 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t43;
t116 = Ifges(6,5) * t44 - Ifges(6,6) * t43 + Ifges(6,3) * t118;
t115 = mrSges(5,3) * t121;
t114 = -t118 / 0.2e1;
t113 = t125 * t100;
t56 = t103 * mrSges(6,1) - t100 * mrSges(6,2);
t36 = -t97 * pkin(2) - t45;
t23 = -t97 * pkin(3) - qJ(4) * t121 + t36;
t18 = -t97 * pkin(4) + t23;
t10 = pkin(5) * t43 - pkin(10) * t44 - t18;
t4 = pkin(10) * t118 + t6;
t1 = t10 * t102 - t4 * t99;
t2 = t10 * t99 + t102 * t4;
t112 = -t1 * t99 + t102 * t2;
t20 = t43 * mrSges(6,1) + t44 * mrSges(6,2);
t111 = -mrSges(7,1) * t99 - mrSges(7,2) * t102;
t89 = t103 * pkin(5);
t54 = pkin(10) * t100 + t89 + t91;
t32 = t102 * t54 - t99 * t119;
t33 = t102 * t119 + t54 * t99;
t110 = t102 * t33 - t32 * t99;
t52 = -mrSges(7,2) * t103 + mrSges(7,3) * t122;
t53 = mrSges(7,1) * t103 + mrSges(7,3) * t117;
t109 = t102 * t52 - t99 * t53;
t5 = -t100 * t22 + t103 * t17;
t37 = -Ifges(7,5) * t117 + Ifges(7,6) * t122 + Ifges(7,3) * t103;
t108 = Ifges(3,6) * t118 + t121 * t143 + t142 * t97;
t106 = qJ(3) ^ 2;
t90 = t98 ^ 2;
t77 = t93 * t90;
t65 = mrSges(4,2) * t121;
t61 = -Ifges(6,1) * t100 - Ifges(6,4) * t103;
t60 = Ifges(7,1) * t99 + t123;
t59 = -Ifges(6,4) * t100 - Ifges(6,2) * t103;
t58 = Ifges(7,2) * t102 + t132;
t55 = -mrSges(7,1) * t102 + mrSges(7,2) * t99;
t51 = mrSges(4,2) * t118 + mrSges(4,3) * t97;
t50 = mrSges(5,2) * t97 - mrSges(5,3) * t118;
t49 = -t97 * mrSges(4,1) + t65;
t48 = -t97 * mrSges(5,1) - t115;
t47 = t111 * t100;
t39 = Ifges(7,5) * t103 + (-Ifges(7,1) * t102 + t132) * t100;
t38 = Ifges(7,6) * t103 + (Ifges(7,2) * t99 - t123) * t100;
t30 = -mrSges(6,2) * t118 - mrSges(6,3) * t43;
t16 = Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t118;
t15 = Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t118;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t43;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t43;
t3 = -pkin(5) * t118 - t5;
t14 = [0.2e1 * t36 * t49 + 0.2e1 * t28 * t50 + 0.2e1 * t34 * t51 + 0.2e1 * t29 * t126 + t44 * t16 + 0.2e1 * t23 * t48 + 0.2e1 * t6 * t30 + 0.2e1 * t5 * t31 - 0.2e1 * t18 * t20 + t26 * t8 + t27 * t9 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + Ifges(2,3) + (t7 - t15) * t43 + (t108 - 0.2e1 * t130 + 0.2e1 * t131) * t97 + m(3) * (t45 ^ 2 + t46 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t18 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (m(3) * pkin(1) ^ 2 * t96 + (-0.2e1 * t45 * mrSges(3,3) + mrSges(4,3) * t140 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,1) + Ifges(3,1) + Ifges(5,1)) * t101) * t96 + (-(2 * Ifges(5,5)) + t143) * t97) * t101 + (mrSges(4,1) * t140 + 0.2e1 * t46 * mrSges(3,3) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2) + Ifges(5,2)) * t104) * t96 + 0.2e1 * (Ifges(3,4) - Ifges(5,4) - Ifges(4,5)) * t121 + (-(2 * Ifges(4,6)) + Ifges(3,6) + (2 * Ifges(5,6))) * t97 + t116) * t104) * t96; (-t59 / 0.2e1 + t37 / 0.2e1) * t43 + t108 + (t98 * t30 - t6 * mrSges(6,3) + Ifges(6,6) * t114 + t7 / 0.2e1 - t15 / 0.2e1) * t103 + t105 * t48 + t91 * t20 - pkin(2) * t49 + t2 * t52 + t1 * t53 - t18 * t56 + t44 * t61 / 0.2e1 + t3 * t47 + t32 * t13 + t33 * t12 + t34 * mrSges(4,3) - t36 * mrSges(4,1) + (t9 * t134 + Ifges(6,5) * t114 + t8 * t136 + t5 * mrSges(6,3) - t16 / 0.2e1 + t127 * t98) * t100 - t23 * mrSges(5,1) + t28 * mrSges(5,2) + m(6) * (-t18 * t91 + (-t100 * t5 + t103 * t6) * t98) - t130 + t131 + t39 * t138 + t38 * t139 + m(5) * (qJ(3) * t28 + t105 * t23) + m(7) * (t100 * t3 * t98 + t1 * t32 + t2 * t33) + m(4) * (-pkin(2) * t36 + qJ(3) * t34) + (t50 + t51) * qJ(3) + (-Ifges(5,5) * t101 + (-Ifges(4,6) + Ifges(5,6)) * t104) * t96; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t105 * mrSges(5,1) + 0.2e1 * t32 * t53 + 0.2e1 * t33 * t52 + 0.2e1 * t91 * t56 + (t37 - t59) * t103 + (-t102 * t39 + t38 * t99 + 0.2e1 * t47 * t98 - t61) * t100 + m(7) * (t32 ^ 2 + t33 ^ 2 + t77) + m(6) * (t90 * t95 + t91 ^ 2 + t77) + m(5) * (t105 ^ 2 + t106) + m(4) * (pkin(2) ^ 2 + t106) + 0.2e1 * (mrSges(4,3) + mrSges(5,2)) * qJ(3) - 0.2e1 * t124 * t98 * mrSges(6,3) + t142; -t115 - t102 * t13 - t99 * t12 + t65 + t128 * t97 + m(7) * (-t1 * t102 - t2 * t99) + m(6) * t18 + m(5) * t23 + m(4) * t36 - t20; -m(4) * pkin(2) - t102 * t53 - t99 * t52 + m(7) * (-t102 * t32 - t33 * t99) - m(6) * t91 + m(5) * t105 - t56 + t128; m(7) * t125 + m(4) + m(5) + m(6); -t127 * t103 + (t30 + t141) * t100 + m(7) * (t112 * t100 - t103 * t3) + m(6) * (t100 * t6 + t103 * t5) + m(5) * t29 + t126; -t103 * t47 + (m(7) * (t110 - t119) + t109) * t100; 0; m(5) + m(6) * t124 + m(7) * (t125 * t93 + t95); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t112 * mrSges(7,3) + t8 * t133 + t9 * t136 + t43 * t137 + t60 * t138 + t58 * t139 + t3 * t55 + t116 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t112 + t141) * pkin(10); t39 * t136 + t38 * t133 - pkin(5) * t47 + t110 * mrSges(7,3) + (m(7) * t110 + t109) * pkin(10) + (t58 * t136 + t60 * t134 - Ifges(6,5) + (-m(7) * pkin(5) - mrSges(6,1) + t55) * t98) * t100 + (-t98 * mrSges(6,2) - Ifges(6,6) + t137) * t103; 0; -t103 * t55 + m(7) * (pkin(10) * t113 + t89) + mrSges(7,3) * t113 + t56; Ifges(6,3) - 0.2e1 * pkin(5) * t55 + t99 * t60 + t102 * t58 + m(7) * (t125 * pkin(10) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t125 * pkin(10) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t32 - mrSges(7,2) * t33 + t37; t55; t47; t111 * pkin(10) + t57; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
