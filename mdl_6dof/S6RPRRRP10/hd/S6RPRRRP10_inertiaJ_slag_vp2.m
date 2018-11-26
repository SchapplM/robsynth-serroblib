% Calculate joint inertia matrix for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2018-11-23 16:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:29:48
% EndTime: 2018-11-23 16:29:49
% DurationCPUTime: 1.04s
% Computational Cost: add. (1223->286), mult. (2340->384), div. (0->0), fcn. (2147->6), ass. (0->101)
t134 = -mrSges(6,1) - mrSges(7,1);
t133 = mrSges(6,2) - mrSges(7,3);
t119 = Ifges(7,2) + Ifges(6,3);
t96 = (-pkin(1) - pkin(7));
t132 = -2 * t96;
t131 = mrSges(7,2) + mrSges(6,3);
t94 = cos(qJ(4));
t95 = cos(qJ(3));
t122 = t94 * t95;
t92 = sin(qJ(3));
t130 = Ifges(5,5) * t122 + Ifges(5,3) * t92;
t129 = 2 * mrSges(7,3);
t128 = 2 * qJ(2);
t127 = -pkin(9) - pkin(8);
t91 = sin(qJ(4));
t126 = Ifges(5,4) * t91;
t125 = Ifges(5,4) * t94;
t124 = t91 * t95;
t123 = t92 * t96;
t121 = t95 * t96;
t66 = pkin(3) * t92 - pkin(8) * t95 + qJ(2);
t59 = t94 * t66;
t17 = -pkin(9) * t122 + t59 + (-t91 * t96 + pkin(4)) * t92;
t37 = t123 * t94 + t66 * t91;
t19 = -pkin(9) * t124 + t37;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t6 = t17 * t90 + t19 * t93;
t62 = t90 * t94 + t91 * t93;
t50 = t62 * t95;
t32 = -mrSges(7,2) * t50 + mrSges(7,3) * t92;
t33 = -mrSges(6,2) * t92 - mrSges(6,3) * t50;
t118 = t32 + t33;
t61 = t90 * t91 - t93 * t94;
t52 = t61 * t95;
t34 = mrSges(6,1) * t92 + mrSges(6,3) * t52;
t35 = -mrSges(7,1) * t92 - mrSges(7,2) * t52;
t117 = -t34 + t35;
t67 = -mrSges(5,1) * t94 + mrSges(5,2) * t91;
t116 = -t67 + mrSges(4,1);
t115 = t91 ^ 2 + t94 ^ 2;
t86 = t92 ^ 2;
t88 = t95 ^ 2;
t114 = t88 + t86;
t112 = t127 * t91;
t70 = t127 * t94;
t29 = -t112 * t93 - t70 * t90;
t31 = t112 * t90 - t93 * t70;
t113 = t29 ^ 2 + t31 ^ 2;
t77 = -pkin(4) * t94 - pkin(3);
t111 = t115 * mrSges(5,3);
t110 = t114 * mrSges(4,3);
t48 = t62 * t92;
t51 = t61 * t92;
t109 = t29 * t48 - t31 * t51;
t60 = pkin(4) * t124 - t121;
t107 = mrSges(5,1) * t91 + mrSges(5,2) * t94;
t5 = t17 * t93 - t19 * t90;
t36 = -t123 * t91 + t59;
t106 = -t36 * t91 + t37 * t94;
t104 = t119 * t92 + (-Ifges(7,4) - Ifges(6,5)) * t52 + (-Ifges(6,6) + Ifges(7,6)) * t50;
t103 = (mrSges(6,1) * t93 - mrSges(6,2) * t90) * pkin(4);
t102 = t133 * t51 + t134 * t48;
t54 = Ifges(7,6) * t61;
t55 = Ifges(6,6) * t61;
t56 = Ifges(6,5) * t62;
t57 = Ifges(7,4) * t62;
t101 = -t133 * t31 + t134 * t29 + t54 - t55 + t56 + t57;
t2 = qJ(6) * t92 + t6;
t3 = -pkin(5) * t92 - t5;
t100 = mrSges(6,1) * t5 - t3 * mrSges(7,1) - t6 * mrSges(6,2) + mrSges(7,3) * t2 + t104;
t97 = qJ(2) ^ 2;
t89 = t96 ^ 2;
t84 = Ifges(5,5) * t91;
t82 = Ifges(5,6) * t94;
t79 = t88 * t96;
t78 = t88 * t89;
t76 = -pkin(4) * t93 - pkin(5);
t74 = pkin(4) * t90 + qJ(6);
t69 = Ifges(5,1) * t91 + t125;
t68 = Ifges(5,2) * t94 + t126;
t65 = mrSges(5,1) * t92 - mrSges(5,3) * t122;
t64 = -mrSges(5,2) * t92 - mrSges(5,3) * t124;
t53 = t107 * t95;
t46 = Ifges(5,5) * t92 + (Ifges(5,1) * t94 - t126) * t95;
t45 = Ifges(5,6) * t92 + (-Ifges(5,2) * t91 + t125) * t95;
t25 = Ifges(6,1) * t62 - Ifges(6,4) * t61;
t24 = Ifges(7,1) * t62 + Ifges(7,5) * t61;
t23 = Ifges(6,4) * t62 - Ifges(6,2) * t61;
t22 = Ifges(7,5) * t62 + Ifges(7,3) * t61;
t21 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t20 = mrSges(7,1) * t61 - mrSges(7,3) * t62;
t16 = pkin(5) * t61 - qJ(6) * t62 + t77;
t14 = mrSges(6,1) * t50 - mrSges(6,2) * t52;
t13 = mrSges(7,1) * t50 + mrSges(7,3) * t52;
t11 = -Ifges(6,1) * t52 - Ifges(6,4) * t50 + Ifges(6,5) * t92;
t10 = -Ifges(7,1) * t52 + Ifges(7,4) * t92 + Ifges(7,5) * t50;
t9 = -Ifges(6,4) * t52 - Ifges(6,2) * t50 + Ifges(6,6) * t92;
t8 = -Ifges(7,5) * t52 + Ifges(7,6) * t92 + Ifges(7,3) * t50;
t7 = pkin(5) * t50 + qJ(6) * t52 + t60;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t128) + 0.2e1 * t7 * t13 + 0.2e1 * t60 * t14 + 0.2e1 * t2 * t32 + 0.2e1 * t3 * t35 + 0.2e1 * t6 * t33 + 0.2e1 * t5 * t34 + 0.2e1 * t36 * t65 + 0.2e1 * t37 * t64 + Ifges(3,1) + Ifges(2,3) - (t10 + t11) * t52 + (t8 - t9) * t50 + t110 * t132 + ((mrSges(4,2) * t128) + Ifges(4,1) * t95 + t53 * t132 - t45 * t91 + t46 * t94) * t95 + (mrSges(4,1) * t128 + Ifges(4,2) * t92 + (-Ifges(5,6) * t91 - (2 * Ifges(4,4))) * t95 + t104 + t130) * t92 + m(4) * (t86 * t89 + t78 + t97) + (m(3) * (pkin(1) ^ 2 + t97)) + m(5) * (t36 ^ 2 + t37 ^ 2 + t78) + m(6) * (t5 ^ 2 + t6 ^ 2 + t60 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2); -(m(3) * pkin(1)) + mrSges(3,2) + (t94 * t64 - t91 * t65) * t92 - t118 * t51 + t117 * t48 - t110 + (-t13 - t14 - t53) * t95 + m(7) * (-t2 * t51 + t3 * t48 - t7 * t95) + m(6) * (-t48 * t5 - t51 * t6 - t60 * t95) + m(5) * (t106 * t92 + t79) + m(4) * (t86 * t96 + t79); m(3) + m(5) * (t115 * t86 + t88) + m(4) * t114 + (m(6) + m(7)) * (t48 ^ 2 + t51 ^ 2 + t88); -pkin(3) * t53 + t16 * t13 + t77 * t14 + t7 * t20 + t60 * t21 + (t84 / 0.2e1 + t82 / 0.2e1 + t56 / 0.2e1 - t55 / 0.2e1 + t57 / 0.2e1 + t54 / 0.2e1 - Ifges(4,6) - (t96 * mrSges(4,2))) * t92 - (t24 / 0.2e1 + t25 / 0.2e1) * t52 + (t22 / 0.2e1 - t23 / 0.2e1) * t50 + t118 * t31 + t117 * t29 + (t45 / 0.2e1 + pkin(8) * t64 + t37 * mrSges(5,3)) * t94 + (t46 / 0.2e1 - pkin(8) * t65 - t36 * mrSges(5,3)) * t91 + (Ifges(4,5) + t94 * t69 / 0.2e1 - t91 * t68 / 0.2e1 + t116 * t96) * t95 + m(5) * (pkin(3) * t121 + pkin(8) * t106) + m(6) * (-t29 * t5 + t31 * t6 + t60 * t77) + m(7) * (t16 * t7 + t2 * t31 + t29 * t3) + (t10 / 0.2e1 + t11 / 0.2e1 + t3 * mrSges(7,2) - t5 * mrSges(6,3)) * t62 + (-t9 / 0.2e1 + t8 / 0.2e1 - t2 * mrSges(7,2) - t6 * mrSges(6,3)) * t61; (-mrSges(4,2) + t111) * t92 + (-t20 - t21 + t116) * t95 + m(6) * (-t77 * t95 + t109) + m(7) * (-t16 * t95 + t109) + m(5) * (pkin(8) * t115 * t92 + pkin(3) * t95) + t131 * (t48 * t62 + t51 * t61); -0.2e1 * pkin(3) * t67 + 0.2e1 * t16 * t20 + 0.2e1 * t77 * t21 + t94 * t68 + t91 * t69 + Ifges(4,3) + m(6) * (t77 ^ 2 + t113) + m(7) * (t16 ^ 2 + t113) + m(5) * (pkin(8) ^ 2 * t115 + pkin(3) ^ 2) + (t24 + t25) * t62 + (t22 - t23) * t61 + 0.2e1 * pkin(8) * t111 + 0.2e1 * (t29 * t62 - t31 * t61) * t131; t100 + m(7) * (t2 * t74 + t3 * t76) + (t90 * t33 + t93 * t34 + m(6) * (t5 * t93 + t6 * t90)) * pkin(4) + t74 * t32 + t76 * t35 - Ifges(5,6) * t124 + t36 * mrSges(5,1) - t37 * mrSges(5,2) + t130; -t107 * t92 + m(7) * (t48 * t76 - t51 * t74) + m(6) * (-t48 * t93 - t51 * t90) * pkin(4) + t102; m(7) * (t29 * t76 + t31 * t74) + t84 + t82 - t107 * pkin(8) + (-t61 * t74 + t62 * t76) * mrSges(7,2) + (m(6) * (-t29 * t93 + t31 * t90) + (-t61 * t90 - t62 * t93) * mrSges(6,3)) * pkin(4) + t101; -0.2e1 * t76 * mrSges(7,1) + t74 * t129 + Ifges(5,3) + 0.2e1 * t103 + m(7) * (t74 ^ 2 + t76 ^ 2) + m(6) * (t90 ^ 2 + t93 ^ 2) * pkin(4) ^ 2 + t119; m(7) * (-pkin(5) * t3 + qJ(6) * t2) - pkin(5) * t35 + qJ(6) * t32 + t100; m(7) * (-pkin(5) * t48 - qJ(6) * t51) + t102; m(7) * (-pkin(5) * t29 + qJ(6) * t31) + (-pkin(5) * t62 - qJ(6) * t61) * mrSges(7,2) + t101; m(7) * (-pkin(5) * t76 + qJ(6) * t74) + t103 + (t74 + qJ(6)) * mrSges(7,3) + (-t76 + pkin(5)) * mrSges(7,1) + t119; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t129 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t119; m(7) * t3 + t35; m(7) * t48; m(7) * t29 + t62 * mrSges(7,2); m(7) * t76 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
