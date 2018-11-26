% Calculate joint inertia matrix for
% S6RPRRRP9
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
% Datum: 2018-11-23 16:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:28:54
% EndTime: 2018-11-23 16:28:55
% DurationCPUTime: 0.97s
% Computational Cost: add. (1216->279), mult. (2349->371), div. (0->0), fcn. (2210->6), ass. (0->106)
t119 = Ifges(6,3) + Ifges(7,3);
t99 = (-pkin(1) - pkin(7));
t137 = -2 * t99;
t97 = cos(qJ(4));
t98 = cos(qJ(3));
t122 = t97 * t98;
t95 = sin(qJ(3));
t136 = Ifges(5,5) * t122 + Ifges(5,3) * t95;
t120 = -mrSges(6,2) - mrSges(7,2);
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t135 = (t96 * mrSges(6,1) + t120 * t93) * pkin(4);
t134 = 2 * mrSges(7,1);
t133 = 2 * qJ(2);
t132 = m(7) * pkin(5);
t131 = -pkin(9) - pkin(8);
t130 = m(7) * t93;
t129 = pkin(4) * t96;
t94 = sin(qJ(4));
t128 = Ifges(5,4) * t94;
t127 = Ifges(5,4) * t97;
t62 = -t93 * t94 + t96 * t97;
t126 = t62 * t93;
t125 = t94 * t98;
t124 = t95 * t99;
t121 = t98 * t99;
t69 = pkin(3) * t95 - pkin(8) * t98 + qJ(2);
t60 = t97 * t69;
t19 = -pkin(9) * t122 + t60 + (-t94 * t99 + pkin(4)) * t95;
t37 = t97 * t124 + t94 * t69;
t22 = -pkin(9) * t125 + t37;
t6 = t93 * t19 + t96 * t22;
t63 = t93 * t97 + t94 * t96;
t50 = t63 * t98;
t32 = -mrSges(7,2) * t95 - mrSges(7,3) * t50;
t33 = -mrSges(6,2) * t95 - mrSges(6,3) * t50;
t118 = t32 + t33;
t73 = t131 * t94;
t74 = t131 * t97;
t31 = t93 * t73 - t96 * t74;
t70 = -mrSges(5,1) * t97 + mrSges(5,2) * t94;
t117 = -t70 + mrSges(4,1);
t116 = t94 ^ 2 + t97 ^ 2;
t89 = t95 ^ 2;
t91 = t98 ^ 2;
t115 = t91 + t89;
t5 = t96 * t19 - t22 * t93;
t52 = t62 * t98;
t2 = pkin(5) * t95 - qJ(6) * t52 + t5;
t34 = mrSges(7,1) * t95 - mrSges(7,3) * t52;
t114 = m(7) * t2 + t34;
t79 = -pkin(4) * t97 - pkin(3);
t112 = t116 * mrSges(5,3);
t111 = t115 * mrSges(4,3);
t15 = t50 * mrSges(7,1) + t52 * mrSges(7,2);
t23 = -t62 * mrSges(7,1) + t63 * mrSges(7,2);
t30 = t96 * t73 + t74 * t93;
t61 = pkin(4) * t125 - t121;
t13 = -qJ(6) * t63 + t30;
t110 = m(7) * t13 - t63 * mrSges(7,3);
t109 = mrSges(5,1) * t94 + mrSges(5,2) * t97;
t36 = -t124 * t94 + t60;
t108 = -t36 * t94 + t37 * t97;
t106 = t119 * t95 + (Ifges(6,5) + Ifges(7,5)) * t52 + (-Ifges(6,6) - Ifges(7,6)) * t50;
t49 = t63 * t95;
t51 = t62 * t95;
t105 = t120 * t51 + (-mrSges(6,1) - mrSges(7,1)) * t49;
t14 = qJ(6) * t62 + t31;
t55 = Ifges(7,6) * t62;
t56 = Ifges(6,6) * t62;
t57 = Ifges(7,5) * t63;
t58 = Ifges(6,5) * t63;
t104 = t30 * mrSges(6,1) + t13 * mrSges(7,1) - t31 * mrSges(6,2) - t14 * mrSges(7,2) + t55 + t56 + t57 + t58;
t3 = -qJ(6) * t50 + t6;
t103 = t5 * mrSges(6,1) + t2 * mrSges(7,1) - t6 * mrSges(6,2) - t3 * mrSges(7,2) + t106;
t102 = pkin(4) ^ 2;
t100 = qJ(2) ^ 2;
t92 = t99 ^ 2;
t87 = t93 ^ 2 * t102;
t86 = Ifges(5,5) * t94;
t85 = Ifges(5,6) * t97;
t81 = t91 * t99;
t80 = t91 * t92;
t78 = pkin(5) + t129;
t72 = Ifges(5,1) * t94 + t127;
t71 = Ifges(5,2) * t97 + t128;
t68 = mrSges(5,1) * t95 - mrSges(5,3) * t122;
t67 = -mrSges(5,2) * t95 - mrSges(5,3) * t125;
t53 = t109 * t98;
t48 = Ifges(5,5) * t95 + (Ifges(5,1) * t97 - t128) * t98;
t47 = Ifges(5,6) * t95 + (-Ifges(5,2) * t94 + t127) * t98;
t39 = -pkin(5) * t62 + t79;
t38 = t93 * pkin(4) * t51;
t35 = mrSges(6,1) * t95 - mrSges(6,3) * t52;
t28 = Ifges(6,1) * t63 + Ifges(6,4) * t62;
t27 = Ifges(7,1) * t63 + Ifges(7,4) * t62;
t26 = Ifges(6,4) * t63 + Ifges(6,2) * t62;
t25 = Ifges(7,4) * t63 + Ifges(7,2) * t62;
t24 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t21 = pkin(5) * t50 + t61;
t16 = mrSges(6,1) * t50 + mrSges(6,2) * t52;
t10 = Ifges(6,1) * t52 - Ifges(6,4) * t50 + Ifges(6,5) * t95;
t9 = Ifges(7,1) * t52 - Ifges(7,4) * t50 + Ifges(7,5) * t95;
t8 = Ifges(6,4) * t52 - Ifges(6,2) * t50 + Ifges(6,6) * t95;
t7 = Ifges(7,4) * t52 - Ifges(7,2) * t50 + Ifges(7,6) * t95;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t133) + 0.2e1 * t21 * t15 + 0.2e1 * t61 * t16 + 0.2e1 * t2 * t34 + 0.2e1 * t3 * t32 + 0.2e1 * t6 * t33 + 0.2e1 * t5 * t35 + 0.2e1 * t36 * t68 + 0.2e1 * t37 * t67 + Ifges(3,1) + Ifges(2,3) + (t9 + t10) * t52 - (t7 + t8) * t50 + t111 * t137 + ((mrSges(4,2) * t133) + Ifges(4,1) * t98 + t53 * t137 - t47 * t94 + t48 * t97) * t98 + (mrSges(4,1) * t133 + Ifges(4,2) * t95 + (-Ifges(5,6) * t94 - (2 * Ifges(4,4))) * t98 + t106 + t136) * t95 + m(4) * (t89 * t92 + t100 + t80) + (m(3) * (pkin(1) ^ 2 + t100)) + m(5) * (t36 ^ 2 + t37 ^ 2 + t80) + m(6) * (t5 ^ 2 + t6 ^ 2 + t61 ^ 2) + m(7) * (t2 ^ 2 + t21 ^ 2 + t3 ^ 2); -(m(3) * pkin(1)) + mrSges(3,2) + (t97 * t67 - t94 * t68) * t95 + t118 * t51 - (t34 + t35) * t49 - t111 + (-t15 - t16 - t53) * t98 + m(7) * (-t2 * t49 - t21 * t98 + t3 * t51) + m(6) * (-t49 * t5 + t51 * t6 - t61 * t98) + m(5) * (t108 * t95 + t81) + m(4) * (t89 * t99 + t81); m(3) + m(5) * (t116 * t89 + t91) + m(4) * t115 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t49 ^ 2 + t51 ^ 2 + t91); -pkin(3) * t53 + t13 * t34 + t14 * t32 + t39 * t15 + t79 * t16 + t21 * t23 + t61 * t24 + t30 * t35 + t31 * t33 + (-(t99 * mrSges(4,2)) + t86 / 0.2e1 + t85 / 0.2e1 + t57 / 0.2e1 + t55 / 0.2e1 + t58 / 0.2e1 + t56 / 0.2e1 - Ifges(4,6)) * t95 + (t27 / 0.2e1 + t28 / 0.2e1) * t52 - (t25 / 0.2e1 + t26 / 0.2e1) * t50 + (pkin(8) * t67 + t37 * mrSges(5,3) + t47 / 0.2e1) * t97 + (-pkin(8) * t68 - t36 * mrSges(5,3) + t48 / 0.2e1) * t94 + (t97 * t72 / 0.2e1 - t94 * t71 / 0.2e1 + Ifges(4,5) + t117 * t99) * t98 + m(5) * (pkin(3) * t121 + pkin(8) * t108) + m(6) * (t30 * t5 + t31 * t6 + t79 * t61) + m(7) * (t13 * t2 + t14 * t3 + t21 * t39) + (-t5 * mrSges(6,3) - t2 * mrSges(7,3) + t9 / 0.2e1 + t10 / 0.2e1) * t63 + (t3 * mrSges(7,3) + t6 * mrSges(6,3) + t7 / 0.2e1 + t8 / 0.2e1) * t62; (-mrSges(4,2) + t112) * t95 + (-t23 - t24 + t117) * t98 + m(6) * (-t30 * t49 + t31 * t51 - t79 * t98) + m(7) * (-t13 * t49 + t14 * t51 - t39 * t98) + m(5) * (pkin(8) * t116 * t95 + pkin(3) * t98) + (mrSges(7,3) + mrSges(6,3)) * (t49 * t63 + t51 * t62); -0.2e1 * pkin(3) * t70 + 0.2e1 * t39 * t23 + 0.2e1 * t79 * t24 + t97 * t71 + t94 * t72 + Ifges(4,3) + 0.2e1 * pkin(8) * t112 + m(6) * (t30 ^ 2 + t31 ^ 2 + t79 ^ 2) + m(7) * (t13 ^ 2 + t14 ^ 2 + t39 ^ 2) + m(5) * (pkin(8) ^ 2 * t116 + pkin(3) ^ 2) + (-0.2e1 * mrSges(6,3) * t30 - 0.2e1 * mrSges(7,3) * t13 + t27 + t28) * t63 + (0.2e1 * mrSges(6,3) * t31 + 0.2e1 * mrSges(7,3) * t14 + t25 + t26) * t62; t103 + t114 * t78 + (t96 * t35 + t118 * t93 + t3 * t130 + m(6) * (t5 * t96 + t6 * t93)) * pkin(4) - Ifges(5,6) * t125 + t36 * mrSges(5,1) - t37 * mrSges(5,2) + t136; -t109 * t95 + m(6) * (-t129 * t49 + t38) + m(7) * (-t49 * t78 + t38) + t105; t85 + t86 + t110 * t78 - t109 * pkin(8) + (mrSges(7,3) * t126 + (-t63 * t96 + t126) * mrSges(6,3) + m(6) * (t30 * t96 + t31 * t93) + t14 * t130) * pkin(4) + t104; t78 * t134 + Ifges(5,3) + m(6) * (t102 * t96 ^ 2 + t87) + m(7) * (t78 ^ 2 + t87) + 0.2e1 * t135 + t119; pkin(5) * t114 + t103; -t132 * t49 + t105; pkin(5) * t110 + t104; t78 * t132 + (pkin(5) + t78) * mrSges(7,1) + t135 + t119; (t134 + t132) * pkin(5) + t119; m(7) * t21 + t15; -m(7) * t98; m(7) * t39 + t23; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
