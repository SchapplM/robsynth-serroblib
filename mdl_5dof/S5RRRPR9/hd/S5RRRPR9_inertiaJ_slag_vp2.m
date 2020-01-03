% Calculate joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:16
% EndTime: 2019-12-31 21:22:19
% DurationCPUTime: 0.79s
% Computational Cost: add. (1259->231), mult. (2573->344), div. (0->0), fcn. (2619->8), ass. (0->97)
t120 = 2 * pkin(6);
t96 = sin(qJ(2));
t98 = cos(qJ(3));
t111 = t96 * t98;
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t95 = sin(qJ(3));
t64 = t92 * t98 + t93 * t95;
t55 = t64 * t96;
t63 = -t92 * t95 + t93 * t98;
t56 = t63 * t96;
t119 = -Ifges(4,5) * t111 - Ifges(5,5) * t56 + Ifges(5,6) * t55;
t118 = pkin(3) * t92;
t99 = cos(qJ(2));
t117 = pkin(6) * t99;
t87 = t96 * pkin(6);
t116 = Ifges(4,4) * t95;
t115 = Ifges(4,4) * t98;
t82 = t93 * pkin(3) + pkin(4);
t94 = sin(qJ(5));
t97 = cos(qJ(5));
t57 = -t94 * t118 + t97 * t82;
t114 = t57 * mrSges(6,1);
t58 = t97 * t118 + t94 * t82;
t113 = t58 * mrSges(6,2);
t112 = t95 * t96;
t110 = -qJ(4) - pkin(7);
t24 = -t97 * t55 - t94 * t56;
t25 = -t94 * t55 + t97 * t56;
t109 = -Ifges(6,5) * t25 - Ifges(6,6) * t24;
t107 = qJ(4) * t96;
t73 = -t99 * pkin(2) - t96 * pkin(7) - pkin(1);
t66 = t98 * t73;
t36 = -t98 * t107 + t66 + (-pkin(6) * t95 - pkin(3)) * t99;
t48 = t98 * t117 + t95 * t73;
t42 = -t95 * t107 + t48;
t16 = t92 * t36 + t93 * t42;
t74 = t110 * t95;
t76 = t110 * t98;
t44 = t92 * t74 - t93 * t76;
t72 = pkin(3) * t112 + t87;
t108 = t95 ^ 2 + t98 ^ 2;
t106 = -Ifges(6,3) - Ifges(5,3) - Ifges(4,3);
t83 = -t98 * pkin(3) - pkin(2);
t28 = t55 * mrSges(5,1) + t56 * mrSges(5,2);
t39 = -t63 * mrSges(5,1) + t64 * mrSges(5,2);
t6 = -t24 * mrSges(6,1) + t25 * mrSges(6,2);
t34 = t97 * t63 - t94 * t64;
t35 = t94 * t63 + t97 * t64;
t12 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t15 = t93 * t36 - t92 * t42;
t43 = t93 * t74 + t92 * t76;
t105 = t95 * mrSges(4,1) + t98 * mrSges(4,2);
t11 = -t55 * pkin(8) + t16;
t8 = -t99 * pkin(4) - t56 * pkin(8) + t15;
t2 = -t94 * t11 + t97 * t8;
t3 = t97 * t11 + t94 * t8;
t104 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t109;
t26 = -t64 * pkin(8) + t43;
t27 = t63 * pkin(8) + t44;
t10 = t94 * t26 + t97 * t27;
t32 = Ifges(6,6) * t34;
t33 = Ifges(6,5) * t35;
t9 = t97 * t26 - t94 * t27;
t103 = t9 * mrSges(6,1) - t10 * mrSges(6,2) + t32 + t33;
t101 = pkin(6) ^ 2;
t91 = t99 ^ 2;
t89 = t96 ^ 2;
t86 = t89 * t101;
t85 = Ifges(4,5) * t95;
t84 = Ifges(4,6) * t98;
t78 = Ifges(4,1) * t95 + t115;
t77 = Ifges(4,2) * t98 + t116;
t75 = -t98 * mrSges(4,1) + t95 * mrSges(4,2);
t71 = -t99 * mrSges(4,1) - mrSges(4,3) * t111;
t70 = t99 * mrSges(4,2) - mrSges(4,3) * t112;
t62 = t105 * t96;
t61 = Ifges(5,5) * t64;
t60 = Ifges(5,6) * t63;
t54 = -Ifges(4,5) * t99 + (Ifges(4,1) * t98 - t116) * t96;
t53 = -Ifges(4,6) * t99 + (-Ifges(4,2) * t95 + t115) * t96;
t49 = -t63 * pkin(4) + t83;
t47 = -t95 * t117 + t66;
t46 = -t99 * mrSges(5,1) - t56 * mrSges(5,3);
t45 = t99 * mrSges(5,2) - t55 * mrSges(5,3);
t41 = Ifges(5,1) * t64 + Ifges(5,4) * t63;
t40 = Ifges(5,4) * t64 + Ifges(5,2) * t63;
t37 = t55 * pkin(4) + t72;
t23 = Ifges(5,1) * t56 - Ifges(5,4) * t55 - Ifges(5,5) * t99;
t22 = Ifges(5,4) * t56 - Ifges(5,2) * t55 - Ifges(5,6) * t99;
t18 = -t99 * mrSges(6,1) - t25 * mrSges(6,3);
t17 = t99 * mrSges(6,2) + t24 * mrSges(6,3);
t14 = Ifges(6,1) * t35 + Ifges(6,4) * t34;
t13 = Ifges(6,4) * t35 + Ifges(6,2) * t34;
t5 = Ifges(6,1) * t25 + Ifges(6,4) * t24 - Ifges(6,5) * t99;
t4 = Ifges(6,4) * t25 + Ifges(6,2) * t24 - Ifges(6,6) * t99;
t1 = [0.2e1 * t15 * t46 + 0.2e1 * t16 * t45 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t55 * t22 + t56 * t23 + t24 * t4 + t25 * t5 + 0.2e1 * t72 * t28 + 0.2e1 * t37 * t6 + 0.2e1 * t47 * t71 + 0.2e1 * t48 * t70 + Ifges(2,3) + (t89 + t91) * mrSges(3,3) * t120 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t96 + t62 * t120 - t95 * t53 + t98 * t54) * t96 + m(3) * (pkin(1) ^ 2 + t91 * t101 + t86) + m(5) * (t15 ^ 2 + t16 ^ 2 + t72 ^ 2) + m(4) * (t47 ^ 2 + t48 ^ 2 + t86) + m(6) * (t2 ^ 2 + t3 ^ 2 + t37 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) - t106) * t99 + (Ifges(4,6) * t95 + (2 * Ifges(3,4))) * t96 + t109 + t119) * t99; (Ifges(3,5) + t98 * t78 / 0.2e1 - t95 * t77 / 0.2e1 + (t75 - mrSges(3,1)) * pkin(6)) * t96 + (t53 / 0.2e1 + pkin(7) * t70 + t48 * mrSges(4,3)) * t98 + (t54 / 0.2e1 - pkin(7) * t71 - t47 * mrSges(4,3)) * t95 + m(5) * (t43 * t15 + t44 * t16 + t83 * t72) + m(6) * (t10 * t3 + t9 * t2 + t49 * t37) + (-t85 / 0.2e1 - t84 / 0.2e1 - t61 / 0.2e1 - t60 / 0.2e1 - t33 / 0.2e1 - t32 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t99 + (-t2 * t35 + t3 * t34) * mrSges(6,3) + (-t15 * t64 + t16 * t63) * mrSges(5,3) + m(4) * (-pkin(2) * t87 + (-t47 * t95 + t48 * t98) * pkin(7)) + t72 * t39 + t83 * t28 - pkin(2) * t62 + t63 * t22 / 0.2e1 + t64 * t23 / 0.2e1 + t43 * t46 + t49 * t6 - t55 * t40 / 0.2e1 + t56 * t41 / 0.2e1 + t34 * t4 / 0.2e1 + t35 * t5 / 0.2e1 + t37 * t12 + t44 * t45 + t24 * t13 / 0.2e1 + t25 * t14 / 0.2e1 + t10 * t17 + t9 * t18; -0.2e1 * pkin(2) * t75 + 0.2e1 * t49 * t12 + t34 * t13 + t35 * t14 + 0.2e1 * t83 * t39 + t63 * t40 + t64 * t41 + t98 * t77 + t95 * t78 + Ifges(3,3) + m(6) * (t10 ^ 2 + t49 ^ 2 + t9 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2 + t83 ^ 2) + m(4) * (t108 * pkin(7) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (t10 * t34 - t9 * t35) * mrSges(6,3) + 0.2e1 * (-t43 * t64 + t44 * t63) * mrSges(5,3) + 0.2e1 * t108 * pkin(7) * mrSges(4,3); -t48 * mrSges(4,2) + t47 * mrSges(4,1) - Ifges(4,6) * t112 + t58 * t17 + t57 * t18 + m(6) * (t57 * t2 + t58 * t3) + t15 * mrSges(5,1) - t16 * mrSges(5,2) + t106 * t99 + (t93 * t46 + t92 * t45 + m(5) * (t15 * t93 + t16 * t92)) * pkin(3) + t104 - t119; m(6) * (t58 * t10 + t57 * t9) - t44 * mrSges(5,2) + t43 * mrSges(5,1) + t61 + t60 + t85 + t84 - t105 * pkin(7) + (t58 * t34 - t57 * t35) * mrSges(6,3) + (m(5) * (t43 * t93 + t44 * t92) + (t92 * t63 - t93 * t64) * mrSges(5,3)) * pkin(3) + t103; 0.2e1 * t114 - 0.2e1 * t113 + m(6) * (t57 ^ 2 + t58 ^ 2) - t106 + (0.2e1 * t93 * mrSges(5,1) - 0.2e1 * t92 * mrSges(5,2) + m(5) * (t92 ^ 2 + t93 ^ 2) * pkin(3)) * pkin(3); m(5) * t72 + m(6) * t37 + t28 + t6; m(5) * t83 + m(6) * t49 + t12 + t39; 0; m(5) + m(6); -Ifges(6,3) * t99 + t104; t103; Ifges(6,3) - t113 + t114; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
