% Calculate joint inertia matrix for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:18
% EndTime: 2019-12-05 17:23:21
% DurationCPUTime: 0.82s
% Computational Cost: add. (1011->240), mult. (2549->365), div. (0->0), fcn. (2691->12), ass. (0->99)
t126 = 2 * pkin(9);
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t58 = -t92 * mrSges(6,1) + t88 * mrSges(6,2);
t125 = -m(6) * pkin(4) - mrSges(5,1) + t58;
t86 = cos(pkin(6));
t95 = cos(qJ(2));
t109 = t86 * t95;
t84 = sin(pkin(6));
t90 = sin(qJ(3));
t111 = t84 * t90;
t85 = sin(pkin(5));
t87 = cos(pkin(5));
t91 = sin(qJ(2));
t94 = cos(qJ(3));
t24 = t87 * t111 + (t90 * t109 + t91 * t94) * t85;
t45 = -t85 * t95 * t84 + t87 * t86;
t89 = sin(qJ(4));
t93 = cos(qJ(4));
t9 = t24 * t89 - t45 * t93;
t124 = t9 ^ 2;
t110 = t84 * t94;
t22 = -t87 * t110 + (-t109 * t94 + t90 * t91) * t85;
t123 = t22 ^ 2;
t47 = t93 * t111 + t89 * t86;
t27 = -t92 * t110 - t88 * t47;
t122 = t27 / 0.2e1;
t28 = -t88 * t110 + t92 * t47;
t121 = t28 / 0.2e1;
t120 = -t88 / 0.2e1;
t119 = t88 / 0.2e1;
t118 = t92 / 0.2e1;
t117 = pkin(2) * t94;
t116 = pkin(9) * t93;
t115 = t9 * t89;
t30 = -mrSges(5,1) * t110 - t47 * mrSges(5,3);
t8 = -t27 * mrSges(6,1) + t28 * mrSges(6,2);
t114 = -t30 + t8;
t113 = Ifges(6,4) * t88;
t112 = Ifges(6,4) * t92;
t108 = t88 * t89;
t107 = t89 * t92;
t50 = t86 * t90 * pkin(2) + pkin(8) * t110;
t37 = t86 * pkin(9) + t50;
t38 = (-pkin(3) * t94 - pkin(9) * t90 - pkin(2)) * t84;
t18 = t93 * t37 + t89 * t38;
t46 = t89 * t111 - t93 * t86;
t106 = -Ifges(5,5) * t47 + Ifges(5,6) * t46;
t60 = Ifges(6,5) * t88 + Ifges(6,6) * t92;
t105 = Ifges(5,5) * t89 + Ifges(5,6) * t93;
t104 = t88 ^ 2 + t92 ^ 2;
t5 = Ifges(6,5) * t28 + Ifges(6,6) * t27 + Ifges(6,3) * t46;
t103 = Ifges(4,5) * t111 + Ifges(4,6) * t110 + Ifges(4,3) * t86;
t68 = pkin(8) * t111;
t36 = t68 + (-pkin(3) - t117) * t86;
t12 = t46 * pkin(4) - t47 * pkin(10) + t36;
t14 = -pkin(10) * t110 + t18;
t1 = t92 * t12 - t88 * t14;
t2 = t88 * t12 + t92 * t14;
t102 = -t1 * t88 + t2 * t92;
t11 = t24 * t93 + t45 * t89;
t100 = t11 * t93 + t115;
t99 = mrSges(6,1) * t88 + mrSges(6,2) * t92;
t57 = -t93 * pkin(4) - t89 * pkin(10) - pkin(3);
t33 = -t88 * t116 + t92 * t57;
t34 = t92 * t116 + t88 * t57;
t98 = -t33 * t88 + t34 * t92;
t17 = -t89 * t37 + t93 * t38;
t39 = Ifges(6,5) * t107 - Ifges(6,6) * t108 - Ifges(6,3) * t93;
t97 = pkin(9) ^ 2;
t83 = t93 ^ 2;
t81 = t89 ^ 2;
t78 = t81 * t97;
t64 = Ifges(5,1) * t89 + Ifges(5,4) * t93;
t63 = Ifges(6,1) * t88 + t112;
t62 = Ifges(5,4) * t89 + Ifges(5,2) * t93;
t61 = Ifges(6,2) * t92 + t113;
t59 = -t93 * mrSges(5,1) + t89 * mrSges(5,2);
t55 = -t93 * mrSges(6,1) - mrSges(6,3) * t107;
t54 = t93 * mrSges(6,2) - mrSges(6,3) * t108;
t53 = -t86 * mrSges(4,2) + mrSges(4,3) * t110;
t52 = t86 * mrSges(4,1) - mrSges(4,3) * t111;
t51 = t99 * t89;
t49 = t86 * t117 - t68;
t48 = (-mrSges(4,1) * t94 + mrSges(4,2) * t90) * t84;
t41 = -Ifges(6,5) * t93 + (Ifges(6,1) * t92 - t113) * t89;
t40 = -Ifges(6,6) * t93 + (-Ifges(6,2) * t88 + t112) * t89;
t29 = mrSges(5,2) * t110 - t46 * mrSges(5,3);
t21 = t46 * mrSges(5,1) + t47 * mrSges(5,2);
t20 = Ifges(5,1) * t47 - Ifges(5,4) * t46 - Ifges(5,5) * t110;
t19 = Ifges(5,4) * t47 - Ifges(5,2) * t46 - Ifges(5,6) * t110;
t16 = t46 * mrSges(6,1) - t28 * mrSges(6,3);
t15 = -t46 * mrSges(6,2) + t27 * mrSges(6,3);
t13 = pkin(4) * t110 - t17;
t7 = Ifges(6,1) * t28 + Ifges(6,4) * t27 + Ifges(6,5) * t46;
t6 = Ifges(6,4) * t28 + Ifges(6,2) * t27 + Ifges(6,6) * t46;
t4 = t11 * t92 + t22 * t88;
t3 = -t11 * t88 + t22 * t92;
t10 = [m(2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t124) + m(5) * (t11 ^ 2 + t123 + t124) + m(4) * (t24 ^ 2 + t45 ^ 2 + t123) + m(3) * (t87 ^ 2 + (t91 ^ 2 + t95 ^ 2) * t85 ^ 2); t11 * t29 + t4 * t15 + t3 * t16 + t24 * t53 + t45 * t48 + t114 * t9 + (t95 * mrSges(3,1) - t91 * mrSges(3,2)) * t85 + (t21 - t52) * t22 + m(6) * (t1 * t3 + t13 * t9 + t2 * t4) + m(5) * (t18 * t11 - t17 * t9 + t36 * t22) + m(4) * (-t84 * pkin(2) * t45 - t49 * t22 + t50 * t24); Ifges(3,3) + 0.2e1 * t1 * t16 + 0.2e1 * t2 * t15 + t28 * t7 + t27 * t6 + 0.2e1 * t13 * t8 + 0.2e1 * t18 * t29 + 0.2e1 * t17 * t30 + 0.2e1 * t36 * t21 + t47 * t20 + 0.2e1 * t50 * t53 + t86 * t103 + 0.2e1 * t49 * t52 + (t5 - t19) * t46 + (-0.2e1 * pkin(2) * t48 + (Ifges(4,1) * t111 + Ifges(4,5) * t86) * t90 + (Ifges(4,6) * t86 + (0.2e1 * Ifges(4,4) * t90 + (Ifges(5,3) + Ifges(4,2)) * t94) * t84 + t106) * t94) * t84 + m(6) * (t1 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t36 ^ 2) + m(4) * (t84 ^ 2 * pkin(2) ^ 2 + t49 ^ 2 + t50 ^ 2); -t24 * mrSges(4,2) + t3 * t55 + t4 * t54 + t9 * t51 + (-mrSges(4,1) + t59) * t22 + t100 * mrSges(5,3) + m(6) * (pkin(9) * t115 + t33 * t3 + t34 * t4) + m(5) * (-pkin(3) * t22 + t100 * pkin(9)); m(5) * (-pkin(3) * t36 + (-t17 * t89 + t18 * t93) * pkin(9)) + (t7 * t118 + t6 * t120 - t17 * mrSges(5,3) + t20 / 0.2e1 + t114 * pkin(9)) * t89 + (t18 * mrSges(5,3) + pkin(9) * t29 + t19 / 0.2e1 - t5 / 0.2e1) * t93 + t103 + m(6) * (t89 * pkin(9) * t13 + t33 * t1 + t34 * t2) + (-t62 / 0.2e1 + t39 / 0.2e1) * t46 - t105 * t110 / 0.2e1 + t49 * mrSges(4,1) - t50 * mrSges(4,2) + t13 * t51 + t2 * t54 + t1 * t55 + t36 * t59 + t47 * t64 / 0.2e1 + t40 * t122 + t41 * t121 + t33 * t16 + t34 * t15 - pkin(3) * t21; -0.2e1 * pkin(3) * t59 + 0.2e1 * t33 * t55 + 0.2e1 * t34 * t54 + Ifges(4,3) + (-t39 + t62) * t93 + (t81 + t83) * mrSges(5,3) * t126 + m(6) * (t33 ^ 2 + t34 ^ 2 + t78) + m(5) * (pkin(3) ^ 2 + t83 * t97 + t78) + (t51 * t126 - t88 * t40 + t92 * t41 + t64) * t89; -t11 * mrSges(5,2) + (m(6) * pkin(10) + mrSges(6,3)) * (-t3 * t88 + t4 * t92) + t125 * t9; -Ifges(5,3) * t110 + t63 * t121 + t46 * t60 / 0.2e1 + t61 * t122 + t6 * t118 + t13 * t58 + t7 * t119 + t17 * mrSges(5,1) - t18 * mrSges(5,2) + t102 * mrSges(6,3) - t106 + (-m(6) * t13 - t8) * pkin(4) + (m(6) * t102 + t92 * t15 - t88 * t16) * pkin(10); t41 * t119 + t40 * t118 - pkin(4) * t51 + (m(6) * t98 + t92 * t54 - t88 * t55) * pkin(10) + (t125 * pkin(9) + t63 * t118 + t61 * t120) * t89 + (-pkin(9) * mrSges(5,2) - t60 / 0.2e1) * t93 + t98 * mrSges(6,3) + t105; Ifges(5,3) + t88 * t63 + t92 * t61 - 0.2e1 * pkin(4) * t58 + m(6) * (t104 * pkin(10) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t104 * pkin(10) * mrSges(6,3); t3 * mrSges(6,1) - t4 * mrSges(6,2); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t5; t33 * mrSges(6,1) - t34 * mrSges(6,2) + t39; -t99 * pkin(10) + t60; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
