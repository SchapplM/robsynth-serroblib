% Calculate joint inertia matrix for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:52:16
% EndTime: 2018-11-23 14:52:17
% DurationCPUTime: 0.83s
% Computational Cost: add. (1273->265), mult. (3244->397), div. (0->0), fcn. (3657->14), ass. (0->111)
t130 = 2 * pkin(9);
t85 = cos(pkin(13));
t86 = cos(pkin(7));
t118 = t85 * t86;
t83 = sin(pkin(7));
t91 = sin(qJ(3));
t120 = t83 * t91;
t82 = sin(pkin(13));
t84 = sin(pkin(6));
t87 = cos(pkin(6));
t95 = cos(qJ(3));
t23 = t87 * t120 + (t91 * t118 + t82 * t95) * t84;
t44 = -t83 * t84 * t85 + t86 * t87;
t90 = sin(qJ(4));
t94 = cos(qJ(4));
t14 = t23 * t90 - t94 * t44;
t13 = t14 ^ 2;
t119 = t83 * t95;
t21 = -t87 * t119 + (-t118 * t95 + t82 * t91) * t84;
t129 = t21 ^ 2;
t46 = t90 * t120 - t94 * t86;
t45 = t46 ^ 2;
t128 = m(7) * pkin(5);
t93 = cos(qJ(5));
t127 = t93 / 0.2e1;
t126 = -pkin(11) - pkin(10);
t125 = pkin(9) * t94;
t89 = sin(qJ(5));
t124 = Ifges(6,4) * t89;
t123 = Ifges(6,4) * t93;
t122 = t14 * t90;
t6 = t46 * t14;
t121 = t46 * t90;
t117 = t89 * t90;
t116 = t90 * t93;
t115 = t94 * mrSges(5,3);
t62 = -mrSges(5,1) * t94 + mrSges(5,2) * t90;
t114 = mrSges(4,1) - t62;
t113 = -Ifges(7,3) - Ifges(6,3);
t88 = sin(qJ(6));
t92 = cos(qJ(6));
t55 = t88 * t93 + t89 * t92;
t42 = t55 * t90;
t54 = -t88 * t89 + t92 * t93;
t43 = t54 * t90;
t112 = -Ifges(7,5) * t43 + Ifges(7,6) * t42;
t61 = -mrSges(6,1) * t93 + mrSges(6,2) * t89;
t111 = t61 - mrSges(5,1);
t60 = -pkin(4) * t94 - pkin(10) * t90 - pkin(3);
t36 = t93 * t125 + t89 * t60;
t110 = t89 ^ 2 + t93 ^ 2;
t24 = -mrSges(7,1) * t54 + mrSges(7,2) * t55;
t109 = t24 + t111;
t16 = t23 * t94 + t44 * t90;
t4 = -t16 * t89 + t21 * t93;
t5 = t16 * t93 + t21 * t89;
t2 = t4 * t92 - t5 * t88;
t3 = t4 * t88 + t5 * t92;
t108 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t48 = t94 * t120 + t86 * t90;
t27 = -t93 * t119 - t48 * t89;
t28 = -t89 * t119 + t48 * t93;
t11 = t27 * t92 - t28 * t88;
t12 = t27 * t88 + t28 * t92;
t107 = t11 * mrSges(7,1) - t12 * mrSges(7,2);
t19 = mrSges(7,1) * t42 + mrSges(7,2) * t43;
t104 = mrSges(6,1) * t89 + mrSges(6,2) * t93;
t49 = t104 * t90;
t106 = t90 * mrSges(5,3) + t19 + t49;
t105 = -t4 * t89 + t5 * t93;
t103 = -t27 * t89 + t28 * t93;
t53 = t93 * t60;
t35 = -t89 * t125 + t53;
t102 = -t35 * t89 + t36 * t93;
t20 = -pkin(11) * t116 + t53 + (-pkin(9) * t89 - pkin(5)) * t94;
t29 = -pkin(11) * t117 + t36;
t8 = t20 * t92 - t29 * t88;
t9 = t20 * t88 + t29 * t92;
t101 = t8 * mrSges(7,1) - t9 * mrSges(7,2) - t112;
t66 = t126 * t89;
t67 = t126 * t93;
t31 = t66 * t92 + t67 * t88;
t32 = t66 * t88 - t67 * t92;
t50 = Ifges(7,6) * t54;
t51 = Ifges(7,5) * t55;
t100 = t31 * mrSges(7,1) - t32 * mrSges(7,2) + t50 + t51;
t99 = (mrSges(7,1) * t92 - mrSges(7,2) * t88) * pkin(5);
t97 = pkin(9) ^ 2;
t81 = t94 ^ 2;
t79 = t90 ^ 2;
t76 = t83 ^ 2;
t75 = t79 * t97;
t74 = Ifges(6,5) * t89;
t73 = Ifges(6,6) * t93;
t72 = -pkin(5) * t93 - pkin(4);
t70 = t76 * t95 ^ 2;
t68 = Ifges(6,5) * t116;
t64 = Ifges(6,1) * t89 + t123;
t63 = Ifges(6,2) * t93 + t124;
t59 = (pkin(5) * t89 + pkin(9)) * t90;
t58 = -mrSges(6,1) * t94 - mrSges(6,3) * t116;
t57 = mrSges(6,2) * t94 - mrSges(6,3) * t117;
t41 = -Ifges(6,5) * t94 + (Ifges(6,1) * t93 - t124) * t90;
t40 = -Ifges(6,6) * t94 + (-Ifges(6,2) * t89 + t123) * t90;
t34 = -mrSges(7,1) * t94 - mrSges(7,3) * t43;
t33 = mrSges(7,2) * t94 - mrSges(7,3) * t42;
t26 = Ifges(7,1) * t55 + Ifges(7,4) * t54;
t25 = Ifges(7,4) * t55 + Ifges(7,2) * t54;
t18 = Ifges(7,1) * t43 - Ifges(7,4) * t42 - Ifges(7,5) * t94;
t17 = Ifges(7,4) * t43 - Ifges(7,2) * t42 - Ifges(7,6) * t94;
t1 = [m(2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t13) + m(6) * (t4 ^ 2 + t5 ^ 2 + t13) + m(5) * (t16 ^ 2 + t129 + t13) + m(4) * (t23 ^ 2 + t44 ^ 2 + t129) + m(3) * (t87 ^ 2 + (t82 ^ 2 + t85 ^ 2) * t84 ^ 2); m(3) * t87 + m(7) * (t11 * t2 + t12 * t3 + t6) + m(6) * (t27 * t4 + t28 * t5 + t6) + m(5) * (-t21 * t119 + t16 * t48 + t6) + m(4) * (t44 * t86 + (-t21 * t95 + t23 * t91) * t83); m(3) + m(7) * (t11 ^ 2 + t12 ^ 2 + t45) + m(6) * (t27 ^ 2 + t28 ^ 2 + t45) + m(5) * (t48 ^ 2 + t45 + t70) + m(4) * (t76 * t91 ^ 2 + t86 ^ 2 + t70); t16 * t115 - t23 * mrSges(4,2) + t2 * t34 + t3 * t33 + t4 * t58 + t5 * t57 - t114 * t21 + t106 * t14 + m(7) * (t14 * t59 + t2 * t8 + t3 * t9) + m(6) * (pkin(9) * t122 + t35 * t4 + t36 * t5) + m(5) * (-pkin(3) * t21 + (t16 * t94 + t122) * pkin(9)); t48 * t115 + t11 * t34 + t12 * t33 + t27 * t58 + t28 * t57 + (-t91 * mrSges(4,2) + t114 * t95) * t83 + t106 * t46 + m(7) * (t11 * t8 + t12 * t9 + t46 * t59) + m(6) * (pkin(9) * t121 + t27 * t35 + t28 * t36) + m(5) * (pkin(3) * t119 + (t48 * t94 + t121) * pkin(9)); -0.2e1 * pkin(3) * t62 - t42 * t17 + t43 * t18 + 0.2e1 * t59 * t19 + 0.2e1 * t9 * t33 + 0.2e1 * t8 * t34 + 0.2e1 * t35 * t58 + 0.2e1 * t36 * t57 + Ifges(4,3) + (t79 + t81) * mrSges(5,3) * t130 + m(7) * (t59 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2 + t75) + m(5) * (pkin(3) ^ 2 + t81 * t97 + t75) + (Ifges(5,1) * t90 + t49 * t130 - t40 * t89 + t41 * t93) * t90 + (-t68 + (Ifges(5,2) - t113) * t94 + (Ifges(6,6) * t89 + (2 * Ifges(5,4))) * t90 + t112) * t94; -t16 * mrSges(5,2) + (-t2 * t55 + t3 * t54) * mrSges(7,3) + t105 * mrSges(6,3) + t109 * t14 + m(7) * (t14 * t72 + t2 * t31 + t3 * t32) + m(6) * (-pkin(4) * t14 + t105 * pkin(10)); -t48 * mrSges(5,2) + (-t11 * t55 + t12 * t54) * mrSges(7,3) + t103 * mrSges(6,3) + t109 * t46 + m(7) * (t11 * t31 + t12 * t32 + t46 * t72) + m(6) * (-pkin(4) * t46 + t103 * pkin(10)); t40 * t127 + t89 * t41 / 0.2e1 + m(7) * (t31 * t8 + t32 * t9 + t72 * t59) + t72 * t19 + t54 * t17 / 0.2e1 + t55 * t18 / 0.2e1 + t59 * t24 - t42 * t25 / 0.2e1 + t43 * t26 / 0.2e1 - pkin(4) * t49 + t32 * t33 + t31 * t34 + (-t74 / 0.2e1 - t73 / 0.2e1 - t51 / 0.2e1 - t50 / 0.2e1 + Ifges(5,6) - pkin(9) * mrSges(5,2)) * t94 + (t54 * t9 - t55 * t8) * mrSges(7,3) + t102 * mrSges(6,3) + (m(6) * t102 + t93 * t57 - t89 * t58) * pkin(10) + (Ifges(5,5) - t89 * t63 / 0.2e1 + t64 * t127 + (-m(6) * pkin(4) + t111) * pkin(9)) * t90; -0.2e1 * pkin(4) * t61 + 0.2e1 * t72 * t24 + t54 * t25 + t55 * t26 + t93 * t63 + t89 * t64 + Ifges(5,3) + m(7) * (t31 ^ 2 + t32 ^ 2 + t72 ^ 2) + m(6) * (t110 * pkin(10) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t31 * t55 + t32 * t54) * mrSges(7,3) + 0.2e1 * t110 * pkin(10) * mrSges(6,3); t4 * mrSges(6,1) - t5 * mrSges(6,2) + (t2 * t92 + t3 * t88) * t128 + t108; t27 * mrSges(6,1) - t28 * mrSges(6,2) + (t11 * t92 + t12 * t88) * t128 + t107; -Ifges(6,6) * t117 + t35 * mrSges(6,1) - t36 * mrSges(6,2) + t68 + t113 * t94 + (m(7) * (t8 * t92 + t88 * t9) + t92 * t34 + t88 * t33) * pkin(5) + t101; t73 + t74 - t104 * pkin(10) + (m(7) * (t31 * t92 + t32 * t88) + (t54 * t88 - t55 * t92) * mrSges(7,3)) * pkin(5) + t100; m(7) * (t88 ^ 2 + t92 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t99 - t113; t108; t107; -Ifges(7,3) * t94 + t101; t100; Ifges(7,3) + t99; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
