% Calculate joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:10:00
% EndTime: 2018-11-23 15:10:01
% DurationCPUTime: 0.94s
% Computational Cost: add. (840->258), mult. (1866->353), div. (0->0), fcn. (1775->10), ass. (0->108)
t86 = sin(pkin(11));
t88 = cos(pkin(11));
t135 = t86 ^ 2 + t88 ^ 2;
t87 = sin(pkin(6));
t95 = cos(qJ(2));
t113 = t87 * t95;
t92 = sin(qJ(2));
t114 = t87 * t92;
t89 = cos(pkin(6));
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t41 = t114 * t94 + t89 * t91;
t15 = t113 * t88 + t41 * t86;
t17 = -t113 * t86 + t41 * t88;
t134 = t15 * t86 + t17 * t88;
t133 = 2 * pkin(8);
t132 = m(5) + m(6);
t112 = t88 * t91;
t115 = t86 * t91;
t44 = mrSges(5,1) * t115 + mrSges(5,2) * t112;
t90 = sin(qJ(6));
t93 = cos(qJ(6));
t49 = t86 * t93 - t88 * t90;
t35 = t49 * t91;
t98 = t86 * t90 + t88 * t93;
t36 = t98 * t91;
t7 = -t35 * mrSges(7,1) + t36 * mrSges(7,2);
t131 = t44 - t7;
t130 = mrSges(5,3) + mrSges(6,2);
t39 = t114 * t91 - t89 * t94;
t129 = t39 ^ 2;
t128 = m(5) * pkin(3);
t127 = pkin(4) + pkin(5);
t105 = qJ(5) * t86 + pkin(3);
t54 = -pkin(4) * t88 - t105;
t126 = m(6) * t54;
t125 = pkin(8) * t94;
t124 = Ifges(5,4) * t86;
t123 = Ifges(5,4) * t88;
t122 = Ifges(6,5) * t86;
t121 = Ifges(6,5) * t88;
t118 = t39 * t91;
t117 = t41 * t94;
t116 = t86 * mrSges(6,3);
t111 = -pkin(9) + qJ(4);
t50 = mrSges(5,2) * t94 - mrSges(5,3) * t115;
t53 = -mrSges(6,2) * t115 - mrSges(6,3) * t94;
t110 = t50 + t53;
t51 = -mrSges(5,1) * t94 - mrSges(5,3) * t112;
t52 = t94 * mrSges(6,1) + mrSges(6,2) * t112;
t109 = -t51 + t52;
t76 = t86 * mrSges(5,2);
t58 = -t88 * mrSges(5,1) + t76;
t108 = t58 - mrSges(4,1);
t55 = -pkin(3) * t94 - qJ(4) * t91 - pkin(2);
t27 = t88 * t125 + t86 * t55;
t107 = t135 * qJ(4) ^ 2;
t106 = Ifges(7,5) * t36 + Ifges(7,6) * t35 + Ifges(7,3) * t94;
t70 = t86 * t125;
t26 = t55 * t88 - t70;
t102 = t134 * qJ(4);
t22 = -qJ(5) * t94 + t27;
t11 = mrSges(7,1) * t98 + t49 * mrSges(7,2);
t80 = t94 * pkin(4);
t23 = -t26 + t80;
t100 = t22 * t88 + t23 * t86;
t99 = -t26 * t86 + t27 * t88;
t97 = pkin(8) ^ 2;
t85 = t94 ^ 2;
t84 = t91 ^ 2;
t82 = t87 ^ 2;
t79 = t84 * t97;
t72 = t82 * t95 ^ 2;
t66 = mrSges(6,1) * t115;
t65 = qJ(5) * t112;
t64 = -mrSges(4,1) * t94 + mrSges(4,2) * t91;
t63 = Ifges(5,1) * t86 + t123;
t62 = Ifges(6,1) * t86 - t121;
t61 = Ifges(5,2) * t88 + t124;
t60 = -Ifges(6,3) * t88 + t122;
t59 = t111 * t88;
t57 = -t88 * mrSges(6,1) - t116;
t56 = t111 * t86;
t46 = Ifges(7,5) * t49;
t45 = Ifges(7,6) * t98;
t43 = -mrSges(6,3) * t112 + t66;
t42 = t127 * t88 + t105;
t34 = -t65 + (pkin(4) * t86 + pkin(8)) * t91;
t33 = -Ifges(5,5) * t94 + (Ifges(5,1) * t88 - t124) * t91;
t32 = -Ifges(6,4) * t94 + (Ifges(6,1) * t88 + t122) * t91;
t31 = -Ifges(5,6) * t94 + (-Ifges(5,2) * t86 + t123) * t91;
t30 = -Ifges(6,6) * t94 + (Ifges(6,3) * t86 + t121) * t91;
t25 = mrSges(7,1) * t94 - mrSges(7,3) * t36;
t24 = -mrSges(7,2) * t94 + mrSges(7,3) * t35;
t20 = -t65 - (-t127 * t86 - pkin(8)) * t91;
t19 = t56 * t90 + t59 * t93;
t18 = t56 * t93 - t59 * t90;
t13 = Ifges(7,1) * t49 - Ifges(7,4) * t98;
t12 = Ifges(7,4) * t49 - Ifges(7,2) * t98;
t10 = pkin(9) * t115 + t22;
t8 = pkin(5) * t94 + t70 + t80 + (-pkin(9) * t91 - t55) * t88;
t6 = Ifges(7,1) * t36 + Ifges(7,4) * t35 + Ifges(7,5) * t94;
t5 = Ifges(7,4) * t36 + Ifges(7,2) * t35 + Ifges(7,6) * t94;
t4 = t15 * t90 + t17 * t93;
t3 = t15 * t93 - t17 * t90;
t2 = t10 * t93 + t8 * t90;
t1 = -t10 * t90 + t8 * t93;
t9 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t129) + m(4) * (t41 ^ 2 + t129 + t72) + m(3) * (t82 * t92 ^ 2 + t89 ^ 2 + t72) + t132 * (t15 ^ 2 + t17 ^ 2 + t129); mrSges(4,3) * t117 + t4 * t24 + t3 * t25 + t110 * t17 + t109 * t15 + (-t92 * mrSges(3,2) + (mrSges(3,1) - t64) * t95) * t87 + (t91 * mrSges(4,3) + t131 + t43) * t39 + m(7) * (t1 * t3 + t2 * t4 + t20 * t39) + m(5) * (pkin(8) * t118 - t15 * t26 + t17 * t27) + m(6) * (t15 * t23 + t17 * t22 + t34 * t39) + m(4) * (pkin(2) * t113 + (t117 + t118) * pkin(8)); -0.2e1 * pkin(2) * t64 + 0.2e1 * t1 * t25 + 0.2e1 * t2 * t24 - 0.2e1 * t20 * t7 + 0.2e1 * t22 * t53 + 0.2e1 * t23 * t52 + 0.2e1 * t26 * t51 + 0.2e1 * t27 * t50 + 0.2e1 * t34 * t43 + t35 * t5 + t36 * t6 + Ifges(3,3) + (t84 + t85) * mrSges(4,3) * t133 + m(4) * (pkin(2) ^ 2 + t85 * t97 + t79) + m(5) * (t26 ^ 2 + t27 ^ 2 + t79) + m(6) * (t22 ^ 2 + t23 ^ 2 + t34 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t20 ^ 2) + ((Ifges(6,2) + Ifges(5,3) + Ifges(4,2)) * t94 + t106) * t94 + (Ifges(4,1) * t91 + t44 * t133 + (t32 + t33) * t88 + (t30 - t31) * t86 + ((2 * Ifges(4,4)) + (-Ifges(6,4) - Ifges(5,5)) * t88 + (Ifges(5,6) - Ifges(6,6)) * t86) * t94) * t91; -t41 * mrSges(4,2) + (-t3 * t49 - t4 * t98) * mrSges(7,3) + (-t11 + t57 + t108) * t39 + m(7) * (t18 * t3 + t19 * t4 - t39 * t42) + m(5) * (-pkin(3) * t39 + t102) + m(6) * (t39 * t54 + t102) + t130 * t134; t42 * t7 - pkin(3) * t44 - t98 * t5 / 0.2e1 + t49 * t6 / 0.2e1 + t54 * t43 + t35 * t12 / 0.2e1 + t36 * t13 / 0.2e1 - t20 * t11 + t19 * t24 + t18 * t25 + (-t30 / 0.2e1 + t31 / 0.2e1) * t88 + (t32 / 0.2e1 + t33 / 0.2e1) * t86 + (-pkin(8) * mrSges(4,2) + t46 / 0.2e1 - t45 / 0.2e1 + Ifges(4,6) + (-Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t88 + (-Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t86) * t94 + (-t1 * t49 - t2 * t98) * mrSges(7,3) + t99 * mrSges(5,3) + t100 * mrSges(6,2) + m(7) * (t1 * t18 + t19 * t2 - t20 * t42) + (m(5) * t99 + m(6) * t100 + t109 * t86 + t110 * t88) * qJ(4) + (Ifges(4,5) + (t62 / 0.2e1 + t63 / 0.2e1) * t88 + (t60 / 0.2e1 - t61 / 0.2e1) * t86 + (t108 - t128) * pkin(8)) * t91 + (t57 + t126) * t34; -0.2e1 * pkin(3) * t58 + 0.2e1 * t42 * t11 - t98 * t12 + t49 * t13 + 0.2e1 * t54 * t57 + Ifges(4,3) + (-t60 + t61) * t88 + (t62 + t63) * t86 + 0.2e1 * (-t18 * t49 - t19 * t98) * mrSges(7,3) + m(7) * (t18 ^ 2 + t19 ^ 2 + t42 ^ 2) + m(5) * (pkin(3) ^ 2 + t107) + m(6) * (t54 ^ 2 + t107) + 0.2e1 * t130 * qJ(4) * t135; 0.2e1 * (m(7) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t39; t66 + (m(5) * pkin(8) - mrSges(6,3) * t88) * t91 + m(6) * t34 + m(7) * t20 + t131; -t128 - t116 + t76 + (-mrSges(6,1) - mrSges(5,1)) * t88 + t126 - m(7) * t42 - t11; m(7) + t132; m(7) * (t3 * t93 + t4 * t90) + m(6) * t15; t90 * t24 + t93 * t25 + m(7) * (t1 * t93 + t2 * t90) + m(6) * t23 + t52; m(7) * (t18 * t93 + t19 * t90) + (m(6) * qJ(4) + mrSges(6,2)) * t86 + (-t49 * t93 - t90 * t98) * mrSges(7,3); 0; m(6) + m(7) * (t90 ^ 2 + t93 ^ 2); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t106; mrSges(7,1) * t18 - mrSges(7,2) * t19 - t45 + t46; 0; mrSges(7,1) * t93 - t90 * mrSges(7,2); Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
