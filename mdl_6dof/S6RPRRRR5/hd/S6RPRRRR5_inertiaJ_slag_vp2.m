% Calculate joint inertia matrix for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:35:02
% EndTime: 2018-11-23 16:35:03
% DurationCPUTime: 0.93s
% Computational Cost: add. (3059->263), mult. (5780->381), div. (0->0), fcn. (6824->10), ass. (0->107)
t113 = cos(qJ(5));
t110 = sin(qJ(4));
t114 = cos(qJ(4));
t106 = sin(pkin(11));
t107 = cos(pkin(11));
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t79 = -t106 * t111 + t107 * t115;
t80 = t106 * t115 + t107 * t111;
t59 = t110 * t79 + t114 * t80;
t141 = t113 * t59;
t58 = t110 * t80 - t114 * t79;
t160 = Ifges(6,5) * t141 + Ifges(6,3) * t58;
t108 = sin(qJ(6));
t109 = sin(qJ(5));
t139 = Ifges(6,5) * t109 + Ifges(6,6) * t113;
t112 = cos(qJ(6));
t83 = -t108 * t109 + t112 * t113;
t150 = t83 * mrSges(7,3);
t159 = t108 * pkin(5) * t150 + t139;
t148 = pkin(7) + qJ(2);
t85 = t148 * t106;
t86 = t148 * t107;
t63 = -t111 * t86 - t115 * t85;
t121 = -pkin(8) * t80 + t63;
t64 = -t111 * t85 + t115 * t86;
t48 = pkin(8) * t79 + t64;
t30 = t110 * t48 - t114 * t121;
t158 = t30 ^ 2;
t157 = 0.2e1 * t30;
t84 = t108 * t113 + t109 * t112;
t60 = -t83 * mrSges(7,1) + t84 * mrSges(7,2);
t156 = 0.2e1 * t60;
t94 = -pkin(2) * t107 - pkin(1);
t68 = -pkin(3) * t79 + t94;
t155 = 0.2e1 * t68;
t154 = 0.2e1 * t79;
t103 = t107 ^ 2;
t152 = Ifges(6,6) * t58;
t151 = pkin(3) * t114;
t149 = t84 * mrSges(7,3);
t32 = t110 * t121 + t114 * t48;
t33 = pkin(4) * t58 - pkin(9) * t59 + t68;
t13 = t109 * t33 + t113 * t32;
t147 = Ifges(7,5) * t84 + Ifges(7,6) * t83;
t146 = Ifges(6,4) * t109;
t145 = Ifges(6,4) * t113;
t12 = -t109 * t32 + t113 * t33;
t144 = t109 * t12;
t143 = t109 * t59;
t142 = t113 * t13;
t95 = pkin(3) * t110 + pkin(9);
t140 = t113 * t95;
t138 = t106 ^ 2 + t103;
t137 = t109 ^ 2 + t113 ^ 2;
t136 = 0.2e1 * mrSges(7,3);
t36 = t84 * t59;
t37 = t83 * t59;
t135 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t58;
t134 = t112 * t149;
t97 = -pkin(5) * t113 - pkin(4);
t133 = -t79 * mrSges(4,1) + t80 * mrSges(4,2);
t132 = -t107 * mrSges(3,1) + t106 * mrSges(3,2);
t131 = t137 * t95;
t61 = Ifges(7,4) * t84 + Ifges(7,2) * t83;
t62 = Ifges(7,1) * t84 + Ifges(7,4) * t83;
t89 = Ifges(6,2) * t113 + t146;
t90 = Ifges(6,1) * t109 + t145;
t130 = t109 * t90 + t113 * t89 + t83 * t61 + t84 * t62 + Ifges(5,3);
t88 = -t113 * mrSges(6,1) + t109 * mrSges(6,2);
t129 = mrSges(6,1) * t109 + mrSges(6,2) * t113;
t128 = t142 - t144;
t75 = (-pkin(10) - t95) * t109;
t101 = t113 * pkin(10);
t76 = t101 + t140;
t56 = -t108 * t76 + t112 * t75;
t57 = t108 * t75 + t112 * t76;
t127 = t56 * mrSges(7,1) - t57 * mrSges(7,2) + t147;
t91 = (-pkin(10) - pkin(9)) * t109;
t92 = pkin(9) * t113 + t101;
t66 = -t108 * t92 + t112 * t91;
t67 = t108 * t91 + t112 * t92;
t126 = t66 * mrSges(7,1) - t67 * mrSges(7,2) + t147;
t125 = 0.2e1 * t137 * mrSges(6,3);
t5 = pkin(5) * t58 - pkin(10) * t141 + t12;
t6 = -pkin(10) * t143 + t13;
t3 = -t108 * t6 + t112 * t5;
t4 = t108 * t5 + t112 * t6;
t124 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t135;
t123 = (mrSges(5,1) * t114 - mrSges(5,2) * t110) * pkin(3);
t122 = (mrSges(7,1) * t112 - mrSges(7,2) * t108) * pkin(5);
t10 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t58;
t16 = pkin(5) * t143 + t30;
t21 = t152 + (-Ifges(6,2) * t109 + t145) * t59;
t22 = Ifges(6,5) * t58 + (Ifges(6,1) * t113 - t146) * t59;
t9 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t58;
t120 = -t32 * mrSges(5,2) - t3 * t149 + t4 * t150 + mrSges(6,3) * t142 + t16 * t60 + t109 * t22 / 0.2e1 + t113 * t21 / 0.2e1 - t36 * t61 / 0.2e1 + t37 * t62 / 0.2e1 - t89 * t143 / 0.2e1 + t90 * t141 / 0.2e1 - Ifges(5,6) * t58 + Ifges(5,5) * t59 + t83 * t9 / 0.2e1 + t84 * t10 / 0.2e1 + (t88 - mrSges(5,1)) * t30 + (t147 + t139) * t58 / 0.2e1;
t96 = -pkin(4) - t151;
t87 = t97 - t151;
t51 = t59 * mrSges(5,2);
t40 = mrSges(6,1) * t58 - mrSges(6,3) * t141;
t39 = -mrSges(6,2) * t58 - mrSges(6,3) * t143;
t38 = t129 * t59;
t18 = mrSges(7,1) * t58 - mrSges(7,3) * t37;
t17 = -mrSges(7,2) * t58 - mrSges(7,3) * t36;
t14 = mrSges(7,1) * t36 + mrSges(7,2) * t37;
t1 = [Ifges(4,2) * t79 ^ 2 + m(3) * (t138 * qJ(2) ^ 2 + pkin(1) ^ 2) - 0.2e1 * pkin(1) * t132 + 0.2e1 * t94 * t133 - t36 * t9 + t37 * t10 + 0.2e1 * t13 * t39 + 0.2e1 * t12 * t40 + 0.2e1 * t16 * t14 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t18 + (mrSges(5,1) * t155 - 0.2e1 * t32 * mrSges(5,3) + Ifges(5,2) * t58 + t135 + t160) * t58 + 0.2e1 * t138 * qJ(2) * mrSges(3,3) + m(4) * (t63 ^ 2 + t64 ^ 2 + t94 ^ 2) + m(7) * (t16 ^ 2 + t3 ^ 2 + t4 ^ 2) + Ifges(3,2) * t103 + (-0.2e1 * t63 * mrSges(4,3) + Ifges(4,1) * t80 + Ifges(4,4) * t154) * t80 + (mrSges(5,3) * t157 + Ifges(5,1) * t59 - 0.2e1 * Ifges(5,4) * t58 + t113 * t22 + (-t21 - t152) * t109) * t59 + m(5) * (t32 ^ 2 + t68 ^ 2 + t158) + m(6) * (t12 ^ 2 + t13 ^ 2 + t158) + Ifges(2,3) + t64 * mrSges(4,3) * t154 + t51 * t155 + t38 * t157 + (Ifges(3,1) * t106 + 0.2e1 * Ifges(3,4) * t107) * t106; -m(3) * pkin(1) + t58 * mrSges(5,1) + t109 * t39 + t113 * t40 + t84 * t17 + t83 * t18 + t51 + m(7) * (t3 * t83 + t4 * t84) + m(6) * (t109 * t13 + t113 * t12) + m(5) * t68 + m(4) * t94 + t132 + t133; m(3) + m(4) + m(5) + m(6) * t137 + m(7) * (t83 ^ 2 + t84 ^ 2); t120 + t39 * t140 + (m(5) * (t110 * t32 - t114 * t30) + (-t110 * t58 - t114 * t59) * mrSges(5,3)) * pkin(3) + m(6) * (t128 * t95 + t30 * t96) + t96 * t38 + t87 * t14 + Ifges(4,6) * t79 + Ifges(4,5) * t80 + t63 * mrSges(4,1) - t64 * mrSges(4,2) + t56 * t18 + t57 * t17 + (-t12 * mrSges(6,3) - t95 * t40) * t109 + m(7) * (t16 * t87 + t3 * t56 + t4 * t57); m(7) * (t56 * t83 + t57 * t84); t87 * t156 + 0.2e1 * t96 * t88 + Ifges(4,3) + 0.2e1 * t123 + (-t56 * t84 + t57 * t83) * t136 + t95 * t125 + m(7) * (t56 ^ 2 + t57 ^ 2 + t87 ^ 2) + m(6) * (t137 * t95 ^ 2 + t96 ^ 2) + m(5) * (t110 ^ 2 + t114 ^ 2) * pkin(3) ^ 2 + t130; (-t109 * t40 + t113 * t39) * pkin(9) + m(7) * (t16 * t97 + t3 * t66 + t4 * t67) + t120 - mrSges(6,3) * t144 + m(6) * (-pkin(4) * t30 + pkin(9) * t128) + t97 * t14 + t66 * t18 + t67 * t17 - pkin(4) * t38; m(7) * (t66 * t83 + t67 * t84); (t96 - pkin(4)) * t88 + (t87 + t97) * t60 + t123 + m(7) * (t56 * t66 + t57 * t67 + t87 * t97) + m(6) * (-pkin(4) * t96 + pkin(9) * t131) + ((-t56 - t66) * t84 + (t57 + t67) * t83) * mrSges(7,3) + (t137 * pkin(9) + t131) * mrSges(6,3) + t130; -0.2e1 * pkin(4) * t88 + t97 * t156 + (-t66 * t84 + t67 * t83) * t136 + pkin(9) * t125 + m(7) * (t66 ^ 2 + t67 ^ 2 + t97 ^ 2) + m(6) * (t137 * pkin(9) ^ 2 + pkin(4) ^ 2) + t130; -Ifges(6,6) * t143 + t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t108 * t4 + t112 * t3) + t108 * t17 + t112 * t18) * pkin(5) + t124 + t160; m(7) * (t108 * t84 + t112 * t83) * pkin(5) - t88 - t60; -t129 * t95 + (m(7) * (t108 * t57 + t112 * t56) - t134) * pkin(5) + t127 + t159; -t129 * pkin(9) + (m(7) * (t108 * t67 + t112 * t66) - t134) * pkin(5) + t126 + t159; Ifges(6,3) + Ifges(7,3) + m(7) * (t108 ^ 2 + t112 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t122; t124; -t60; t127; t126; Ifges(7,3) + t122; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
