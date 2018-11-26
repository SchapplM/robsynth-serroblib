% Calculate joint inertia matrix for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:15:03
% EndTime: 2018-11-23 17:15:04
% DurationCPUTime: 1.27s
% Computational Cost: add. (1278->252), mult. (2295->323), div. (0->0), fcn. (2084->6), ass. (0->94)
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t109 = t78 ^ 2 + t81 ^ 2;
t48 = -mrSges(6,1) * t81 + mrSges(6,2) * t78;
t115 = t48 - mrSges(5,1);
t143 = (mrSges(7,2) + mrSges(6,3)) * t109;
t47 = mrSges(7,1) * t81 + mrSges(7,3) * t78;
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t153 = (-mrSges(5,2) + t143) * t79 + (t47 - t115) * t82;
t152 = Ifges(6,5) + Ifges(7,4);
t151 = Ifges(7,2) + Ifges(6,3);
t80 = sin(qJ(2));
t83 = cos(qJ(2));
t150 = t80 ^ 2 + t83 ^ 2;
t148 = t109 * t79;
t144 = -m(4) * pkin(2) - mrSges(4,1);
t134 = pkin(7) - pkin(8);
t107 = t134 * t80;
t55 = t134 * t83;
t19 = -t107 * t82 + t55 * t79;
t21 = t107 * t79 + t82 * t55;
t38 = -t79 * t83 + t80 * t82;
t92 = pkin(5) * t78 - qJ(6) * t81;
t5 = t38 * t92 + t19;
t142 = -t21 * mrSges(5,2) + t115 * t19 - t5 * t47;
t64 = Ifges(7,6) * t81;
t141 = (qJ(6) * mrSges(7,2) + Ifges(6,6)) * t81 - (mrSges(7,2) * pkin(5) - t152) * t78 - t64;
t130 = Ifges(7,5) * t78;
t49 = -Ifges(7,3) * t81 + t130;
t132 = Ifges(6,4) * t78;
t52 = Ifges(6,2) * t81 + t132;
t129 = Ifges(7,5) * t81;
t53 = Ifges(7,1) * t78 - t129;
t131 = Ifges(6,4) * t81;
t54 = Ifges(6,1) * t78 + t131;
t140 = (t53 / 0.2e1 + t54 / 0.2e1) * t81 + (t49 / 0.2e1 - t52 / 0.2e1) * t78 + Ifges(5,5);
t88 = (t53 + t54) * t78 + (t52 - t49) * t81 + Ifges(5,3);
t139 = t19 ^ 2;
t138 = 0.2e1 * t19;
t46 = -pkin(2) * t83 - qJ(3) * t80 - pkin(1);
t34 = pkin(3) * t83 - t46;
t137 = 0.2e1 * t34;
t136 = -0.2e1 * t46;
t135 = -0.2e1 * t48;
t37 = t79 * t80 + t82 * t83;
t128 = Ifges(6,6) * t37;
t127 = t19 * t82;
t126 = t38 * t78;
t125 = t38 * t81;
t84 = -pkin(2) - pkin(3);
t43 = -qJ(3) * t79 + t82 * t84;
t124 = t43 * mrSges(5,1);
t44 = qJ(3) * t82 + t79 * t84;
t123 = t44 * mrSges(5,2);
t14 = -mrSges(6,2) * t37 - mrSges(6,3) * t126;
t17 = -mrSges(7,2) * t126 + mrSges(7,3) * t37;
t120 = t14 + t17;
t15 = mrSges(6,1) * t37 - mrSges(6,3) * t125;
t16 = -mrSges(7,1) * t37 + mrSges(7,2) * t125;
t119 = -t15 + t16;
t11 = pkin(4) * t37 - pkin(9) * t38 + t34;
t4 = t11 * t78 + t21 * t81;
t42 = -pkin(9) + t44;
t118 = t148 * t42;
t117 = t109 * pkin(9) * t42;
t116 = t109 * t42 ^ 2;
t112 = t148 * pkin(9);
t111 = t109 * pkin(9) ^ 2;
t110 = t150 * pkin(7) ^ 2;
t103 = Ifges(6,6) * t81 / 0.2e1 - t64 / 0.2e1 - Ifges(5,6) + t152 * t78 / 0.2e1;
t98 = Ifges(7,6) * t126 + t125 * t152 + t151 * t37;
t1 = qJ(6) * t37 + t4;
t3 = t11 * t81 - t21 * t78;
t2 = -pkin(5) * t37 - t3;
t97 = t1 * t81 + t2 * t78;
t96 = -t3 * t78 + t4 * t81;
t94 = t78 * mrSges(6,1) + t81 * mrSges(6,2);
t93 = t78 * mrSges(7,1) - t81 * mrSges(7,3);
t45 = -pkin(5) * t81 - qJ(6) * t78 - pkin(4);
t6 = Ifges(7,6) * t37 + (Ifges(7,3) * t78 + t129) * t38;
t7 = t128 + (-Ifges(6,2) * t78 + t131) * t38;
t91 = -t1 * mrSges(7,2) - t4 * mrSges(6,3) + t6 / 0.2e1 - t7 / 0.2e1;
t8 = Ifges(7,4) * t37 + (Ifges(7,1) * t81 + t130) * t38;
t9 = Ifges(6,5) * t37 + (Ifges(6,1) * t81 - t132) * t38;
t90 = -t2 * mrSges(7,2) + t3 * mrSges(6,3) - t8 / 0.2e1 - t9 / 0.2e1;
t87 = -m(7) * t92 - t93 - t94;
t76 = t82 ^ 2;
t73 = t79 ^ 2;
t41 = pkin(4) - t43;
t22 = -t43 - t45;
t13 = t94 * t38;
t12 = t93 * t38;
t10 = [0.2e1 * t1 * t17 + 0.2e1 * t5 * t12 + t13 * t138 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t136 + (Ifges(4,3) + Ifges(3,2)) * t83) * t83 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t136 + (Ifges(4,1) + Ifges(3,1)) * t80 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t83) * t80 + (mrSges(5,1) * t137 - 0.2e1 * mrSges(5,3) * t21 + Ifges(5,2) * t37 + t98) * t37 + (mrSges(5,2) * t137 + mrSges(5,3) * t138 + Ifges(5,1) * t38 - 0.2e1 * Ifges(5,4) * t37 + (t8 + t9) * t81 + (t6 - t7 - t128) * t78) * t38 + m(4) * (t46 ^ 2 + t110) + m(3) * (pkin(1) ^ 2 + t110) + m(5) * (t21 ^ 2 + t34 ^ 2 + t139) + m(6) * (t3 ^ 2 + t4 ^ 2 + t139) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t150; t22 * t12 + t41 * t13 + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t83 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t80 + (t120 * t42 + t91) * t81 + (t119 * t42 + t90) * t78 + m(6) * (t19 * t41 + t42 * t96) + m(7) * (t22 * t5 + t42 * t97) + m(5) * (-t19 * t43 + t21 * t44) + (-t44 * mrSges(5,3) - t103) * t37 + (-t43 * mrSges(5,3) - t140) * t38 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t83 + (-mrSges(3,1) + t144) * t80) * pkin(7) - t142; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t124 + 0.2e1 * t123 + 0.2e1 * qJ(3) * mrSges(4,3) + 0.2e1 * t22 * t47 + t41 * t135 + Ifges(4,2) + Ifges(3,3) + m(7) * (t22 ^ 2 + t116) + m(6) * (t41 ^ 2 + t116) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t88 - 0.2e1 * t42 * t143; (m(4) * pkin(7) + mrSges(4,2)) * t80 + (-t38 * mrSges(5,3) - t12 - t13) * t82 + (-t37 * mrSges(5,3) + t119 * t78 + t120 * t81) * t79 + m(7) * (-t5 * t82 + t79 * t97) + m(6) * (t79 * t96 - t127) + m(5) * (t21 * t79 - t127); m(7) * (-t22 * t82 + t118) + m(6) * (-t41 * t82 + t118) + m(5) * (t43 * t82 + t44 * t79) + t144 - t153; m(4) + m(5) * (t73 + t76) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t109 * t73 + t76); -pkin(4) * t13 + t45 * t12 + (pkin(9) * t120 - t91) * t81 + (pkin(9) * t119 - t90) * t78 + m(7) * (pkin(9) * t97 + t45 * t5) + m(6) * (-pkin(4) * t19 + pkin(9) * t96) + t103 * t37 + t140 * t38 + t142; t124 - t123 + (pkin(4) + t41) * t48 + (t45 - t22) * t47 + m(7) * (t22 * t45 + t117) + m(6) * (-pkin(4) * t41 + t117) + (-pkin(9) + t42) * t143 - t88; m(6) * (pkin(4) * t82 + t112) + m(7) * (-t45 * t82 + t112) + t153; pkin(4) * t135 - 0.2e1 * t45 * t47 + m(6) * (pkin(4) ^ 2 + t111) + m(7) * (t45 ^ 2 + t111) + t88 + 0.2e1 * pkin(9) * t143; -Ifges(6,6) * t126 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t17 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - pkin(5) * t16 - t2 * mrSges(7,1) + t98; t42 * t87 - t141; t87 * t79; pkin(9) * t87 + t141; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t151; m(7) * t2 + t16; (m(7) * t42 - mrSges(7,2)) * t78; m(7) * t78 * t79; (m(7) * pkin(9) + mrSges(7,2)) * t78; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
