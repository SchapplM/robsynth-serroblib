% Calculate joint inertia matrix for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2018-11-23 16:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:11:16
% EndTime: 2018-11-23 16:11:17
% DurationCPUTime: 0.99s
% Computational Cost: add. (1159->260), mult. (2295->359), div. (0->0), fcn. (2135->8), ass. (0->98)
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t100 = mrSges(5,1) * t90 + mrSges(5,2) * t92;
t91 = sin(qJ(3));
t56 = t100 * t91;
t86 = sin(pkin(10));
t88 = cos(pkin(10));
t59 = t86 * t92 + t88 * t90;
t44 = t59 * t91;
t58 = t86 * t90 - t88 * t92;
t46 = t58 * t91;
t14 = mrSges(7,1) * t44 + t46 * mrSges(7,3);
t15 = t44 * mrSges(6,1) - mrSges(6,2) * t46;
t97 = -t14 - t15;
t128 = -t56 + t97;
t87 = sin(pkin(9));
t75 = pkin(1) * t87 + pkin(7);
t127 = 0.2e1 * t75;
t126 = m(6) * pkin(4);
t123 = m(6) + m(7);
t125 = mrSges(7,2) + mrSges(6,3);
t18 = mrSges(7,1) * t58 - t59 * mrSges(7,3);
t19 = t58 * mrSges(6,1) + mrSges(6,2) * t59;
t124 = -t18 - t19;
t93 = cos(qJ(3));
t122 = pkin(3) * t93;
t121 = pkin(8) * t91;
t120 = Ifges(5,4) * t90;
t119 = Ifges(5,4) * t92;
t118 = t75 * t93;
t117 = t90 * t91;
t69 = t91 * t75;
t116 = t91 * t92;
t64 = -mrSges(5,1) * t92 + mrSges(5,2) * t90;
t115 = -mrSges(4,1) + t64;
t113 = -qJ(5) - pkin(8);
t108 = qJ(5) * t91;
t89 = cos(pkin(9));
t77 = -pkin(1) * t89 - pkin(2);
t57 = -t121 + t77 - t122;
t48 = t92 * t57;
t12 = -t92 * t108 + t48 + (-t75 * t90 - pkin(4)) * t93;
t25 = t118 * t92 + t57 * t90;
t16 = -t108 * t90 + t25;
t4 = t12 * t86 + t16 * t88;
t30 = -mrSges(7,2) * t44 - mrSges(7,3) * t93;
t31 = mrSges(6,2) * t93 - mrSges(6,3) * t44;
t112 = t30 + t31;
t32 = -mrSges(6,1) * t93 + mrSges(6,3) * t46;
t33 = mrSges(7,1) * t93 - mrSges(7,2) * t46;
t111 = -t32 + t33;
t49 = pkin(4) * t117 + t69;
t110 = t90 ^ 2 + t92 ^ 2;
t83 = t91 ^ 2;
t85 = t93 ^ 2;
t109 = t83 + t85;
t104 = t113 * t90;
t65 = t113 * t92;
t27 = -t104 * t88 - t65 * t86;
t29 = t104 * t86 - t88 * t65;
t107 = t27 ^ 2 + t29 ^ 2;
t106 = -Ifges(6,3) - Ifges(5,3) - Ifges(7,2);
t78 = -pkin(4) * t92 - pkin(3);
t105 = t27 * t44 - t29 * t46;
t103 = t110 * mrSges(5,3);
t101 = -Ifges(5,5) * t116 + (Ifges(7,4) + Ifges(6,5)) * t46 + (Ifges(6,6) - Ifges(7,6)) * t44;
t3 = t12 * t88 - t16 * t86;
t24 = -t118 * t90 + t48;
t99 = -t24 * t90 + t25 * t92;
t81 = Ifges(5,5) * t90;
t80 = Ifges(5,6) * t92;
t76 = -pkin(4) * t88 - pkin(5);
t73 = pkin(4) * t86 + qJ(6);
t72 = t75 ^ 2;
t68 = t83 * t72;
t67 = Ifges(5,1) * t90 + t119;
t66 = Ifges(5,2) * t92 + t120;
t62 = -mrSges(5,1) * t93 - mrSges(5,3) * t116;
t61 = mrSges(5,2) * t93 - mrSges(5,3) * t117;
t55 = Ifges(7,4) * t59;
t54 = Ifges(6,5) * t59;
t53 = Ifges(6,6) * t58;
t52 = Ifges(7,6) * t58;
t43 = -Ifges(5,5) * t93 + (Ifges(5,1) * t92 - t120) * t91;
t42 = -Ifges(5,6) * t93 + (-Ifges(5,2) * t90 + t119) * t91;
t23 = Ifges(6,1) * t59 - Ifges(6,4) * t58;
t22 = Ifges(7,1) * t59 + Ifges(7,5) * t58;
t21 = Ifges(6,4) * t59 - Ifges(6,2) * t58;
t20 = Ifges(7,5) * t59 + Ifges(7,3) * t58;
t17 = pkin(5) * t58 - qJ(6) * t59 + t78;
t10 = -Ifges(6,1) * t46 - Ifges(6,4) * t44 - Ifges(6,5) * t93;
t9 = -Ifges(7,1) * t46 - Ifges(7,4) * t93 + Ifges(7,5) * t44;
t8 = -Ifges(6,4) * t46 - Ifges(6,2) * t44 - Ifges(6,6) * t93;
t7 = -Ifges(7,5) * t46 - Ifges(7,6) * t93 + Ifges(7,3) * t44;
t5 = pkin(5) * t44 + qJ(6) * t46 + t49;
t2 = pkin(5) * t93 - t3;
t1 = -qJ(6) * t93 + t4;
t6 = [0.2e1 * t1 * t30 + 0.2e1 * t5 * t14 + 0.2e1 * t49 * t15 + 0.2e1 * t2 * t33 + 0.2e1 * t24 * t62 + 0.2e1 * t25 * t61 + 0.2e1 * t3 * t32 + 0.2e1 * t4 * t31 + Ifges(2,3) + Ifges(3,3) - (t9 + t10) * t46 + (t7 - t8) * t44 + (0.2e1 * mrSges(4,2) * t77 + Ifges(4,1) * t91 + t56 * t127 - t42 * t90 + t43 * t92) * t91 + (-0.2e1 * t77 * mrSges(4,1) + (Ifges(4,2) - t106) * t93 + (Ifges(5,6) * t90 + (2 * Ifges(4,4))) * t91 + t101) * t93 + m(4) * (t72 * t85 + t77 ^ 2 + t68) + m(5) * (t24 ^ 2 + t25 ^ 2 + t68) + m(6) * (t3 ^ 2 + t4 ^ 2 + t49 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (t87 ^ 2 + t89 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t89 - mrSges(3,2) * t87) * pkin(1) + t109 * mrSges(4,3) * t127; -t112 * t46 + t111 * t44 + t128 * t93 + m(7) * (-t1 * t46 + t2 * t44 - t5 * t93) + m(6) * (-t3 * t44 - t4 * t46 - t49 * t93) + (-t90 * t62 + t92 * t61 + m(5) * (t99 - t118)) * t91; m(3) + m(5) * (t110 * t83 + t85) + m(4) * t109 + t123 * (t44 ^ 2 + t46 ^ 2 + t85); -pkin(3) * t56 + t17 * t14 + t78 * t15 + t5 * t18 + t49 * t19 + (-t81 / 0.2e1 - t80 / 0.2e1 - t54 / 0.2e1 + t53 / 0.2e1 - t55 / 0.2e1 - t52 / 0.2e1 + Ifges(4,6) - t75 * mrSges(4,2)) * t93 - (t22 / 0.2e1 + t23 / 0.2e1) * t46 + (t20 / 0.2e1 - t21 / 0.2e1) * t44 + t112 * t29 + t111 * t27 + (t42 / 0.2e1 + pkin(8) * t61 + t25 * mrSges(5,3)) * t92 + (t43 / 0.2e1 - pkin(8) * t62 - t24 * mrSges(5,3)) * t90 + (Ifges(4,5) + t92 * t67 / 0.2e1 - t90 * t66 / 0.2e1 + t115 * t75) * t91 + m(5) * (-pkin(3) * t69 + pkin(8) * t99) + m(6) * (-t27 * t3 + t29 * t4 + t49 * t78) + m(7) * (t1 * t29 + t17 * t5 + t2 * t27) + (t9 / 0.2e1 + t10 / 0.2e1 - t3 * mrSges(6,3) + t2 * mrSges(7,2)) * t59 + (t7 / 0.2e1 - t8 / 0.2e1 - t1 * mrSges(7,2) - t4 * mrSges(6,3)) * t58; (-mrSges(4,2) + t103) * t91 + (-t115 + t124) * t93 + m(7) * (-t17 * t93 + t105) + m(5) * (t110 * t121 + t122) + m(6) * (-t78 * t93 + t105) + t125 * (t44 * t59 + t46 * t58); -0.2e1 * pkin(3) * t64 + 0.2e1 * t17 * t18 + 0.2e1 * t78 * t19 + t92 * t66 + t90 * t67 + Ifges(4,3) + m(6) * (t78 ^ 2 + t107) + m(7) * (t17 ^ 2 + t107) + m(5) * (pkin(8) ^ 2 * t110 + pkin(3) ^ 2) + (t22 + t23) * t59 + (t20 - t21) * t58 + 0.2e1 * pkin(8) * t103 + 0.2e1 * (t27 * t59 - t29 * t58) * t125; -Ifges(5,6) * t117 - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t76 * t33 + t1 * mrSges(7,3) + t73 * t30 + m(7) * (t1 * t73 + t2 * t76) + t24 * mrSges(5,1) - t25 * mrSges(5,2) + t106 * t93 + (m(6) * (t3 * t88 + t4 * t86) + t86 * t31 + t88 * t32) * pkin(4) - t101; m(7) * (t44 * t76 - t46 * t73) + (-t44 * t88 - t46 * t86) * t126 + t128; t29 * mrSges(7,3) - t27 * mrSges(7,1) - t29 * mrSges(6,2) - t27 * mrSges(6,1) + m(7) * (t27 * t76 + t29 * t73) + t55 + t52 + t54 - t53 + t81 + t80 - t100 * pkin(8) + (-t58 * t73 + t59 * t76) * mrSges(7,2) + (m(6) * (-t27 * t88 + t29 * t86) + (-t58 * t86 - t59 * t88) * mrSges(6,3)) * pkin(4); -0.2e1 * t76 * mrSges(7,1) + 0.2e1 * t73 * mrSges(7,3) + m(7) * (t73 ^ 2 + t76 ^ 2) - t106 + (0.2e1 * mrSges(6,1) * t88 - 0.2e1 * mrSges(6,2) * t86 + (t86 ^ 2 + t88 ^ 2) * t126) * pkin(4); m(6) * t49 + m(7) * t5 - t97; -t123 * t93; m(6) * t78 + m(7) * t17 - t124; 0; t123; m(7) * t2 + t33; m(7) * t44; m(7) * t27 + t59 * mrSges(7,2); m(7) * t76 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
