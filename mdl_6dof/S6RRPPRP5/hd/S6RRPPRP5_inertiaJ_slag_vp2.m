% Calculate joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2018-11-23 16:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:47:30
% EndTime: 2018-11-23 16:47:31
% DurationCPUTime: 0.96s
% Computational Cost: add. (1106->247), mult. (2004->316), div. (0->0), fcn. (1842->6), ass. (0->90)
t127 = m(6) + m(7);
t131 = pkin(3) + pkin(7);
t130 = (-mrSges(7,2) - mrSges(6,3));
t129 = Ifges(7,2) + Ifges(6,3);
t95 = sin(qJ(2));
t96 = cos(qJ(2));
t128 = t95 ^ 2 + t96 ^ 2;
t126 = -mrSges(6,2) + mrSges(7,3);
t121 = cos(qJ(5));
t92 = cos(pkin(9));
t108 = t121 * t92;
t91 = sin(pkin(9));
t94 = sin(qJ(5));
t58 = t91 * t94 - t108;
t116 = t94 * t92;
t99 = -t121 * t91 - t116;
t110 = t58 ^ 2 + t99 ^ 2;
t125 = -m(4) * pkin(2) + mrSges(4,2);
t124 = -m(7) * pkin(5) - mrSges(7,1);
t123 = -t91 / 0.2e1;
t93 = -pkin(2) - qJ(4);
t122 = -pkin(8) + t93;
t106 = -qJ(3) * t95 - pkin(1);
t55 = t93 * t96 + t106;
t70 = t131 * t95;
t62 = t92 * t70;
t12 = pkin(4) * t95 + t62 + (pkin(8) * t96 - t55) * t91;
t117 = t92 * t96;
t19 = t92 * t55 + t91 * t70;
t16 = -pkin(8) * t117 + t19;
t4 = t94 * t12 + t121 * t16;
t120 = Ifges(5,4) * t91;
t119 = Ifges(5,4) * t92;
t118 = t91 * t96;
t43 = -t108 * t96 + t118 * t94;
t30 = mrSges(7,2) * t43 + mrSges(7,3) * t95;
t31 = -mrSges(6,2) * t95 + mrSges(6,3) * t43;
t115 = t30 + t31;
t44 = t99 * t96;
t32 = mrSges(6,1) * t95 - mrSges(6,3) * t44;
t33 = -t95 * mrSges(7,1) + t44 * mrSges(7,2);
t114 = -t32 + t33;
t67 = t91 * mrSges(5,1) + t92 * mrSges(5,2);
t113 = t128 * pkin(7) ^ 2;
t71 = t131 * t96;
t112 = t91 ^ 2 + t92 ^ 2;
t75 = t91 * pkin(4) + qJ(3);
t65 = t122 * t91;
t21 = -t108 * t122 + t65 * t94;
t23 = t116 * t122 + t121 * t65;
t111 = t21 ^ 2 + t23 ^ 2;
t46 = pkin(4) * t117 + t71;
t109 = m(5) * t112;
t107 = t112 * mrSges(5,3);
t15 = -t43 * mrSges(6,1) + t44 * mrSges(6,2);
t25 = -mrSges(6,1) * t99 - t58 * mrSges(6,2);
t14 = -t43 * mrSges(7,1) - t44 * mrSges(7,3);
t24 = -mrSges(7,1) * t99 + t58 * mrSges(7,3);
t104 = 2 * t130;
t47 = mrSges(5,1) * t117 - mrSges(5,2) * t118;
t102 = pkin(5) * t58 + qJ(6) * t99;
t18 = -t55 * t91 + t62;
t101 = t18 * t92 + t19 * t91;
t100 = t129 * t95 + (Ifges(7,4) + Ifges(6,5)) * t44 + (Ifges(6,6) - Ifges(7,6)) * t43;
t3 = t12 * t121 - t94 * t16;
t97 = qJ(3) ^ 2;
t69 = Ifges(5,1) * t92 - t120;
t68 = -Ifges(5,2) * t91 + t119;
t66 = -pkin(2) * t96 + t106;
t64 = -mrSges(5,2) * t95 - mrSges(5,3) * t117;
t63 = mrSges(5,1) * t95 + mrSges(5,3) * t118;
t54 = Ifges(7,4) * t58;
t53 = Ifges(6,5) * t58;
t52 = Ifges(6,6) * t99;
t51 = Ifges(7,6) * t99;
t42 = Ifges(5,5) * t95 + (-Ifges(5,1) * t91 - t119) * t96;
t41 = Ifges(5,6) * t95 + (-Ifges(5,2) * t92 - t120) * t96;
t29 = -Ifges(6,1) * t58 + Ifges(6,4) * t99;
t28 = -Ifges(7,1) * t58 - Ifges(7,5) * t99;
t27 = -Ifges(6,4) * t58 + Ifges(6,2) * t99;
t26 = -Ifges(7,5) * t58 - Ifges(7,3) * t99;
t17 = -pkin(5) * t99 + qJ(6) * t58 + t75;
t10 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t95;
t9 = Ifges(7,1) * t44 + Ifges(7,4) * t95 - Ifges(7,5) * t43;
t8 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t95;
t7 = Ifges(7,5) * t44 + Ifges(7,6) * t95 - Ifges(7,3) * t43;
t5 = -pkin(5) * t43 - qJ(6) * t44 + t46;
t2 = -t95 * pkin(5) - t3;
t1 = qJ(6) * t95 + t4;
t6 = [0.2e1 * t1 * t30 + 0.2e1 * t5 * t14 + 0.2e1 * t46 * t15 + 0.2e1 * t18 * t63 + 0.2e1 * t19 * t64 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t32 + 0.2e1 * t4 * t31 + 0.2e1 * t71 * t47 + Ifges(2,3) + (t9 + t10) * t44 + (t8 - t7) * t43 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t66 * mrSges(4,2) - t92 * t41 - t91 * t42 + (Ifges(3,2) + Ifges(4,3)) * t96) * t96 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t66 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2) + Ifges(3,1)) * t95 + (-Ifges(5,5) * t91 - Ifges(5,6) * t92 + (2 * Ifges(3,4)) + (2 * Ifges(4,6))) * t96 + t100) * t95 + m(5) * (t18 ^ 2 + t19 ^ 2 + t71 ^ 2) + m(4) * (t66 ^ 2 + t113) + m(3) * (pkin(1) ^ 2 + t113) + m(6) * (t3 ^ 2 + t4 ^ 2 + t46 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t128; qJ(3) * t47 + t17 * t14 + t75 * t15 + t5 * t24 + t46 * t25 + t71 * t67 + (t28 / 0.2e1 + t29 / 0.2e1) * t44 + (t27 / 0.2e1 - t26 / 0.2e1) * t43 + t115 * t23 + t114 * t21 + (t42 / 0.2e1 - t18 * mrSges(5,3) + t93 * t63) * t92 + (-t41 / 0.2e1 - t19 * mrSges(5,3) + t93 * t64) * t91 + (Ifges(5,5) * t92 / 0.2e1 + Ifges(5,6) * t123 - t53 / 0.2e1 + t52 / 0.2e1 - t54 / 0.2e1 - t51 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1)) * t95 + m(5) * (qJ(3) * t71 + t101 * t93) + m(6) * (-t21 * t3 + t23 * t4 + t75 * t46) + m(7) * (t23 * t1 + t17 * t5 + t21 * t2) - (t7 / 0.2e1 - t8 / 0.2e1 - t1 * mrSges(7,2) - t4 * mrSges(6,3)) * t99 + (-t9 / 0.2e1 - t10 / 0.2e1 - t2 * mrSges(7,2) + t3 * mrSges(6,3)) * t58 + (-Ifges(4,5) + Ifges(3,6) + t69 * t123 - t92 * t68 / 0.2e1 + qJ(3) * mrSges(4,1)) * t96 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t96 + (-mrSges(3,1) + t125) * t95) * pkin(7); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t17 * t24 + 0.2e1 * t75 * t25 - t91 * t68 + t92 * t69 + Ifges(4,1) + Ifges(3,3) - (t104 * t23 + t26 - t27) * t99 + (t104 * t21 - t28 - t29) * t58 + m(6) * (t75 ^ 2 + t111) + m(7) * (t17 ^ 2 + t111) + m(5) * (t112 * t93 ^ 2 + t97) + m(4) * (pkin(2) ^ 2 + t97) + 0.2e1 * (mrSges(4,3) + t67) * qJ(3) - 0.2e1 * t93 * t107; t92 * t63 + t91 * t64 + (m(4) * pkin(7) + mrSges(4,1)) * t95 - t115 * t99 + t114 * t58 + m(7) * (-t1 * t99 + t2 * t58) + m(6) * (-t3 * t58 - t4 * t99) + m(5) * t101; -t107 + t93 * t109 + t127 * (t21 * t58 - t23 * t99) + t125 + t130 * t110; t127 * t110 + m(4) + t109; m(5) * t71 + m(6) * t46 + m(7) * t5 + t14 + t15 + t47; m(5) * qJ(3) + m(6) * t75 + m(7) * t17 + t24 + t25 + t67; 0; m(5) + t127; m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) + qJ(6) * t30 - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t33 + t100; t102 * mrSges(7,2) - t51 + t52 - t53 - t54 + (m(7) * qJ(6) + t126) * t23 + (-mrSges(6,1) + t124) * t21; -m(7) * t102 - t126 * t99 + (-mrSges(6,1) - mrSges(7,1)) * t58; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t129; m(7) * t2 + t33; m(7) * t21 - t58 * mrSges(7,2); m(7) * t58; 0; t124; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
