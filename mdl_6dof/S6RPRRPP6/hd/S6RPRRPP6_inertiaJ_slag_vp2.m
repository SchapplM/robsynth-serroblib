% Calculate joint inertia matrix for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2018-11-23 16:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:14:15
% EndTime: 2018-11-23 16:14:16
% DurationCPUTime: 0.93s
% Computational Cost: add. (1095->268), mult. (2137->364), div. (0->0), fcn. (1969->6), ass. (0->98)
t129 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t94 = (-pkin(1) - pkin(7));
t128 = -2 * t94;
t127 = m(6) * pkin(4);
t122 = m(6) + m(7);
t126 = mrSges(7,2) + mrSges(6,3);
t88 = sin(pkin(9));
t89 = cos(pkin(9));
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t59 = t88 * t92 + t89 * t90;
t93 = cos(qJ(3));
t46 = t59 * t93;
t58 = t88 * t90 - t89 * t92;
t48 = t58 * t93;
t11 = t46 * mrSges(7,1) + t48 * mrSges(7,3);
t12 = t46 * mrSges(6,1) - t48 * mrSges(6,2);
t125 = -t11 - t12;
t17 = t58 * mrSges(7,1) - t59 * mrSges(7,3);
t18 = t58 * mrSges(6,1) + t59 * mrSges(6,2);
t124 = -t17 - t18;
t123 = 2 * qJ(2);
t121 = Ifges(5,4) * t90;
t120 = Ifges(5,4) * t92;
t119 = t90 * t93;
t91 = sin(qJ(3));
t118 = t91 * t94;
t117 = t92 * t93;
t116 = t93 * t94;
t114 = -qJ(5) - pkin(8);
t108 = qJ(5) * t93;
t64 = pkin(3) * t91 - pkin(8) * t93 + qJ(2);
t57 = t92 * t64;
t15 = -t92 * t108 + t57 + (-t90 * t94 + pkin(4)) * t91;
t33 = t92 * t118 + t90 * t64;
t23 = -t108 * t90 + t33;
t4 = t88 * t15 + t89 * t23;
t28 = -mrSges(7,2) * t46 + mrSges(7,3) * t91;
t29 = -mrSges(6,2) * t91 - mrSges(6,3) * t46;
t113 = t28 + t29;
t30 = mrSges(6,1) * t91 + mrSges(6,3) * t48;
t31 = -t91 * mrSges(7,1) - t48 * mrSges(7,2);
t112 = -t30 + t31;
t65 = -mrSges(5,1) * t92 + mrSges(5,2) * t90;
t111 = -t65 + mrSges(4,1);
t110 = t90 ^ 2 + t92 ^ 2;
t84 = t91 ^ 2;
t86 = t93 ^ 2;
t109 = t86 + t84;
t105 = t114 * t90;
t66 = t114 * t92;
t25 = -t105 * t89 - t66 * t88;
t27 = t105 * t88 - t89 * t66;
t107 = t25 ^ 2 + t27 ^ 2;
t75 = -pkin(4) * t92 - pkin(3);
t44 = t59 * t91;
t47 = t58 * t91;
t106 = t25 * t44 - t27 * t47;
t104 = t110 * mrSges(5,3);
t103 = t109 * mrSges(4,3);
t60 = pkin(4) * t119 - t116;
t101 = mrSges(5,1) * t90 + mrSges(5,2) * t92;
t3 = t15 * t89 - t23 * t88;
t32 = -t118 * t90 + t57;
t100 = -t32 * t90 + t33 * t92;
t98 = Ifges(5,5) * t117 + (-Ifges(7,4) - Ifges(6,5)) * t48 + (-Ifges(6,6) + Ifges(7,6)) * t46 + t129 * t91;
t95 = qJ(2) ^ 2;
t87 = t94 ^ 2;
t82 = Ifges(5,5) * t90;
t80 = Ifges(5,6) * t92;
t77 = t86 * t94;
t76 = t86 * t87;
t74 = -pkin(4) * t89 - pkin(5);
t72 = pkin(4) * t88 + qJ(6);
t68 = Ifges(5,1) * t90 + t120;
t67 = Ifges(5,2) * t92 + t121;
t63 = mrSges(5,1) * t91 - mrSges(5,3) * t117;
t62 = -mrSges(5,2) * t91 - mrSges(5,3) * t119;
t55 = t101 * t93;
t54 = Ifges(7,4) * t59;
t53 = Ifges(6,5) * t59;
t52 = Ifges(6,6) * t58;
t51 = Ifges(7,6) * t58;
t43 = Ifges(5,5) * t91 + (Ifges(5,1) * t92 - t121) * t93;
t42 = Ifges(5,6) * t91 + (-Ifges(5,2) * t90 + t120) * t93;
t22 = Ifges(6,1) * t59 - Ifges(6,4) * t58;
t21 = Ifges(7,1) * t59 + Ifges(7,5) * t58;
t20 = Ifges(6,4) * t59 - Ifges(6,2) * t58;
t19 = Ifges(7,5) * t59 + Ifges(7,3) * t58;
t14 = pkin(5) * t58 - qJ(6) * t59 + t75;
t10 = -Ifges(6,1) * t48 - Ifges(6,4) * t46 + Ifges(6,5) * t91;
t9 = -Ifges(7,1) * t48 + Ifges(7,4) * t91 + Ifges(7,5) * t46;
t8 = -Ifges(6,4) * t48 - Ifges(6,2) * t46 + Ifges(6,6) * t91;
t7 = -Ifges(7,5) * t48 + Ifges(7,6) * t91 + Ifges(7,3) * t46;
t5 = pkin(5) * t46 + qJ(6) * t48 + t60;
t2 = -pkin(5) * t91 - t3;
t1 = qJ(6) * t91 + t4;
t6 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t123) + 0.2e1 * t1 * t28 + 0.2e1 * t5 * t11 + 0.2e1 * t60 * t12 + 0.2e1 * t2 * t31 + 0.2e1 * t4 * t29 + 0.2e1 * t3 * t30 + 0.2e1 * t32 * t63 + 0.2e1 * t33 * t62 + Ifges(3,1) + Ifges(2,3) - (t9 + t10) * t48 + (t7 - t8) * t46 + t103 * t128 + ((mrSges(4,2) * t123) + Ifges(4,1) * t93 + t55 * t128 - t42 * t90 + t43 * t92) * t93 + (mrSges(4,1) * t123 + Ifges(4,2) * t91 + (-Ifges(5,6) * t90 - (2 * Ifges(4,4))) * t93 + t98) * t91 + m(4) * (t84 * t87 + t76 + t95) + (m(3) * (pkin(1) ^ 2 + t95)) + m(5) * (t32 ^ 2 + t33 ^ 2 + t76) + m(6) * (t3 ^ 2 + t4 ^ 2 + t60 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); -(m(3) * pkin(1)) + mrSges(3,2) + (t92 * t62 - t90 * t63) * t91 - t113 * t47 + t112 * t44 - t103 + (-t55 + t125) * t93 + m(6) * (-t3 * t44 - t4 * t47 - t60 * t93) + m(7) * (-t1 * t47 + t2 * t44 - t5 * t93) + m(5) * (t100 * t91 + t77) + m(4) * (t84 * t94 + t77); m(3) + m(5) * (t110 * t84 + t86) + m(4) * t109 + t122 * (t44 ^ 2 + t47 ^ 2 + t86); -pkin(3) * t55 + t14 * t11 + t75 * t12 + t5 * t17 + t60 * t18 + (-(t94 * mrSges(4,2)) + t82 / 0.2e1 + t80 / 0.2e1 + t53 / 0.2e1 - t52 / 0.2e1 + t54 / 0.2e1 + t51 / 0.2e1 - Ifges(4,6)) * t91 - (t21 / 0.2e1 + t22 / 0.2e1) * t48 + (t19 / 0.2e1 - t20 / 0.2e1) * t46 + t113 * t27 + t112 * t25 + (pkin(8) * t62 + t33 * mrSges(5,3) + t42 / 0.2e1) * t92 + (-pkin(8) * t63 - t32 * mrSges(5,3) + t43 / 0.2e1) * t90 + (t92 * t68 / 0.2e1 - t90 * t67 / 0.2e1 + Ifges(4,5) + t111 * t94) * t93 + m(5) * (pkin(3) * t116 + pkin(8) * t100) + m(6) * (-t25 * t3 + t27 * t4 + t60 * t75) + m(7) * (t1 * t27 + t14 * t5 + t25 * t2) + (t9 / 0.2e1 + t10 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t59 + (-t8 / 0.2e1 + t7 / 0.2e1 - t1 * mrSges(7,2) - t4 * mrSges(6,3)) * t58; (-mrSges(4,2) + t104) * t91 + (t111 + t124) * t93 + m(6) * (-t75 * t93 + t106) + m(7) * (-t14 * t93 + t106) + m(5) * (pkin(8) * t110 * t91 + pkin(3) * t93) + t126 * (t44 * t59 + t47 * t58); -0.2e1 * pkin(3) * t65 + 0.2e1 * t14 * t17 + 0.2e1 * t75 * t18 + t92 * t67 + t90 * t68 + Ifges(4,3) + m(6) * (t75 ^ 2 + t107) + m(7) * (t14 ^ 2 + t107) + m(5) * (pkin(8) ^ 2 * t110 + pkin(3) ^ 2) + (t21 + t22) * t59 + (t19 - t20) * t58 + 0.2e1 * pkin(8) * t104 + 0.2e1 * (t25 * t59 - t27 * t58) * t126; m(7) * (t1 * t72 + t2 * t74) + t98 + (t88 * t29 + t89 * t30 + m(6) * (t3 * t89 + t4 * t88)) * pkin(4) - Ifges(5,6) * t119 + t72 * t28 + t74 * t31 + t32 * mrSges(5,1) - t33 * mrSges(5,2) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2); -t101 * t91 - (-mrSges(6,2) + mrSges(7,3)) * t47 + (-mrSges(6,1) - mrSges(7,1)) * t44 + m(7) * (t44 * t74 - t47 * t72) + (-t44 * t89 - t47 * t88) * t127; -t25 * mrSges(7,1) + t27 * mrSges(7,3) - t27 * mrSges(6,2) - t25 * mrSges(6,1) + m(7) * (t25 * t74 + t27 * t72) + t80 + t54 + t51 + t53 - t52 + t82 - t101 * pkin(8) + (-t58 * t72 + t59 * t74) * mrSges(7,2) + (m(6) * (-t25 * t89 + t27 * t88) + (-t58 * t88 - t59 * t89) * mrSges(6,3)) * pkin(4); -0.2e1 * t74 * mrSges(7,1) + 0.2e1 * t72 * mrSges(7,3) + m(7) * (t72 ^ 2 + t74 ^ 2) + (0.2e1 * mrSges(6,1) * t89 - 0.2e1 * mrSges(6,2) * t88 + (t88 ^ 2 + t89 ^ 2) * t127) * pkin(4) + t129; m(6) * t60 + m(7) * t5 - t125; -t122 * t93; m(6) * t75 + m(7) * t14 - t124; 0; t122; m(7) * t2 + t31; m(7) * t44; m(7) * t25 + t59 * mrSges(7,2); m(7) * t74 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
