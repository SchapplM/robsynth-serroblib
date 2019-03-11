% Calculate joint inertia matrix for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:31
% EndTime: 2019-03-09 09:28:34
% DurationCPUTime: 1.12s
% Computational Cost: add. (1290->290), mult. (2746->386), div. (0->0), fcn. (2574->8), ass. (0->115)
t151 = Ifges(3,5) + Ifges(5,5);
t150 = Ifges(4,1) + Ifges(5,1) + Ifges(3,3);
t96 = (qJ(3) - pkin(9));
t149 = -2 * t96;
t97 = (pkin(2) + qJ(4));
t148 = 2 * t97;
t102 = cos(qJ(6));
t104 = cos(qJ(2));
t94 = sin(pkin(6));
t120 = t94 * t104;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t101 = sin(qJ(2));
t123 = t101 * t94;
t95 = cos(pkin(6));
t44 = t100 * t123 + t103 * t95;
t99 = sin(qJ(6));
t26 = t102 * t120 - t44 * t99;
t43 = t100 * t95 - t103 * t123;
t12 = -mrSges(7,2) * t43 + mrSges(7,3) * t26;
t27 = t102 * t44 + t99 * t120;
t13 = mrSges(7,1) * t43 - mrSges(7,3) * t27;
t147 = t102 * t12 - t99 * t13;
t55 = -mrSges(7,1) * t102 + mrSges(7,2) * t99;
t146 = m(7) * pkin(5) + mrSges(6,1) - t55;
t145 = t97 ^ 2;
t112 = -t97 * t104 - pkin(1);
t119 = qJ(3) * t101;
t28 = (t112 - t119) * t94;
t144 = -0.2e1 * t28;
t143 = t26 / 0.2e1;
t142 = t27 / 0.2e1;
t57 = Ifges(7,5) * t99 + Ifges(7,6) * t102;
t141 = t57 / 0.2e1;
t140 = -t99 / 0.2e1;
t139 = t99 / 0.2e1;
t138 = pkin(3) + pkin(8);
t137 = pkin(1) * t95;
t136 = t102 / 0.2e1;
t135 = Ifges(7,4) * t99;
t46 = pkin(8) * t120 + t101 * t137;
t34 = -t95 * qJ(3) - t46;
t29 = pkin(3) * t120 - t34;
t19 = pkin(4) * t120 - pkin(9) * t95 + t29;
t22 = (-t96 * t101 + t112) * t94;
t5 = -t100 * t22 + t103 * t19;
t3 = -pkin(5) * t120 - t5;
t134 = t103 * t3;
t73 = t104 * t137;
t45 = -pkin(8) * t123 + t73;
t133 = t45 * mrSges(3,1);
t132 = t46 * mrSges(3,2);
t11 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t31 = mrSges(6,1) * t120 - mrSges(6,3) * t44;
t130 = -t11 + t31;
t6 = t100 * t19 + t103 * t22;
t21 = t43 * mrSges(6,1) + t44 * mrSges(6,2);
t48 = mrSges(5,1) * t123 - t95 * mrSges(5,3);
t129 = t21 - t48;
t56 = t100 * mrSges(6,1) + t103 * mrSges(6,2);
t128 = t56 + mrSges(5,3);
t50 = mrSges(5,1) * t120 + t95 * mrSges(5,2);
t51 = mrSges(4,1) * t123 + t95 * mrSges(4,2);
t127 = t102 ^ 2 + t99 ^ 2;
t91 = t100 ^ 2;
t93 = t103 ^ 2;
t126 = t91 + t93;
t125 = Ifges(7,4) * t102;
t124 = t100 * t96;
t121 = t103 * t99;
t118 = t102 * t103;
t7 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t43;
t117 = Ifges(6,5) * t44 - Ifges(6,6) * t43 + Ifges(6,3) * t120;
t84 = t95 * pkin(2);
t116 = -t95 * qJ(4) - t73 - t84;
t115 = t126 * mrSges(6,3);
t17 = -(-pkin(4) - t138) * t123 + t116;
t10 = pkin(5) * t43 - pkin(10) * t44 - t17;
t4 = pkin(10) * t120 + t6;
t1 = t10 * t102 - t4 * t99;
t2 = t10 * t99 + t102 * t4;
t113 = -t1 * t99 + t102 * t2;
t111 = t6 * t100 + t5 * t103;
t110 = mrSges(7,1) * t99 + mrSges(7,2) * t102;
t54 = pkin(5) * t100 - pkin(10) * t103 + t97;
t32 = t102 * t54 - t99 * t124;
t33 = t102 * t124 + t54 * t99;
t109 = t102 * t33 - t32 * t99;
t52 = -mrSges(7,2) * t100 - mrSges(7,3) * t121;
t53 = mrSges(7,1) * t100 - mrSges(7,3) * t118;
t108 = t102 * t52 - t99 * t53;
t37 = Ifges(7,5) * t118 - Ifges(7,6) * t121 + Ifges(7,3) * t100;
t107 = Ifges(3,6) * t120 + t151 * t123 + t150 * t95;
t105 = qJ(3) ^ 2;
t89 = t96 ^ 2;
t88 = Ifges(6,5) * t103;
t76 = t93 * t96;
t75 = t93 * t89;
t61 = Ifges(6,1) * t103 - Ifges(6,4) * t100;
t60 = Ifges(7,1) * t99 + t125;
t59 = Ifges(6,4) * t103 - Ifges(6,2) * t100;
t58 = Ifges(7,2) * t102 + t135;
t49 = -mrSges(4,1) * t120 - mrSges(4,3) * t95;
t47 = t110 * t103;
t39 = Ifges(7,5) * t100 + (Ifges(7,1) * t102 - t135) * t103;
t38 = Ifges(7,6) * t100 + (-Ifges(7,2) * t99 + t125) * t103;
t36 = -t45 - t84;
t35 = (-pkin(2) * t104 - pkin(1) - t119) * t94;
t30 = -mrSges(6,2) * t120 - mrSges(6,3) * t43;
t23 = t138 * t123 + t116;
t15 = Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t120;
t14 = Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t120;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t43;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t43;
t16 = [0.2e1 * t1 * t13 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + t44 * t15 - 0.2e1 * t17 * t21 + 0.2e1 * t23 * t48 + t26 * t8 + t27 * t9 + 0.2e1 * t29 * t50 + 0.2e1 * t6 * t30 + 0.2e1 * t5 * t31 + 0.2e1 * t34 * t49 + 0.2e1 * t36 * t51 + Ifges(2,3) + (t7 - t14) * t43 + (t107 - 0.2e1 * t132 + 0.2e1 * t133) * t95 + m(3) * (t45 ^ 2 + t46 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (m(3) * pkin(1) ^ 2 * t94 + (mrSges(5,2) * t144 - 0.2e1 * t45 * mrSges(3,3) - 0.2e1 * t35 * mrSges(4,3) + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(5,3) + Ifges(4,2) + Ifges(3,1)) * t101) * t94 + (-(2 * Ifges(4,4)) + t151) * t95) * t101 + (0.2e1 * t35 * mrSges(4,2) + 0.2e1 * t46 * mrSges(3,3) + mrSges(5,3) * t144 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(5,2) + Ifges(4,3)) * t104) * t94 + 0.2e1 * (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) * t123 + (Ifges(3,6) - (2 * Ifges(5,4)) - (2 * Ifges(4,5))) * t95 + t117) * t104) * t94; (-t59 / 0.2e1 + t37 / 0.2e1) * t43 + t107 + m(6) * (t111 * t96 - t17 * t97) - t132 + t133 + (t9 * t136 + t8 * t140 - t5 * mrSges(6,3) + t15 / 0.2e1 + t130 * t96) * t103 + m(7) * (t1 * t32 - t96 * t134 + t33 * t2) + t129 * t97 - t17 * t56 + t44 * t61 / 0.2e1 + t3 * t47 - pkin(2) * t51 + t2 * t52 + t1 * t53 - t34 * mrSges(4,3) + t36 * mrSges(4,2) + t29 * mrSges(5,2) + t32 * t13 + t33 * t12 - t23 * mrSges(5,3) + t38 * t143 + t39 * t142 + m(5) * (qJ(3) * t29 - t23 * t97) + m(4) * (-pkin(2) * t36 - qJ(3) * t34) + (-t49 + t50) * qJ(3) + (t96 * t30 - t6 * mrSges(6,3) - t14 / 0.2e1 + t7 / 0.2e1) * t100 + (-Ifges(4,4) * t101 + (-Ifges(6,6) * t100 / 0.2e1 + t88 / 0.2e1 - Ifges(5,4) - Ifges(4,5)) * t104) * t94; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t32 * t53 + 0.2e1 * t33 * t52 + (t37 - t59) * t100 + (t102 * t39 + t149 * t47 - t38 * t99 + t61) * t103 + m(7) * (t32 ^ 2 + t33 ^ 2 + t75) + m(6) * (t89 * t91 + t145 + t75) + m(5) * (t105 + t145) + m(4) * (pkin(2) ^ 2 + t105) + t128 * t148 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * qJ(3) + t115 * t149 + t150; -t102 * t13 - t99 * t12 + m(7) * (-t1 * t102 - t2 * t99) + m(6) * t17 + m(5) * t23 + m(4) * t36 + t51 - t129; -m(4) * pkin(2) - t102 * t53 - t99 * t52 + mrSges(4,2) + m(7) * (-t102 * t32 - t33 * t99) + (-m(6) / 0.2e1 - m(5) / 0.2e1) * t148 - t128; m(7) * t127 + m(4) + m(5) + m(6); t130 * t103 + (t30 + t147) * t100 + m(7) * (t113 * t100 - t134) + m(6) * t111 + m(5) * t29 + t50; m(5) * qJ(3) - t103 * t47 + mrSges(5,2) + t108 * t100 - t115 + m(7) * (t109 * t100 + t76) + m(6) * (t91 * t96 + t76); 0; m(5) + m(6) * t126 + m(7) * (t127 * t91 + t93); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t113 * mrSges(7,3) + t8 * t136 + t9 * t139 + t43 * t141 + t60 * t142 + t58 * t143 + t3 * t55 + t117 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t113 + t147) * pkin(10); -pkin(5) * t47 + t39 * t139 + t38 * t136 + t88 + (m(7) * t109 + t108) * pkin(10) + (t60 * t136 + t58 * t140 + t146 * t96) * t103 + t109 * mrSges(7,3) + (-t96 * mrSges(6,2) - Ifges(6,6) + t141) * t100; 0; t146 * t103 + (-mrSges(6,2) + (m(7) * pkin(10) + mrSges(7,3)) * t127) * t100; Ifges(6,3) + t102 * t58 + m(7) * (t127 * pkin(10) ^ 2 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t55 + t99 * t60 + 0.2e1 * t127 * pkin(10) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t32 - mrSges(7,2) * t33 + t37; t55; -t110 * t100; -t110 * pkin(10) + t57; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
