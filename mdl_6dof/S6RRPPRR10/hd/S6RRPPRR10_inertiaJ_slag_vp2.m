% Calculate joint inertia matrix for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:12
% EndTime: 2019-03-09 09:34:15
% DurationCPUTime: 1.07s
% Computational Cost: add. (1740->262), mult. (3150->362), div. (0->0), fcn. (3214->8), ass. (0->102)
t151 = pkin(3) + pkin(7);
t111 = sin(qJ(2));
t113 = cos(qJ(2));
t150 = t111 ^ 2 + t113 ^ 2;
t109 = sin(qJ(6));
t112 = cos(qJ(6));
t106 = sin(pkin(10));
t107 = cos(pkin(10));
t110 = sin(qJ(5));
t138 = cos(qJ(5));
t72 = -t138 * t106 - t110 * t107;
t74 = -t110 * t106 + t138 * t107;
t125 = t109 * t72 + t112 * t74;
t41 = t109 * t74 - t112 * t72;
t149 = t125 ^ 2 + t41 ^ 2;
t108 = -pkin(2) - qJ(4);
t137 = -pkin(8) + t108;
t79 = t137 * t106;
t80 = t137 * t107;
t45 = -t110 * t79 + t138 * t80;
t24 = -pkin(9) * t74 + t45;
t46 = t110 * t80 + t138 * t79;
t25 = pkin(9) * t72 + t46;
t7 = -t109 * t25 + t112 * t24;
t8 = t109 * t24 + t112 * t25;
t148 = t125 * t7 + t41 * t8;
t126 = -qJ(3) * t111 - pkin(1);
t70 = t108 * t113 + t126;
t85 = t151 * t111;
t76 = t107 * t85;
t30 = pkin(4) * t111 + t76 + (pkin(8) * t113 - t70) * t106;
t132 = t107 * t113;
t44 = t106 * t85 + t107 * t70;
t33 = -pkin(8) * t132 + t44;
t12 = -t110 * t33 + t138 * t30;
t59 = t72 * t113;
t4 = pkin(5) * t111 - pkin(9) * t59 + t12;
t13 = t110 * t30 + t138 * t33;
t58 = t74 * t113;
t5 = -pkin(9) * t58 + t13;
t2 = -t109 * t5 + t112 * t4;
t3 = t109 * t4 + t112 * t5;
t147 = t125 * t2 + t3 * t41;
t146 = t109 * t41 + t112 * t125;
t141 = Ifges(6,5) * t59 - Ifges(6,6) * t58 + Ifges(6,3) * t111;
t140 = -m(4) * pkin(2) + mrSges(4,2);
t139 = -t106 / 0.2e1;
t82 = t106 * mrSges(5,1) + t107 * mrSges(5,2);
t136 = t150 * pkin(7) ^ 2;
t135 = Ifges(5,4) * t106;
t134 = Ifges(5,4) * t107;
t89 = t106 * pkin(4) + qJ(3);
t86 = t151 * t113;
t133 = t106 * t113;
t131 = t106 ^ 2 + t107 ^ 2;
t130 = t72 ^ 2 + t74 ^ 2;
t28 = -t109 * t59 - t112 * t58;
t29 = -t109 * t58 + t112 * t59;
t129 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t111;
t62 = pkin(4) * t132 + t86;
t32 = t58 * mrSges(6,1) + t59 * mrSges(6,2);
t47 = -t72 * mrSges(6,1) + t74 * mrSges(6,2);
t11 = -t28 * mrSges(7,1) + t29 * mrSges(7,2);
t14 = mrSges(7,1) * t41 + mrSges(7,2) * t125;
t128 = mrSges(7,1) * t125 - t41 * mrSges(7,2);
t127 = m(5) * t131;
t124 = t131 * mrSges(5,3);
t63 = mrSges(5,1) * t132 - mrSges(5,2) * t133;
t122 = t12 * t74 - t13 * t72;
t121 = t45 * t74 - t46 * t72;
t37 = Ifges(7,6) * t41;
t38 = Ifges(7,5) * t125;
t120 = t7 * mrSges(7,1) - t8 * mrSges(7,2) - t37 + t38;
t43 = -t106 * t70 + t76;
t119 = t44 * t106 + t43 * t107;
t118 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t129;
t117 = (mrSges(7,1) * t112 - mrSges(7,2) * t109) * pkin(5);
t114 = qJ(3) ^ 2;
t84 = Ifges(5,1) * t107 - t135;
t83 = -Ifges(5,2) * t106 + t134;
t81 = -pkin(2) * t113 + t126;
t78 = -mrSges(5,2) * t111 - mrSges(5,3) * t132;
t77 = mrSges(5,1) * t111 + mrSges(5,3) * t133;
t69 = Ifges(6,5) * t74;
t68 = Ifges(6,6) * t72;
t57 = Ifges(5,5) * t111 + (-Ifges(5,1) * t106 - t134) * t113;
t56 = Ifges(5,6) * t111 + (-Ifges(5,2) * t107 - t135) * t113;
t52 = -pkin(5) * t72 + t89;
t51 = mrSges(6,1) * t111 - mrSges(6,3) * t59;
t50 = -mrSges(6,2) * t111 - mrSges(6,3) * t58;
t49 = Ifges(6,1) * t74 + Ifges(6,4) * t72;
t48 = Ifges(6,4) * t74 + Ifges(6,2) * t72;
t34 = pkin(5) * t58 + t62;
t27 = Ifges(6,1) * t59 - Ifges(6,4) * t58 + Ifges(6,5) * t111;
t26 = Ifges(6,4) * t59 - Ifges(6,2) * t58 + Ifges(6,6) * t111;
t18 = mrSges(7,1) * t111 - mrSges(7,3) * t29;
t17 = -mrSges(7,2) * t111 + mrSges(7,3) * t28;
t16 = Ifges(7,1) * t125 - Ifges(7,4) * t41;
t15 = Ifges(7,4) * t125 - Ifges(7,2) * t41;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t111;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t111;
t1 = [t29 * t10 + 0.2e1 * t34 * t11 + 0.2e1 * t12 * t51 + 0.2e1 * t13 * t50 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t58 * t26 + t59 * t27 + t28 * t9 + 0.2e1 * t62 * t32 + 0.2e1 * t43 * t77 + 0.2e1 * t44 * t78 + 0.2e1 * t86 * t63 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t81 * mrSges(4,2) - t106 * t57 - t107 * t56 + (Ifges(3,2) + Ifges(4,3)) * t113) * t113 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t81 * mrSges(4,3) + (Ifges(5,3) + Ifges(3,1) + Ifges(4,2)) * t111 + (-Ifges(5,5) * t106 - Ifges(5,6) * t107 + (2 * Ifges(3,4)) + (2 * Ifges(4,6))) * t113 + t129 + t141) * t111 + m(4) * (t81 ^ 2 + t136) + m(3) * (pkin(1) ^ 2 + t136) + m(5) * (t43 ^ 2 + t44 ^ 2 + t86 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2 + t62 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t150; -t147 * mrSges(7,3) + (-t43 * mrSges(5,3) + t108 * t77 + t57 / 0.2e1) * t107 + (t108 * t78 - t44 * mrSges(5,3) - t56 / 0.2e1) * t106 + t125 * t10 / 0.2e1 - t41 * t9 / 0.2e1 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t107 / 0.2e1 + Ifges(5,6) * t139 + t69 / 0.2e1 + t68 / 0.2e1 + t38 / 0.2e1 - t37 / 0.2e1 - Ifges(4,4) + Ifges(3,5)) * t111 + t86 * t82 + t89 * t32 + t72 * t26 / 0.2e1 + t74 * t27 / 0.2e1 + t59 * t49 / 0.2e1 + t62 * t47 + qJ(3) * t63 + t46 * t50 + t45 * t51 + t52 * t11 - t58 * t48 / 0.2e1 + t34 * t14 + t28 * t15 / 0.2e1 + t29 * t16 / 0.2e1 + t8 * t17 + t7 * t18 + (-t107 * t83 / 0.2e1 + t84 * t139 + qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6)) * t113 + m(5) * (qJ(3) * t86 + t119 * t108) - t122 * mrSges(6,3) + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t113 + (-mrSges(3,1) + t140) * t111) * pkin(7) + m(6) * (t12 * t45 + t13 * t46 + t62 * t89) + m(7) * (t2 * t7 + t3 * t8 + t34 * t52); -0.2e1 * pkin(2) * mrSges(4,2) - t106 * t83 + t107 * t84 + 0.2e1 * t52 * t14 - t41 * t15 + t125 * t16 + 0.2e1 * t89 * t47 + t72 * t48 + t74 * t49 + Ifges(4,1) + Ifges(3,3) + m(7) * (t52 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(6) * (t45 ^ 2 + t46 ^ 2 + t89 ^ 2) + m(5) * (t131 * t108 ^ 2 + t114) + m(4) * (pkin(2) ^ 2 + t114) + 0.2e1 * (mrSges(4,3) + t82) * qJ(3) - 0.2e1 * t148 * mrSges(7,3) - 0.2e1 * t121 * mrSges(6,3) - 0.2e1 * t108 * t124; t106 * t78 + t107 * t77 + t41 * t17 + t125 * t18 - t72 * t50 + t74 * t51 + (m(4) * pkin(7) + mrSges(4,1)) * t111 + m(7) * t147 + m(6) * t122 + m(5) * t119; m(6) * t121 + m(7) * t148 - t130 * mrSges(6,3) - t149 * mrSges(7,3) + t108 * t127 - t124 + t140; m(6) * t130 + m(7) * t149 + m(4) + t127; m(5) * t86 + m(6) * t62 + m(7) * t34 + t11 + t32 + t63; m(5) * qJ(3) + m(6) * t89 + m(7) * t52 + t14 + t47 + t82; 0; m(5) + m(6) + m(7); t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t109 * t3 + t112 * t2) + t109 * t17 + t112 * t18) * pkin(5) + t118 + t141; t45 * mrSges(6,1) - t46 * mrSges(6,2) + t68 + t69 + (m(7) * (t109 * t8 + t112 * t7) - t146 * mrSges(7,3)) * pkin(5) + t120; m(7) * t146 * pkin(5) + t74 * mrSges(6,1) + t72 * mrSges(6,2) + t128; 0; Ifges(6,3) + Ifges(7,3) + m(7) * (t109 ^ 2 + t112 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t117; t118; t120; t128; 0; Ifges(7,3) + t117; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
