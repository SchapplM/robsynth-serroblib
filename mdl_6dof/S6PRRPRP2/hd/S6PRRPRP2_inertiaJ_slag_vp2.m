% Calculate joint inertia matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:31
% EndTime: 2019-03-08 21:28:33
% DurationCPUTime: 0.82s
% Computational Cost: add. (1035->232), mult. (2185->319), div. (0->0), fcn. (2287->10), ass. (0->94)
t127 = Ifges(7,2) + Ifges(6,3);
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t126 = t76 ^ 2 + t79 ^ 2;
t73 = sin(pkin(6));
t81 = cos(qJ(2));
t103 = t73 * t81;
t78 = sin(qJ(2));
t104 = t73 * t78;
t75 = cos(pkin(6));
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t38 = -t77 * t104 + t75 * t80;
t39 = t80 * t104 + t75 * t77;
t72 = sin(pkin(11));
t74 = cos(pkin(11));
t19 = t38 * t72 + t39 * t74;
t11 = -t76 * t103 + t19 * t79;
t9 = t79 * t103 + t19 * t76;
t125 = t11 * t79 + t76 * t9;
t124 = 0.2e1 * t126;
t118 = m(5) * pkin(3);
t123 = m(7) + m(6);
t122 = mrSges(6,3) + mrSges(7,2);
t121 = -m(7) * pkin(5) - mrSges(7,1);
t17 = -t74 * t38 + t39 * t72;
t16 = t17 ^ 2;
t102 = -qJ(4) - pkin(8);
t50 = t102 * t80;
t93 = t102 * t77;
t30 = -t50 * t72 - t74 * t93;
t120 = t30 ^ 2;
t71 = t80 ^ 2;
t119 = 0.2e1 * t30;
t117 = pkin(3) * t72;
t116 = pkin(3) * t74;
t114 = Ifges(6,4) * t76;
t113 = Ifges(6,4) * t79;
t112 = Ifges(7,5) * t76;
t111 = Ifges(7,5) * t79;
t44 = t72 * t77 - t74 * t80;
t110 = Ifges(6,6) * t44;
t109 = Ifges(7,6) * t44;
t107 = t17 * t30;
t45 = t72 * t80 + t74 * t77;
t106 = t45 * t76;
t105 = t45 * t79;
t24 = -mrSges(6,2) * t44 - mrSges(6,3) * t106;
t27 = -mrSges(7,2) * t106 + mrSges(7,3) * t44;
t101 = t24 + t27;
t25 = mrSges(6,1) * t44 - mrSges(6,3) * t105;
t26 = -t44 * mrSges(7,1) + mrSges(7,2) * t105;
t100 = t25 - t26;
t63 = -pkin(3) * t80 - pkin(2);
t23 = pkin(4) * t44 - pkin(9) * t45 + t63;
t32 = -t74 * t50 + t72 * t93;
t4 = t76 * t23 + t79 * t32;
t48 = -t79 * mrSges(6,1) + t76 * mrSges(6,2);
t99 = t48 - mrSges(5,1);
t61 = pkin(9) + t117;
t98 = t126 * t61 ^ 2;
t96 = t77 ^ 2 + t71;
t62 = -pkin(4) - t116;
t28 = t44 * mrSges(5,1) + t45 * mrSges(5,2);
t92 = t125 * t61;
t90 = Ifges(7,6) * t106 + t127 * t44 + (Ifges(7,4) + Ifges(6,5)) * t105;
t88 = t76 * mrSges(6,1) + t79 * mrSges(6,2);
t47 = -t79 * mrSges(7,1) - t76 * mrSges(7,3);
t87 = t76 * mrSges(7,1) - t79 * mrSges(7,3);
t86 = pkin(5) * t79 + qJ(6) * t76;
t85 = pkin(5) * t76 - qJ(6) * t79;
t3 = t23 * t79 - t32 * t76;
t84 = -t38 * t77 + t39 * t80;
t67 = t73 ^ 2;
t66 = Ifges(7,4) * t76;
t65 = Ifges(6,5) * t76;
t64 = Ifges(6,6) * t79;
t59 = t67 * t81 ^ 2;
t54 = Ifges(6,1) * t76 + t113;
t53 = Ifges(7,1) * t76 - t111;
t52 = Ifges(6,2) * t79 + t114;
t51 = -Ifges(7,3) * t79 + t112;
t49 = -mrSges(4,1) * t80 + mrSges(4,2) * t77;
t43 = t62 - t86;
t22 = t88 * t45;
t21 = t87 * t45;
t15 = Ifges(6,5) * t44 + (Ifges(6,1) * t79 - t114) * t45;
t14 = Ifges(7,4) * t44 + (Ifges(7,1) * t79 + t112) * t45;
t13 = t110 + (-Ifges(6,2) * t76 + t113) * t45;
t12 = t109 + (Ifges(7,3) * t76 + t111) * t45;
t5 = t85 * t45 + t30;
t2 = -pkin(5) * t44 - t3;
t1 = qJ(6) * t44 + t4;
t6 = [m(2) + m(5) * (t19 ^ 2 + t16 + t59) + m(4) * (t38 ^ 2 + t39 ^ 2 + t59) + m(3) * (t67 * t78 ^ 2 + t75 ^ 2 + t59) + t123 * (t11 ^ 2 + t9 ^ 2 + t16); -t19 * t44 * mrSges(5,3) - t100 * t9 + t101 * t11 + t84 * mrSges(4,3) + (t45 * mrSges(5,3) + t21 + t22) * t17 + (-t78 * mrSges(3,2) + (mrSges(3,1) - t28 - t49) * t81) * t73 + m(7) * (t1 * t11 + t17 * t5 + t2 * t9) + m(6) * (t11 * t4 - t3 * t9 + t107) + m(5) * (-t63 * t103 + t19 * t32 + t107) + m(4) * (pkin(2) * t103 + t84 * pkin(8)); Ifges(4,2) * t71 - 0.2e1 * pkin(2) * t49 + 0.2e1 * t1 * t27 + 0.2e1 * t2 * t26 + 0.2e1 * t5 * t21 + t22 * t119 + 0.2e1 * t4 * t24 + 0.2e1 * t3 * t25 + 0.2e1 * t63 * t28 + Ifges(3,3) + (Ifges(4,1) * t77 + 0.2e1 * Ifges(4,4) * t80) * t77 + 0.2e1 * t96 * pkin(8) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t32 + Ifges(5,2) * t44 + t90) * t44 + m(4) * (t96 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t32 ^ 2 + t63 ^ 2 + t120) + m(6) * (t3 ^ 2 + t4 ^ 2 + t120) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (mrSges(5,3) * t119 + Ifges(5,1) * t45 - 0.2e1 * Ifges(5,4) * t44 + (t14 + t15) * t79 + (t12 - t13 - t110) * t76) * t45; t38 * mrSges(4,1) - t39 * mrSges(4,2) - t19 * mrSges(5,2) + (t47 + t99) * t17 + m(7) * (t17 * t43 + t92) + m(6) * (t17 * t62 + t92) + (-t17 * t74 + t19 * t72) * t118 + t122 * t125; -t32 * mrSges(5,2) + Ifges(4,5) * t77 + Ifges(4,6) * t80 + t43 * t21 + t62 * t22 + t5 * t47 + t99 * t30 + (-t77 * mrSges(4,1) - t80 * mrSges(4,2)) * pkin(8) + (t66 / 0.2e1 + t65 / 0.2e1 + t64 / 0.2e1 - Ifges(5,6) - mrSges(5,3) * t117) * t44 + (-t109 / 0.2e1 - t12 / 0.2e1 + t13 / 0.2e1 + t1 * mrSges(7,2) + t4 * mrSges(6,3) + t101 * t61) * t79 + (t14 / 0.2e1 + t15 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3) - t100 * t61) * t76 + m(7) * (t43 * t5 + (t1 * t79 + t2 * t76) * t61) + m(6) * (t30 * t62 + (-t3 * t76 + t4 * t79) * t61) + (-t30 * t74 + t32 * t72) * t118 + (-mrSges(5,3) * t116 + Ifges(5,5) + (t53 / 0.2e1 + t54 / 0.2e1) * t79 + (t51 / 0.2e1 - t52 / 0.2e1) * t76) * t45; 0.2e1 * t43 * t47 + 0.2e1 * t62 * t48 + Ifges(4,3) + Ifges(5,3) + (t52 - t51) * t79 + (t53 + t54) * t76 + m(7) * (t43 ^ 2 + t98) + m(6) * (t62 ^ 2 + t98) + t122 * t61 * t124 + (0.2e1 * mrSges(5,1) * t74 - 0.2e1 * mrSges(5,2) * t72 + (t72 ^ 2 + t74 ^ 2) * t118) * pkin(3); -m(5) * t103 + t123 * (t76 * t11 - t79 * t9); t100 * t79 + t101 * t76 + m(7) * (t1 * t76 - t2 * t79) + m(6) * (t3 * t79 + t4 * t76) + m(5) * t63 + t28; 0; m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t124; (-mrSges(6,1) + t121) * t9 + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t11; -Ifges(6,6) * t106 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) + qJ(6) * t27 + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) - pkin(5) * t26 + t90; -Ifges(7,6) * t79 + t64 + t65 + t66 - t85 * mrSges(7,2) + (-m(7) * t85 - t87 - t88) * t61; m(7) * t86 - t47 - t48; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t127; m(7) * t9; m(7) * t2 + t26; (m(7) * t61 + mrSges(7,2)) * t76; -m(7) * t79; t121; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
