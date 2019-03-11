% Calculate joint inertia matrix for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:05
% EndTime: 2019-03-08 23:11:07
% DurationCPUTime: 0.78s
% Computational Cost: add. (1029->214), mult. (2092->293), div. (0->0), fcn. (2115->10), ass. (0->89)
t126 = mrSges(6,2) - mrSges(5,1);
t72 = sin(qJ(6));
t76 = cos(qJ(6));
t51 = mrSges(7,1) * t72 + mrSges(7,2) * t76;
t125 = mrSges(6,3) + t51;
t110 = t72 ^ 2 + t76 ^ 2;
t100 = t110 * mrSges(7,3);
t70 = sin(pkin(6));
t75 = sin(qJ(2));
t114 = t70 * t75;
t71 = cos(pkin(6));
t74 = sin(qJ(3));
t78 = cos(qJ(3));
t39 = -t74 * t114 + t71 * t78;
t40 = t78 * t114 + t71 * t74;
t73 = sin(qJ(4));
t77 = cos(qJ(4));
t19 = t39 * t73 + t40 * t77;
t14 = t19 ^ 2;
t58 = pkin(3) * t73 + qJ(5);
t123 = t58 ^ 2;
t68 = t78 ^ 2;
t122 = pkin(4) + pkin(10);
t121 = -pkin(9) - pkin(8);
t120 = mrSges(7,1) * t76;
t119 = Ifges(7,4) * t72;
t118 = Ifges(7,4) * t76;
t117 = t19 * t58;
t47 = t73 * t74 - t77 * t78;
t116 = t47 * t72;
t115 = t47 * t76;
t79 = cos(qJ(2));
t113 = t70 * t79;
t112 = mrSges(6,1) + mrSges(5,3);
t111 = -mrSges(5,2) + mrSges(6,3);
t109 = t74 ^ 2 + t68;
t108 = qJ(5) * t19;
t107 = qJ(5) * t58;
t48 = t73 * t78 + t74 * t77;
t104 = Ifges(7,5) * t116 + Ifges(7,6) * t115 + Ifges(7,3) * t48;
t103 = t121 * t74;
t54 = t121 * t78;
t31 = -t77 * t103 - t54 * t73;
t33 = t73 * t103 - t77 * t54;
t102 = t31 ^ 2 + t33 ^ 2;
t62 = -pkin(3) * t78 - pkin(2);
t61 = -pkin(3) * t77 - pkin(4);
t101 = m(7) * t110;
t99 = t110 * t122;
t63 = Ifges(7,5) * t76;
t98 = -Ifges(7,6) * t72 + t63;
t97 = 0.2e1 * t125;
t87 = -qJ(5) * t48 + t62;
t11 = t122 * t47 + t87;
t15 = pkin(5) * t48 + t31;
t1 = -t11 * t72 + t15 * t76;
t2 = t11 * t76 + t15 * t72;
t96 = t1 * t76 + t2 * t72;
t17 = -t77 * t39 + t40 * t73;
t5 = t72 * t113 + t17 * t76;
t6 = -t76 * t113 + t17 * t72;
t95 = t5 * t76 + t6 * t72;
t94 = mrSges(6,2) - t100;
t93 = -t72 * mrSges(7,2) + t120;
t92 = t17 * t31 + t19 * t33;
t23 = mrSges(7,1) * t48 - mrSges(7,3) * t116;
t24 = -mrSges(7,2) * t48 + mrSges(7,3) * t115;
t91 = t76 * t23 + t72 * t24;
t90 = -t39 * t74 + t40 * t78;
t89 = -0.2e1 * t100;
t52 = -Ifges(7,2) * t72 + t118;
t53 = Ifges(7,1) * t76 - t119;
t88 = -t72 * t52 + t76 * t53 + Ifges(6,1) + Ifges(5,3);
t86 = (mrSges(5,1) * t77 - mrSges(5,2) * t73) * pkin(3);
t85 = -t95 * mrSges(7,3) + (t111 + t51) * t19 + t126 * t17;
t10 = Ifges(7,5) * t48 + (Ifges(7,1) * t72 + t118) * t47;
t16 = -t47 * pkin(5) + t33;
t9 = Ifges(7,6) * t48 + (Ifges(7,2) * t76 + t119) * t47;
t84 = -t96 * mrSges(7,3) + t111 * t33 + t53 * t116 / 0.2e1 - t72 * t9 / 0.2e1 + t52 * t115 / 0.2e1 + t76 * t10 / 0.2e1 + t16 * t51 + (-Ifges(5,6) + Ifges(6,5)) * t47 + t126 * t31 + (t98 / 0.2e1 + Ifges(5,5) - Ifges(6,4)) * t48;
t81 = qJ(5) ^ 2;
t64 = t70 ^ 2;
t57 = t64 * t79 ^ 2;
t56 = -pkin(10) + t61;
t50 = -mrSges(4,1) * t78 + mrSges(4,2) * t74;
t26 = -mrSges(6,2) * t47 - mrSges(6,3) * t48;
t25 = mrSges(5,1) * t47 + mrSges(5,2) * t48;
t22 = pkin(4) * t47 + t87;
t21 = t93 * t47;
t3 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t14) + m(4) * (t39 ^ 2 + t40 ^ 2 + t57) + m(3) * (t64 * t75 ^ 2 + t71 ^ 2 + t57) + (m(6) + m(5)) * (t17 ^ 2 + t14 + t57); t5 * t23 + t6 * t24 + t112 * t48 * t17 + t90 * mrSges(4,3) + (-t112 * t47 - t21) * t19 + (-t75 * mrSges(3,2) + (mrSges(3,1) - t25 - t26 - t50) * t79) * t70 + m(6) * (-t22 * t113 + t92) + m(7) * (t1 * t5 + t16 * t19 + t2 * t6) + m(5) * (-t62 * t113 + t92) + m(4) * (pkin(2) * t113 + t90 * pkin(8)); Ifges(4,2) * t68 - 0.2e1 * pkin(2) * t50 + 0.2e1 * t1 * t23 - 0.2e1 * t16 * t21 + 0.2e1 * t2 * t24 + 0.2e1 * t22 * t26 + 0.2e1 * t62 * t25 + Ifges(3,3) + (Ifges(4,1) * t74 + 0.2e1 * Ifges(4,4) * t78) * t74 + 0.2e1 * t109 * pkin(8) * mrSges(4,3) + m(4) * (t109 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t62 ^ 2 + t102) + m(6) * (t22 ^ 2 + t102) + m(7) * (t1 ^ 2 + t16 ^ 2 + t2 ^ 2) + ((Ifges(6,2) + Ifges(5,1)) * t48 + 0.2e1 * t112 * t31 + t104) * t48 + (t72 * t10 + t76 * t9 + (Ifges(5,2) + Ifges(6,3)) * t47 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t48 - 0.2e1 * t112 * t33) * t47; t39 * mrSges(4,1) - t40 * mrSges(4,2) + m(6) * (t17 * t61 + t117) + m(7) * (t95 * t56 + t117) + m(5) * (-t17 * t77 + t19 * t73) * pkin(3) + t85; m(6) * (t31 * t61 + t33 * t58) + (-t74 * mrSges(4,1) - t78 * mrSges(4,2)) * pkin(8) + (-t58 * t47 + t61 * t48) * mrSges(6,1) + Ifges(4,6) * t78 + Ifges(4,5) * t74 - t58 * t21 + t84 + t91 * t56 + m(7) * (t16 * t58 + t96 * t56) + (m(5) * (-t31 * t77 + t33 * t73) + (-t73 * t47 - t77 * t48) * mrSges(5,3)) * pkin(3); 0.2e1 * t61 * mrSges(6,2) + Ifges(4,3) + t58 * t97 + 0.2e1 * t86 + t56 * t89 + m(7) * (t110 * t56 ^ 2 + t123) + m(6) * (t61 ^ 2 + t123) + m(5) * (t73 ^ 2 + t77 ^ 2) * pkin(3) ^ 2 + t88; m(6) * (-pkin(4) * t17 + t108) + m(7) * (-t122 * t95 + t108) + t85; -t91 * t122 + m(6) * (-pkin(4) * t31 + qJ(5) * t33) + (-pkin(4) * t48 - qJ(5) * t47) * mrSges(6,1) - qJ(5) * t21 + t84 + m(7) * (qJ(5) * t16 - t122 * t96); t86 + (-pkin(4) + t61) * mrSges(6,2) + m(7) * (-t56 * t99 + t107) + m(6) * (-pkin(4) * t61 + t107) + t88 + (-t56 + t122) * t100 + t125 * (qJ(5) + t58); -0.2e1 * pkin(4) * mrSges(6,2) + qJ(5) * t97 - t122 * t89 + m(7) * (t110 * t122 ^ 2 + t81) + m(6) * (pkin(4) ^ 2 + t81) + t88; m(6) * t17 + m(7) * t95; m(6) * t31 + m(7) * t96 + t48 * mrSges(6,1) + t91; m(6) * t61 + t56 * t101 + t94; -m(6) * pkin(4) - m(7) * t99 + t94; m(6) + t101; mrSges(7,1) * t5 - mrSges(7,2) * t6; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t104; t93 * t56 + t98; -t122 * t120 + t63 + (mrSges(7,2) * t122 - Ifges(7,6)) * t72; t93; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
