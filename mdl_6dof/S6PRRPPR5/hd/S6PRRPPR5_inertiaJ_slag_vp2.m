% Calculate joint inertia matrix for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:12
% EndTime: 2019-03-08 21:18:14
% DurationCPUTime: 0.78s
% Computational Cost: add. (834->226), mult. (1714->311), div. (0->0), fcn. (1635->10), ass. (0->92)
t122 = pkin(4) + pkin(8);
t83 = sin(qJ(3));
t86 = cos(qJ(3));
t121 = t83 ^ 2 + t86 ^ 2;
t79 = cos(pkin(11));
t109 = t79 * t86;
t77 = sin(pkin(11));
t112 = t77 * t86;
t35 = mrSges(6,1) * t109 - mrSges(6,2) * t112;
t82 = sin(qJ(6));
t85 = cos(qJ(6));
t91 = t77 * t82 - t79 * t85;
t26 = t91 * t86;
t41 = -t85 * t77 - t82 * t79;
t27 = t41 * t86;
t8 = -t26 * mrSges(7,1) + t27 * mrSges(7,2);
t120 = t35 + t8;
t108 = mrSges(5,1) + mrSges(4,3);
t118 = -m(5) * pkin(3) + mrSges(5,2);
t78 = sin(pkin(6));
t84 = sin(qJ(2));
t111 = t78 * t84;
t80 = cos(pkin(6));
t33 = t86 * t111 + t80 * t83;
t29 = t33 ^ 2;
t117 = -t77 / 0.2e1;
t81 = -pkin(3) - qJ(5);
t116 = -pkin(9) + t81;
t115 = Ifges(6,4) * t77;
t114 = Ifges(6,4) * t79;
t31 = t83 * t111 - t80 * t86;
t113 = t31 * t83;
t87 = cos(qJ(2));
t110 = t78 * t87;
t51 = t77 * mrSges(6,1) + t79 * mrSges(6,2);
t107 = mrSges(5,3) + t51;
t98 = -t83 * qJ(4) - pkin(2);
t39 = t81 * t86 + t98;
t56 = t122 * t83;
t11 = t79 * t39 + t77 * t56;
t106 = t121 * pkin(8) ^ 2;
t57 = t122 * t86;
t105 = t77 ^ 2 + t79 ^ 2;
t104 = t33 * qJ(4);
t102 = t41 ^ 2 + t91 ^ 2;
t101 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t83;
t100 = m(6) * t105;
t99 = t105 * mrSges(6,3);
t14 = -t41 * mrSges(7,1) - mrSges(7,2) * t91;
t45 = t79 * t56;
t7 = t83 * pkin(5) + t45 + (pkin(9) * t86 - t39) * t77;
t9 = -pkin(9) * t109 + t11;
t1 = t7 * t85 - t82 * t9;
t2 = t7 * t82 + t85 * t9;
t96 = -t1 * t91 - t2 * t41;
t17 = t77 * t110 + t31 * t79;
t18 = -t79 * t110 + t31 * t77;
t3 = t17 * t85 - t18 * t82;
t4 = t17 * t82 + t18 * t85;
t95 = -t3 * t91 - t4 * t41;
t10 = -t39 * t77 + t45;
t94 = t10 * t79 + t11 * t77;
t48 = t116 * t77;
t49 = t116 * t79;
t12 = -t48 * t82 + t49 * t85;
t13 = t48 * t85 + t49 * t82;
t93 = -t12 * t91 - t13 * t41;
t92 = t17 * t79 + t18 * t77;
t90 = (t33 * t86 + t113) * pkin(8);
t88 = qJ(4) ^ 2;
t72 = t78 ^ 2;
t61 = t72 * t87 ^ 2;
t59 = pkin(5) * t77 + qJ(4);
t55 = -mrSges(4,1) * t86 + mrSges(4,2) * t83;
t54 = mrSges(5,2) * t86 - mrSges(5,3) * t83;
t53 = Ifges(6,1) * t79 - t115;
t52 = -Ifges(6,2) * t77 + t114;
t50 = -pkin(3) * t86 + t98;
t47 = -mrSges(6,2) * t83 - mrSges(6,3) * t109;
t46 = mrSges(6,1) * t83 + mrSges(6,3) * t112;
t38 = Ifges(7,5) * t91;
t37 = Ifges(7,6) * t41;
t30 = pkin(5) * t109 + t57;
t25 = Ifges(6,5) * t83 + (-Ifges(6,1) * t77 - t114) * t86;
t24 = Ifges(6,6) * t83 + (-Ifges(6,2) * t79 - t115) * t86;
t20 = mrSges(7,1) * t83 - mrSges(7,3) * t27;
t19 = -mrSges(7,2) * t83 + mrSges(7,3) * t26;
t16 = -Ifges(7,1) * t91 + Ifges(7,4) * t41;
t15 = -Ifges(7,4) * t91 + Ifges(7,2) * t41;
t6 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t83;
t5 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t83;
t21 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t29) + m(6) * (t17 ^ 2 + t18 ^ 2 + t29) + m(3) * (t72 * t84 ^ 2 + t80 ^ 2 + t61) + (m(5) + m(4)) * (t31 ^ 2 + t29 + t61); t17 * t46 + t18 * t47 + t4 * t19 + t3 * t20 + t108 * t113 + (-t84 * mrSges(3,2) + (mrSges(3,1) - t54 - t55) * t87) * t78 + (t108 * t86 + t120) * t33 + m(7) * (t1 * t3 + t2 * t4 + t30 * t33) + m(5) * (-t50 * t110 + t90) + m(6) * (t10 * t17 + t11 * t18 + t33 * t57) + m(4) * (pkin(2) * t110 + t90); -0.2e1 * pkin(2) * t55 + 0.2e1 * t1 * t20 + 0.2e1 * t10 * t46 + 0.2e1 * t11 * t47 + 0.2e1 * t2 * t19 + t26 * t5 + t27 * t6 + 0.2e1 * t30 * t8 + 0.2e1 * t57 * t35 + 0.2e1 * t50 * t54 + Ifges(3,3) + (-t79 * t24 - t77 * t25 + (Ifges(5,3) + Ifges(4,2)) * t86) * t86 + m(5) * (t50 ^ 2 + t106) + m(4) * (pkin(2) ^ 2 + t106) + m(6) * (t10 ^ 2 + t11 ^ 2 + t57 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t30 ^ 2) + ((Ifges(6,3) + Ifges(4,1) + Ifges(5,2)) * t83 + (-Ifges(6,5) * t77 - Ifges(6,6) * t79 + (2 * Ifges(4,4)) + (2 * Ifges(5,6))) * t86 + t101) * t83 + 0.2e1 * t108 * pkin(8) * t121; (-mrSges(4,1) + mrSges(5,2)) * t31 - t95 * mrSges(7,3) - t92 * mrSges(6,3) + (-mrSges(4,2) + t14 + t107) * t33 + m(7) * (t12 * t3 + t13 * t4 + t33 * t59) + m(5) * (-pkin(3) * t31 + t104) + m(6) * (t92 * t81 + t104); t57 * t51 + t59 * t8 + t41 * t5 / 0.2e1 - t91 * t6 / 0.2e1 + t26 * t15 / 0.2e1 + t27 * t16 / 0.2e1 + t30 * t14 + qJ(4) * t35 + t13 * t19 + t12 * t20 - t96 * mrSges(7,3) + (t25 / 0.2e1 + t81 * t46 - t10 * mrSges(6,3)) * t79 + (-t24 / 0.2e1 + t81 * t47 - t11 * mrSges(6,3)) * t77 + m(6) * (qJ(4) * t57 + t94 * t81) + m(7) * (t1 * t12 + t13 * t2 + t30 * t59) + (Ifges(6,5) * t79 / 0.2e1 + Ifges(6,6) * t117 - t38 / 0.2e1 + t37 / 0.2e1 + Ifges(4,5) - Ifges(5,4) - pkin(3) * mrSges(5,1)) * t83 + (Ifges(4,6) - Ifges(5,5) + t53 * t117 - t79 * t52 / 0.2e1 + qJ(4) * mrSges(5,1)) * t86 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t86 + (-mrSges(4,1) + t118) * t83) * pkin(8); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t59 * t14 + t41 * t15 - t91 * t16 - t77 * t52 + t79 * t53 + Ifges(5,1) + Ifges(4,3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t59 ^ 2) + m(6) * (t105 * t81 ^ 2 + t88) + m(5) * (pkin(3) ^ 2 + t88) - 0.2e1 * t93 * mrSges(7,3) + 0.2e1 * t107 * qJ(4) - 0.2e1 * t81 * t99; m(5) * t31 + m(6) * t92 + m(7) * t95; -t41 * t19 - t91 * t20 + t79 * t46 + t77 * t47 + (m(5) * pkin(8) + mrSges(5,1)) * t83 + m(7) * t96 + m(6) * t94; m(7) * t93 - t102 * mrSges(7,3) + t81 * t100 + t118 - t99; m(7) * t102 + m(5) + t100; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t33; m(6) * t57 + m(7) * t30 + t120; m(6) * qJ(4) + m(7) * t59 + t14 + t51; 0; m(6) + m(7); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t101; mrSges(7,1) * t12 - mrSges(7,2) * t13 + t37 - t38; -mrSges(7,1) * t91 + mrSges(7,2) * t41; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
