% Calculate joint inertia matrix for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:18
% EndTime: 2019-03-08 20:09:20
% DurationCPUTime: 0.73s
% Computational Cost: add. (970->215), mult. (2069->294), div. (0->0), fcn. (2186->10), ass. (0->91)
t124 = Ifges(7,2) + Ifges(6,3);
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t123 = t74 ^ 2 + t77 ^ 2;
t71 = sin(pkin(6));
t78 = cos(qJ(2));
t100 = t71 * t78;
t112 = cos(qJ(4));
t76 = sin(qJ(2));
t101 = t71 * t76;
t70 = sin(pkin(11));
t72 = cos(pkin(11));
t73 = cos(pkin(6));
t38 = -t101 * t70 + t72 * t73;
t39 = t101 * t72 + t70 * t73;
t75 = sin(qJ(4));
t19 = t112 * t39 + t38 * t75;
t11 = -t100 * t74 + t19 * t77;
t9 = t100 * t77 + t19 * t74;
t122 = t11 * t77 + t74 * t9;
t121 = 0.2e1 * t123;
t120 = m(6) + m(7);
t119 = mrSges(6,3) + mrSges(7,2);
t44 = -t112 * t72 + t70 * t75;
t45 = t112 * t70 + t72 * t75;
t29 = t44 * mrSges(5,1) + mrSges(5,2) * t45;
t47 = -t72 * mrSges(4,1) + mrSges(4,2) * t70;
t118 = -t29 - t47;
t117 = -m(7) * pkin(5) - mrSges(7,1);
t17 = -t112 * t38 + t39 * t75;
t16 = t17 ^ 2;
t99 = pkin(8) + qJ(3);
t48 = t99 * t72;
t90 = t99 * t70;
t30 = t112 * t90 + t75 * t48;
t116 = t30 ^ 2;
t67 = t72 ^ 2;
t115 = 0.2e1 * t30;
t114 = -m(4) - m(5);
t111 = Ifges(6,4) * t74;
t110 = Ifges(6,4) * t77;
t109 = Ifges(7,5) * t74;
t108 = Ifges(7,5) * t77;
t107 = Ifges(6,6) * t44;
t106 = Ifges(7,6) * t44;
t104 = t30 * t17;
t103 = t45 * t74;
t102 = t45 * t77;
t24 = -mrSges(6,2) * t44 - mrSges(6,3) * t103;
t27 = -mrSges(7,2) * t103 + mrSges(7,3) * t44;
t98 = t24 + t27;
t25 = mrSges(6,1) * t44 - mrSges(6,3) * t102;
t26 = -mrSges(7,1) * t44 + mrSges(7,2) * t102;
t97 = t25 - t26;
t58 = -pkin(3) * t72 - pkin(2);
t23 = pkin(4) * t44 - pkin(9) * t45 + t58;
t32 = t112 * t48 - t75 * t90;
t4 = t23 * t74 + t32 * t77;
t50 = -t77 * mrSges(6,1) + t74 * mrSges(6,2);
t96 = t50 - mrSges(5,1);
t95 = t123 * pkin(9) ^ 2;
t94 = t70 ^ 2 + t67;
t89 = t122 * pkin(9);
t87 = Ifges(7,6) * t103 + t124 * t44 + (Ifges(7,4) + Ifges(6,5)) * t102;
t85 = t74 * mrSges(6,1) + t77 * mrSges(6,2);
t49 = -t77 * mrSges(7,1) - t74 * mrSges(7,3);
t84 = t74 * mrSges(7,1) - t77 * mrSges(7,3);
t83 = pkin(5) * t77 + qJ(6) * t74;
t82 = pkin(5) * t74 - qJ(6) * t77;
t3 = t23 * t77 - t32 * t74;
t81 = -t38 * t70 + t39 * t72;
t66 = t71 ^ 2;
t62 = Ifges(7,4) * t74;
t61 = Ifges(6,5) * t74;
t60 = Ifges(6,6) * t77;
t56 = t66 * t78 ^ 2;
t54 = Ifges(6,1) * t74 + t110;
t53 = Ifges(7,1) * t74 - t108;
t52 = Ifges(6,2) * t77 + t111;
t51 = -Ifges(7,3) * t77 + t109;
t46 = -pkin(4) - t83;
t22 = t85 * t45;
t21 = t84 * t45;
t15 = Ifges(6,5) * t44 + (Ifges(6,1) * t77 - t111) * t45;
t14 = Ifges(7,4) * t44 + (Ifges(7,1) * t77 + t109) * t45;
t13 = t107 + (-Ifges(6,2) * t74 + t110) * t45;
t12 = t106 + (Ifges(7,3) * t74 + t108) * t45;
t5 = t45 * t82 + t30;
t2 = -pkin(5) * t44 - t3;
t1 = qJ(6) * t44 + t4;
t6 = [m(2) + m(5) * (t19 ^ 2 + t16 + t56) + m(4) * (t38 ^ 2 + t39 ^ 2 + t56) + m(3) * (t66 * t76 ^ 2 + t73 ^ 2 + t56) + t120 * (t11 ^ 2 + t9 ^ 2 + t16); -t19 * t44 * mrSges(5,3) - t97 * t9 + t98 * t11 + t81 * mrSges(4,3) + (t45 * mrSges(5,3) + t21 + t22) * t17 + (-t76 * mrSges(3,2) + (mrSges(3,1) + t118) * t78) * t71 + m(6) * (t11 * t4 - t3 * t9 + t104) + m(7) * (t1 * t11 + t17 * t5 + t2 * t9) + m(5) * (-t100 * t58 + t19 * t32 + t104) + m(4) * (pkin(2) * t100 + qJ(3) * t81); Ifges(4,2) * t67 - 0.2e1 * pkin(2) * t47 + 0.2e1 * t1 * t27 + 0.2e1 * t2 * t26 + 0.2e1 * t5 * t21 + t22 * t115 + 0.2e1 * t4 * t24 + 0.2e1 * t3 * t25 + 0.2e1 * t58 * t29 + Ifges(3,3) + (Ifges(4,1) * t70 + 0.2e1 * Ifges(4,4) * t72) * t70 + 0.2e1 * t94 * qJ(3) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t32 + Ifges(5,2) * t44 + t87) * t44 + m(4) * (qJ(3) ^ 2 * t94 + pkin(2) ^ 2) + m(5) * (t32 ^ 2 + t58 ^ 2 + t116) + m(6) * (t3 ^ 2 + t4 ^ 2 + t116) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (mrSges(5,3) * t115 + Ifges(5,1) * t45 - 0.2e1 * Ifges(5,4) * t44 + (t14 + t15) * t77 + (t12 - t13 - t107) * t74) * t45; t114 * t100 + t120 * (t11 * t74 - t77 * t9); -m(4) * pkin(2) + t97 * t77 + t98 * t74 + m(7) * (t1 * t74 - t2 * t77) + m(6) * (t3 * t77 + t4 * t74) + m(5) * t58 - t118; (m(6) / 0.2e1 + m(7) / 0.2e1) * t121 - t114; -t19 * mrSges(5,2) + (t49 + t96) * t17 + m(6) * (-pkin(4) * t17 + t89) + m(7) * (t17 * t46 + t89) + t119 * t122; -t32 * mrSges(5,2) - pkin(4) * t22 + t46 * t21 + t5 * t49 + (t62 / 0.2e1 + t61 / 0.2e1 + t60 / 0.2e1 - Ifges(5,6)) * t44 + t96 * t30 + (-t106 / 0.2e1 - t12 / 0.2e1 + t13 / 0.2e1 + t4 * mrSges(6,3) + t1 * mrSges(7,2) + t98 * pkin(9)) * t77 + (t14 / 0.2e1 + t15 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3) - t97 * pkin(9)) * t74 + m(7) * (t46 * t5 + (t1 * t77 + t2 * t74) * pkin(9)) + m(6) * (-pkin(4) * t30 + (-t3 * t74 + t4 * t77) * pkin(9)) + (Ifges(5,5) + (t53 / 0.2e1 + t54 / 0.2e1) * t77 + (t51 / 0.2e1 - t52 / 0.2e1) * t74) * t45; 0; -0.2e1 * pkin(4) * t50 + 0.2e1 * t46 * t49 + Ifges(5,3) + (-t51 + t52) * t77 + (t54 + t53) * t74 + m(7) * (t46 ^ 2 + t95) + m(6) * (pkin(4) ^ 2 + t95) + t119 * pkin(9) * t121; (-mrSges(6,1) + t117) * t9 + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t11; -Ifges(6,6) * t103 - pkin(5) * t26 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t27 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t87; m(7) * t83 - t49 - t50; -Ifges(7,6) * t77 + t60 + t61 + t62 - t82 * mrSges(7,2) + (-m(7) * t82 - t84 - t85) * pkin(9); 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t124; m(7) * t9; m(7) * t2 + t26; -m(7) * t77; (m(7) * pkin(9) + mrSges(7,2)) * t74; t117; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
