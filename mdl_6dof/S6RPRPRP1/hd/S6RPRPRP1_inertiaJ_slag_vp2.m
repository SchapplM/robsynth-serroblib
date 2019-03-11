% Calculate joint inertia matrix for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:22
% EndTime: 2019-03-09 03:01:24
% DurationCPUTime: 0.78s
% Computational Cost: add. (994->212), mult. (1838->292), div. (0->0), fcn. (1801->8), ass. (0->89)
t112 = Ifges(6,3) + Ifges(7,3);
t71 = cos(qJ(3));
t111 = t71 ^ 2;
t110 = 2 * mrSges(7,3);
t109 = m(6) + m(7);
t65 = sin(pkin(9));
t52 = t65 * pkin(1) + pkin(7);
t83 = qJ(4) + t52;
t34 = t83 * t71;
t64 = sin(pkin(10));
t66 = cos(pkin(10));
t69 = sin(qJ(3));
t78 = t83 * t69;
t15 = t64 * t34 + t66 * t78;
t108 = t15 ^ 2;
t36 = t64 * t69 - t66 * t71;
t107 = t36 ^ 2;
t106 = 0.2e1 * t15;
t67 = cos(pkin(9));
t54 = -t67 * pkin(1) - pkin(2);
t40 = -t71 * pkin(3) + t54;
t105 = 0.2e1 * t40;
t104 = m(5) * pkin(3);
t103 = m(7) * pkin(5);
t102 = t64 * pkin(3);
t101 = t66 * pkin(3);
t70 = cos(qJ(5));
t100 = mrSges(6,2) * t70;
t68 = sin(qJ(5));
t99 = Ifges(6,4) * t68;
t98 = Ifges(6,4) * t70;
t97 = Ifges(7,4) * t68;
t96 = Ifges(7,4) * t70;
t95 = t15 * t36;
t38 = t64 * t71 + t66 * t69;
t94 = t38 * t68;
t93 = t38 * t70;
t92 = t68 * mrSges(6,2);
t91 = t68 * mrSges(7,3);
t90 = -Ifges(6,6) - Ifges(7,6);
t14 = t36 * pkin(4) - t38 * pkin(8) + t40;
t17 = t66 * t34 - t64 * t78;
t4 = t68 * t14 + t70 * t17;
t20 = -t36 * mrSges(7,2) - t38 * t91;
t21 = -t36 * mrSges(6,2) - mrSges(6,3) * t94;
t89 = t20 + t21;
t22 = t36 * mrSges(7,1) - mrSges(7,3) * t93;
t23 = t36 * mrSges(6,1) - mrSges(6,3) * t93;
t88 = t22 + t23;
t18 = mrSges(7,1) * t94 + mrSges(7,2) * t93;
t42 = -t70 * mrSges(6,1) + t92;
t87 = t42 - mrSges(5,1);
t86 = t68 ^ 2 + t70 ^ 2;
t85 = t69 ^ 2 + t111;
t84 = qJ(6) * t38;
t51 = pkin(8) + t102;
t82 = qJ(6) + t51;
t80 = -mrSges(6,1) - t103;
t53 = -pkin(4) - t101;
t79 = t86 * t51;
t55 = t68 * mrSges(7,2);
t41 = -t70 * mrSges(7,1) + t55;
t3 = t70 * t14 - t68 * t17;
t77 = (Ifges(6,5) + Ifges(7,5)) * t93 + t112 * t36;
t76 = -t3 * t68 + t4 * t70;
t75 = -t71 * mrSges(4,1) + t69 * mrSges(4,2);
t74 = mrSges(6,1) * t68 + t100;
t59 = Ifges(6,5) * t68;
t58 = Ifges(7,5) * t68;
t57 = Ifges(6,6) * t70;
t56 = Ifges(7,6) * t70;
t46 = Ifges(6,1) * t68 + t98;
t45 = Ifges(7,1) * t68 + t96;
t44 = Ifges(6,2) * t70 + t99;
t43 = Ifges(7,2) * t70 + t97;
t39 = -t70 * pkin(5) + t53;
t35 = t38 ^ 2;
t33 = t82 * t70;
t32 = t82 * t68;
t29 = t38 * mrSges(5,2);
t19 = t74 * t38;
t9 = Ifges(6,5) * t36 + (Ifges(6,1) * t70 - t99) * t38;
t8 = Ifges(7,5) * t36 + (Ifges(7,1) * t70 - t97) * t38;
t7 = Ifges(6,6) * t36 + (-Ifges(6,2) * t68 + t98) * t38;
t6 = Ifges(7,6) * t36 + (-Ifges(7,2) * t68 + t96) * t38;
t5 = pkin(5) * t94 + t15;
t2 = -t68 * t84 + t4;
t1 = t36 * pkin(5) - t70 * t84 + t3;
t10 = [0.2e1 * t54 * t75 + Ifges(4,2) * t111 + t29 * t105 + 0.2e1 * t5 * t18 + t19 * t106 + 0.2e1 * t2 * t20 + 0.2e1 * t4 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t3 * t23 + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t105 - 0.2e1 * t17 * mrSges(5,3) + Ifges(5,2) * t36 + t77) * t36 + (mrSges(5,3) * t106 + Ifges(5,1) * t38 - 0.2e1 * Ifges(5,4) * t36 + (t8 + t9) * t70 + (t90 * t36 - t6 - t7) * t68) * t38 + m(4) * (t85 * t52 ^ 2 + t54 ^ 2) + m(5) * (t17 ^ 2 + t40 ^ 2 + t108) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t108) + m(3) * (t65 ^ 2 + t67 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t69 + 0.2e1 * Ifges(4,4) * t71) * t69 + 0.2e1 * (t67 * mrSges(3,1) - t65 * mrSges(3,2)) * pkin(1) + 0.2e1 * t85 * t52 * mrSges(4,3); (t18 + t19) * t36 + (-t88 * t68 + t89 * t70) * t38 + m(7) * (t5 * t36 + (-t1 * t68 + t2 * t70) * t38) + m(5) * (t17 * t38 + t95) + m(6) * (t76 * t38 + t95); m(3) + m(5) * (t35 + t107) + m(4) * t85 + (t86 * t35 + t107) * t109; -t17 * mrSges(5,2) + Ifges(4,5) * t69 + Ifges(4,6) * t71 + t39 * t18 + t53 * t19 + t33 * t20 - t32 * t22 + t5 * t41 + (-mrSges(4,1) * t69 - mrSges(4,2) * t71) * t52 + (-mrSges(5,3) * t102 + t58 / 0.2e1 + t56 / 0.2e1 + t59 / 0.2e1 + t57 / 0.2e1 - Ifges(5,6)) * t36 + t87 * t15 + (t6 / 0.2e1 + t7 / 0.2e1 + t51 * t21 + t2 * mrSges(7,3) + t4 * mrSges(6,3)) * t70 + (t8 / 0.2e1 + t9 / 0.2e1 - t51 * t23 - t1 * mrSges(7,3) - t3 * mrSges(6,3)) * t68 + m(6) * (t53 * t15 + t76 * t51) + m(7) * (-t1 * t32 + t2 * t33 + t39 * t5) + (-t15 * t66 + t17 * t64) * t104 + (-mrSges(5,3) * t101 + Ifges(5,5) + (t45 / 0.2e1 + t46 / 0.2e1) * t70 + (-t43 / 0.2e1 - t44 / 0.2e1) * t68) * t38; -t29 + (t41 + t87) * t36 + m(7) * (t39 * t36 + (t32 * t68 + t33 * t70) * t38) + m(6) * (t53 * t36 + t38 * t79) + (-t36 * t66 + t38 * t64) * t104 - t75 + (mrSges(6,3) + mrSges(7,3)) * t38 * t86; 0.2e1 * t39 * t41 + 0.2e1 * t53 * t42 + Ifges(4,3) + Ifges(5,3) + (t33 * t110 + t43 + t44) * t70 + (t32 * t110 + t45 + t46) * t68 + m(7) * (t32 ^ 2 + t33 ^ 2 + t39 ^ 2) + m(6) * (t86 * t51 ^ 2 + t53 ^ 2) + m(5) * (t64 ^ 2 + t66 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t66 - mrSges(5,2) * t64) * pkin(3) + 0.2e1 * mrSges(6,3) * t79; t36 * mrSges(5,1) + t29 + t88 * t70 + t89 * t68 + m(7) * (t1 * t70 + t2 * t68) + m(6) * (t3 * t70 + t4 * t68) + m(5) * t40; 0; m(7) * (-t32 * t70 + t33 * t68); t86 * t109 + m(5); t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + t90 * t94 + (m(7) * t1 + t22) * pkin(5) + t77; (t80 * t68 - t100) * t38 - t18; -t32 * mrSges(7,1) - t33 * mrSges(7,2) + t56 + t57 + t58 + t59 - t74 * t51 + (-m(7) * t32 - t91) * pkin(5); -t92 - t55 + (mrSges(7,1) - t80) * t70; (0.2e1 * mrSges(7,1) + t103) * pkin(5) + t112; m(7) * t5 + t18; m(7) * t36; m(7) * t39 + t41; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
