% Calculate joint inertia matrix for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:04
% EndTime: 2019-03-09 03:18:05
% DurationCPUTime: 0.66s
% Computational Cost: add. (957->207), mult. (1750->263), div. (0->0), fcn. (1770->6), ass. (0->83)
t105 = Ifges(6,6) + Ifges(7,6);
t104 = Ifges(6,3) + Ifges(7,3);
t66 = cos(pkin(9));
t61 = t66 ^ 2;
t103 = -2 * mrSges(7,3);
t65 = sin(pkin(9));
t68 = sin(qJ(3));
t98 = cos(qJ(3));
t37 = t68 * t65 - t98 * t66;
t38 = t98 * t65 + t68 * t66;
t54 = -t66 * pkin(2) - pkin(1);
t73 = -t38 * qJ(4) + t54;
t16 = t37 * pkin(3) + t73;
t102 = -0.2e1 * t16;
t101 = 0.2e1 * t54;
t100 = m(7) * pkin(5);
t99 = pkin(3) + pkin(8);
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t51 = t67 ^ 2 + t69 ^ 2;
t40 = m(6) * t51;
t86 = pkin(7) + qJ(2);
t43 = t86 * t65;
t44 = t86 * t66;
t23 = t98 * t43 + t68 * t44;
t14 = t38 * pkin(4) + t23;
t7 = t99 * t37 + t73;
t4 = t67 * t14 + t69 * t7;
t97 = Ifges(6,4) * t67;
t96 = Ifges(6,4) * t69;
t95 = Ifges(7,4) * t67;
t94 = Ifges(7,4) * t69;
t93 = t37 * t67;
t92 = t37 * t69;
t91 = t69 * mrSges(6,1);
t90 = t69 * mrSges(6,2);
t89 = t69 * mrSges(7,3);
t88 = mrSges(5,1) + mrSges(4,3);
t87 = mrSges(5,2) - mrSges(4,1);
t19 = t38 * mrSges(7,1) - mrSges(7,3) * t93;
t20 = t38 * mrSges(6,1) - mrSges(6,3) * t93;
t85 = t19 + t20;
t21 = -t38 * mrSges(7,2) + t37 * t89;
t22 = -t38 * mrSges(6,2) + mrSges(6,3) * t92;
t84 = t21 + t22;
t45 = t67 * mrSges(7,1) + t69 * mrSges(7,2);
t83 = t65 ^ 2 + t61;
t82 = qJ(6) * t37;
t81 = -qJ(6) - t99;
t80 = m(7) * t51 + m(5) + t40;
t25 = -t68 * t43 + t98 * t44;
t79 = t23 ^ 2 + t25 ^ 2;
t78 = -mrSges(6,1) - t100;
t77 = t51 * mrSges(6,3);
t76 = -t66 * mrSges(3,1) + t65 * mrSges(3,2);
t17 = -mrSges(7,1) * t92 + mrSges(7,2) * t93;
t13 = t69 * t14;
t3 = -t67 * t7 + t13;
t75 = t69 * t3 + t67 * t4;
t74 = (Ifges(6,5) + Ifges(7,5)) * t93 + t105 * t92 + t104 * t38;
t71 = qJ(4) ^ 2;
t59 = Ifges(6,5) * t69;
t58 = Ifges(7,5) * t69;
t53 = t67 * pkin(5) + qJ(4);
t50 = Ifges(6,1) * t69 - t97;
t49 = Ifges(7,1) * t69 - t95;
t48 = -Ifges(6,2) * t67 + t96;
t47 = -Ifges(7,2) * t67 + t94;
t46 = t67 * mrSges(6,1) + t90;
t42 = t81 * t69;
t41 = t81 * t67;
t33 = t38 * mrSges(5,3);
t32 = t38 * mrSges(4,2);
t18 = (t67 * mrSges(6,2) - t91) * t37;
t15 = -t37 * pkin(4) + t25;
t11 = Ifges(6,5) * t38 + (Ifges(6,1) * t67 + t96) * t37;
t10 = Ifges(7,5) * t38 + (Ifges(7,1) * t67 + t94) * t37;
t9 = Ifges(6,6) * t38 + (Ifges(6,2) * t69 + t97) * t37;
t8 = Ifges(7,6) * t38 + (Ifges(7,2) * t69 + t95) * t37;
t5 = (-pkin(5) * t69 - pkin(4)) * t37 + t25;
t2 = t69 * t82 + t4;
t1 = t38 * pkin(5) + t13 + (-t7 - t82) * t67;
t6 = [t33 * t102 + t32 * t101 + Ifges(3,2) * t61 - 0.2e1 * pkin(1) * t76 + 0.2e1 * t5 * t17 + 0.2e1 * t15 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + 0.2e1 * t4 * t22 + Ifges(2,3) + (Ifges(3,1) * t65 + 0.2e1 * Ifges(3,4) * t66) * t65 + 0.2e1 * t83 * qJ(2) * mrSges(3,3) + m(3) * (t83 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t54 ^ 2 + t79) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t15 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t16 ^ 2 + t79) + ((Ifges(4,1) + Ifges(5,2)) * t38 + 0.2e1 * t88 * t23 + t74) * t38 + (mrSges(4,1) * t101 + mrSges(5,2) * t102 + (Ifges(4,2) + Ifges(5,3)) * t37 + (t8 + t9) * t69 + (t10 + t11) * t67 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t38 - 0.2e1 * t88 * t25) * t37; -m(3) * pkin(1) + t32 - t33 + t84 * t69 - t85 * t67 - t87 * t37 + m(7) * (-t67 * t1 + t69 * t2) + m(6) * (-t67 * t3 + t69 * t4) + m(5) * t16 + m(4) * t54 + t76; m(3) + m(4) + t80; qJ(4) * t18 + t15 * t46 + t53 * t17 + t42 * t19 + t41 * t21 + t5 * t45 + (mrSges(5,3) - mrSges(4,2)) * t25 + t87 * t23 + (-t99 * t22 - t2 * mrSges(7,3) - t4 * mrSges(6,3) - t8 / 0.2e1 - t9 / 0.2e1) * t67 + (-pkin(3) * mrSges(5,1) + t58 / 0.2e1 + t59 / 0.2e1 + Ifges(4,5) - Ifges(5,4) + (-Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t67) * t38 + (-t99 * t20 - t1 * mrSges(7,3) - t3 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1) * t69 + m(6) * (qJ(4) * t15 - t75 * t99) + m(7) * (t42 * t1 + t41 * t2 + t53 * t5) + m(5) * (-pkin(3) * t23 + qJ(4) * t25) + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + (t47 / 0.2e1 + t48 / 0.2e1) * t69 + (t49 / 0.2e1 + t50 / 0.2e1) * t67) * t37; m(7) * (t41 * t69 - t42 * t67); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t53 * t45 + Ifges(5,1) + Ifges(4,3) + (t42 * t103 + t49 + t50) * t69 + (t41 * t103 - t47 - t48) * t67 + m(7) * (t41 ^ 2 + t42 ^ 2 + t53 ^ 2) + m(6) * (t51 * t99 ^ 2 + t71) + m(5) * (pkin(3) ^ 2 + t71) + 0.2e1 * (t46 + mrSges(5,3)) * qJ(4) + 0.2e1 * t99 * t77; t38 * mrSges(5,1) + t85 * t69 + t84 * t67 + m(7) * (t69 * t1 + t67 * t2) + m(6) * t75 + m(5) * t23; 0; -m(5) * pkin(3) + mrSges(5,2) - t51 * mrSges(7,3) - t77 + m(7) * (t67 * t41 + t69 * t42) - t99 * t40; t80; t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (m(7) * t1 + t19) * pkin(5) + t74; t78 * t67 - t45 - t90; -t99 * t91 + t42 * mrSges(7,1) - t41 * mrSges(7,2) + t58 + t59 + (m(7) * t42 - t89) * pkin(5) + (mrSges(6,2) * t99 - t105) * t67; (-mrSges(6,2) - mrSges(7,2)) * t67 + (mrSges(7,1) - t78) * t69; (0.2e1 * mrSges(7,1) + t100) * pkin(5) + t104; m(7) * t5 + t17; 0; m(7) * t53 + t45; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
