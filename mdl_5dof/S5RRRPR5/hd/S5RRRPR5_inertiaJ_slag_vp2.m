% Calculate joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:02
% EndTime: 2019-12-31 21:13:04
% DurationCPUTime: 0.59s
% Computational Cost: add. (1165->168), mult. (2210->251), div. (0->0), fcn. (2381->8), ass. (0->69)
t66 = sin(qJ(3));
t67 = sin(qJ(2));
t69 = cos(qJ(3));
t70 = cos(qJ(2));
t43 = -t66 * t67 + t69 * t70;
t44 = t66 * t70 + t69 * t67;
t63 = sin(pkin(9));
t64 = cos(pkin(9));
t27 = -t64 * t43 + t63 * t44;
t28 = t63 * t43 + t64 * t44;
t65 = sin(qJ(5));
t91 = t65 * mrSges(6,3);
t15 = -t27 * mrSges(6,2) - t28 * t91;
t68 = cos(qJ(5));
t94 = t28 * t68;
t16 = t27 * mrSges(6,1) - mrSges(6,3) * t94;
t105 = t68 * t15 - t65 * t16;
t100 = -pkin(7) - pkin(6);
t83 = t100 * t67;
t84 = t100 * t70;
t31 = t66 * t83 - t69 * t84;
t21 = t43 * qJ(4) + t31;
t30 = t66 * t84 + t69 * t83;
t75 = -t44 * qJ(4) + t30;
t11 = t63 * t21 - t64 * t75;
t104 = t11 ^ 2;
t103 = 0.2e1 * t11;
t56 = -t70 * pkin(2) - pkin(1);
t32 = -t43 * pkin(3) + t56;
t102 = 0.2e1 * t32;
t47 = -t68 * mrSges(6,1) + t65 * mrSges(6,2);
t101 = 0.2e1 * t47;
t99 = pkin(2) * t66;
t10 = t27 * pkin(4) - t28 * pkin(8) + t32;
t13 = t64 * t21 + t63 * t75;
t3 = t65 * t10 + t68 * t13;
t98 = t3 * t68;
t97 = Ifges(6,4) * t65;
t96 = Ifges(6,4) * t68;
t95 = t28 * t65;
t55 = t69 * pkin(2) + pkin(3);
t36 = t64 * t55 - t63 * t99;
t93 = t36 * mrSges(5,1);
t37 = t63 * t55 + t64 * t99;
t92 = t37 * mrSges(5,2);
t88 = Ifges(6,5) * t94 + Ifges(6,3) * t27;
t87 = Ifges(6,5) * t65 + Ifges(6,6) * t68;
t86 = t65 ^ 2 + t68 ^ 2;
t85 = t67 ^ 2 + t70 ^ 2;
t53 = t63 * pkin(3) + pkin(8);
t82 = t86 * t53;
t48 = Ifges(6,2) * t68 + t97;
t49 = Ifges(6,1) * t65 + t96;
t81 = t68 * t48 + t65 * t49 + Ifges(4,3) + Ifges(5,3);
t2 = t68 * t10 - t65 * t13;
t80 = -t2 * t65 + t98;
t79 = t64 * mrSges(5,1) - t63 * mrSges(5,2);
t78 = mrSges(6,1) * t65 + mrSges(6,2) * t68;
t77 = 0.2e1 * t86 * mrSges(6,3);
t76 = (t69 * mrSges(4,1) - t66 * mrSges(4,2)) * pkin(2);
t6 = Ifges(6,6) * t27 + (-Ifges(6,2) * t65 + t96) * t28;
t7 = Ifges(6,5) * t27 + (Ifges(6,1) * t68 - t97) * t28;
t74 = -t31 * mrSges(4,2) - t13 * mrSges(5,2) + mrSges(6,3) * t98 - t2 * t91 + t65 * t7 / 0.2e1 - t48 * t95 / 0.2e1 + t49 * t94 / 0.2e1 + Ifges(5,5) * t28 + t30 * mrSges(4,1) + Ifges(4,6) * t43 + Ifges(4,5) * t44 + t68 * t6 / 0.2e1 + (t87 / 0.2e1 - Ifges(5,6)) * t27 + (t47 - mrSges(5,1)) * t11;
t54 = -t64 * pkin(3) - pkin(4);
t35 = pkin(8) + t37;
t34 = -pkin(4) - t36;
t23 = t28 * mrSges(5,2);
t14 = t78 * t28;
t1 = [t67 * (Ifges(3,1) * t67 + Ifges(3,4) * t70) + t70 * (Ifges(3,4) * t67 + Ifges(3,2) * t70) - 0.2e1 * pkin(1) * (-t70 * mrSges(3,1) + t67 * mrSges(3,2)) + t43 * (Ifges(4,4) * t44 + Ifges(4,2) * t43) + 0.2e1 * t56 * (-t43 * mrSges(4,1) + t44 * mrSges(4,2)) + t23 * t102 + t44 * (Ifges(4,1) * t44 + Ifges(4,4) * t43) + t14 * t103 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (mrSges(5,1) * t102 - 0.2e1 * t13 * mrSges(5,3) + Ifges(5,2) * t27 + t88) * t27 + (mrSges(5,3) * t103 + Ifges(5,1) * t28 - t65 * t6 + t68 * t7 + (-Ifges(6,6) * t65 - (2 * Ifges(5,4))) * t27) * t28 + m(3) * (t85 * pkin(6) ^ 2 + pkin(1) ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2 + t56 ^ 2) + m(5) * (t13 ^ 2 + t32 ^ 2 + t104) + m(6) * (t2 ^ 2 + t3 ^ 2 + t104) + 0.2e1 * (-t30 * t44 + t31 * t43) * mrSges(4,3) + 0.2e1 * t85 * pkin(6) * mrSges(3,3); (m(4) * (t30 * t69 + t31 * t66) + (t66 * t43 - t69 * t44) * mrSges(4,3)) * pkin(2) + m(6) * (t34 * t11 + t80 * t35) + t105 * t35 + m(5) * (-t36 * t11 + t37 * t13) + (-t67 * mrSges(3,1) - t70 * mrSges(3,2)) * pkin(6) + (-t37 * t27 - t36 * t28) * mrSges(5,3) + t74 + Ifges(3,5) * t67 + Ifges(3,6) * t70 + t34 * t14; 0.2e1 * t93 - 0.2e1 * t92 + t34 * t101 + Ifges(3,3) + 0.2e1 * t76 + t35 * t77 + m(5) * (t36 ^ 2 + t37 ^ 2) + m(6) * (t86 * t35 ^ 2 + t34 ^ 2) + m(4) * (t66 ^ 2 + t69 ^ 2) * pkin(2) ^ 2 + t81; (m(5) * (-t11 * t64 + t13 * t63) + (-t63 * t27 - t64 * t28) * mrSges(5,3)) * pkin(3) + t74 + (m(6) * t11 + t14) * t54 + (m(6) * t80 + t105) * t53; m(6) * (t54 * t34 + t35 * t82) - t92 + t93 + (t34 + t54) * t47 + t76 + (m(5) * (t36 * t64 + t37 * t63) + t79) * pkin(3) + (t86 * t35 + t82) * mrSges(6,3) + t81; t54 * t101 + t53 * t77 + m(6) * (t86 * t53 ^ 2 + t54 ^ 2) + t81 + (0.2e1 * t79 + m(5) * (t63 ^ 2 + t64 ^ 2) * pkin(3)) * pkin(3); t27 * mrSges(5,1) + t65 * t15 + t68 * t16 + t23 + m(6) * (t68 * t2 + t65 * t3) + m(5) * t32; 0; 0; m(6) * t86 + m(5); t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t95 + t88; -t78 * t35 + t87; -t78 * t53 + t87; -t47; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
