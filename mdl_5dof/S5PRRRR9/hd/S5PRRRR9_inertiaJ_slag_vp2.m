% Calculate joint inertia matrix for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:05
% EndTime: 2019-12-05 17:19:07
% DurationCPUTime: 0.60s
% Computational Cost: add. (661->196), mult. (1496->294), div. (0->0), fcn. (1469->10), ass. (0->87)
t99 = 2 * pkin(7);
t65 = cos(pkin(5));
t68 = sin(qJ(3));
t72 = cos(qJ(3));
t64 = sin(pkin(5));
t69 = sin(qJ(2));
t91 = t64 * t69;
t31 = -t65 * t72 + t68 * t91;
t30 = t31 ^ 2;
t71 = cos(qJ(4));
t98 = t71 / 0.2e1;
t97 = -pkin(9) - pkin(8);
t96 = pkin(7) * t72;
t67 = sin(qJ(4));
t95 = Ifges(5,4) * t67;
t94 = Ifges(5,4) * t71;
t93 = t31 * t68;
t33 = t65 * t68 + t72 * t91;
t92 = t33 * t72;
t73 = cos(qJ(2));
t90 = t64 * t73;
t89 = t67 * t68;
t88 = t68 * t71;
t87 = -Ifges(6,3) - Ifges(5,3);
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t40 = t66 * t71 + t70 * t67;
t28 = t40 * t68;
t39 = -t66 * t67 + t70 * t71;
t29 = t39 * t68;
t86 = -Ifges(6,5) * t29 + Ifges(6,6) * t28;
t45 = -t71 * mrSges(5,1) + t67 * mrSges(5,2);
t85 = t45 - mrSges(4,1);
t44 = -t72 * pkin(3) - t68 * pkin(8) - pkin(2);
t23 = t67 * t44 + t71 * t96;
t84 = t67 ^ 2 + t71 ^ 2;
t14 = -t33 * t67 - t71 * t90;
t15 = t33 * t71 - t67 * t90;
t5 = t70 * t14 - t66 * t15;
t6 = t66 * t14 + t70 * t15;
t83 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t82 = t67 * mrSges(5,1) + t71 * mrSges(5,2);
t81 = -t14 * t67 + t15 * t71;
t38 = t71 * t44;
t22 = -t67 * t96 + t38;
t80 = -t22 * t67 + t23 * t71;
t10 = -pkin(9) * t88 + t38 + (-pkin(7) * t67 - pkin(4)) * t72;
t16 = -pkin(9) * t89 + t23;
t2 = t70 * t10 - t66 * t16;
t3 = t66 * t10 + t70 * t16;
t79 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t86;
t49 = t97 * t67;
t50 = t97 * t71;
t18 = t70 * t49 + t66 * t50;
t19 = t66 * t49 - t70 * t50;
t35 = Ifges(6,6) * t39;
t36 = Ifges(6,5) * t40;
t78 = t18 * mrSges(6,1) - t19 * mrSges(6,2) + t35 + t36;
t77 = (t70 * mrSges(6,1) - t66 * mrSges(6,2)) * pkin(4);
t75 = pkin(7) ^ 2;
t63 = t72 ^ 2;
t61 = t68 ^ 2;
t59 = t64 ^ 2;
t58 = t61 * t75;
t57 = Ifges(5,5) * t67;
t56 = Ifges(5,6) * t71;
t55 = -t71 * pkin(4) - pkin(3);
t53 = t59 * t73 ^ 2;
t51 = Ifges(5,5) * t88;
t48 = Ifges(5,1) * t67 + t94;
t47 = Ifges(5,2) * t71 + t95;
t46 = -t72 * mrSges(4,1) + t68 * mrSges(4,2);
t43 = (pkin(4) * t67 + pkin(7)) * t68;
t42 = -t72 * mrSges(5,1) - mrSges(5,3) * t88;
t41 = t72 * mrSges(5,2) - mrSges(5,3) * t89;
t34 = t82 * t68;
t27 = -Ifges(5,5) * t72 + (Ifges(5,1) * t71 - t95) * t68;
t26 = -Ifges(5,6) * t72 + (-Ifges(5,2) * t67 + t94) * t68;
t21 = -t72 * mrSges(6,1) - t29 * mrSges(6,3);
t20 = t72 * mrSges(6,2) - t28 * mrSges(6,3);
t13 = Ifges(6,1) * t40 + Ifges(6,4) * t39;
t12 = Ifges(6,4) * t40 + Ifges(6,2) * t39;
t11 = -t39 * mrSges(6,1) + t40 * mrSges(6,2);
t9 = t28 * mrSges(6,1) + t29 * mrSges(6,2);
t8 = Ifges(6,1) * t29 - Ifges(6,4) * t28 - Ifges(6,5) * t72;
t7 = Ifges(6,4) * t29 - Ifges(6,2) * t28 - Ifges(6,6) * t72;
t1 = [m(2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t30) + m(5) * (t14 ^ 2 + t15 ^ 2 + t30) + m(4) * (t33 ^ 2 + t30 + t53) + m(3) * (t59 * t69 ^ 2 + t65 ^ 2 + t53); mrSges(4,3) * t92 + t14 * t42 + t15 * t41 + t6 * t20 + t5 * t21 + (-t69 * mrSges(3,2) + (mrSges(3,1) - t46) * t73) * t64 + (t68 * mrSges(4,3) + t34 + t9) * t31 + m(6) * (t2 * t5 + t3 * t6 + t43 * t31) + m(5) * (pkin(7) * t93 + t22 * t14 + t23 * t15) + m(4) * (pkin(2) * t90 + (t92 + t93) * pkin(7)); -0.2e1 * pkin(2) * t46 + 0.2e1 * t2 * t21 + 0.2e1 * t3 * t20 + 0.2e1 * t22 * t42 + 0.2e1 * t23 * t41 - t28 * t7 + t29 * t8 + 0.2e1 * t43 * t9 + Ifges(3,3) + (t61 + t63) * mrSges(4,3) * t99 + m(6) * (t2 ^ 2 + t3 ^ 2 + t43 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2 + t58) + m(4) * (pkin(2) ^ 2 + t63 * t75 + t58) + (Ifges(4,1) * t68 - t67 * t26 + t71 * t27 + t34 * t99) * t68 + (-t51 + (Ifges(4,2) - t87) * t72 + (Ifges(5,6) * t67 + (2 * Ifges(4,4))) * t68 + t86) * t72; -t33 * mrSges(4,2) + (t6 * t39 - t5 * t40) * mrSges(6,3) + t81 * mrSges(5,3) + (t11 + t85) * t31 + m(6) * (t18 * t5 + t19 * t6 + t55 * t31) + m(5) * (-pkin(3) * t31 + t81 * pkin(8)); t26 * t98 + t67 * t27 / 0.2e1 + t55 * t9 + t39 * t7 / 0.2e1 + t40 * t8 / 0.2e1 + t43 * t11 - t28 * t12 / 0.2e1 + t29 * t13 / 0.2e1 - pkin(3) * t34 + t19 * t20 + t18 * t21 + m(6) * (t18 * t2 + t19 * t3 + t55 * t43) + (-t57 / 0.2e1 - t56 / 0.2e1 - t36 / 0.2e1 - t35 / 0.2e1 + Ifges(4,6) - pkin(7) * mrSges(4,2)) * t72 + (-t2 * t40 + t3 * t39) * mrSges(6,3) + t80 * mrSges(5,3) + (m(5) * t80 + t71 * t41 - t67 * t42) * pkin(8) + (Ifges(4,5) - t67 * t47 / 0.2e1 + t48 * t98 + (-m(5) * pkin(3) + t85) * pkin(7)) * t68; -0.2e1 * pkin(3) * t45 + 0.2e1 * t55 * t11 + t39 * t12 + t40 * t13 + t71 * t47 + t67 * t48 + Ifges(4,3) + m(6) * (t18 ^ 2 + t19 ^ 2 + t55 ^ 2) + m(5) * (t84 * pkin(8) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t18 * t40 + t19 * t39) * mrSges(6,3) + 0.2e1 * t84 * pkin(8) * mrSges(5,3); t14 * mrSges(5,1) - t15 * mrSges(5,2) + m(6) * (t5 * t70 + t6 * t66) * pkin(4) + t83; -Ifges(5,6) * t89 + t22 * mrSges(5,1) - t23 * mrSges(5,2) + t51 + t87 * t72 + (m(6) * (t2 * t70 + t3 * t66) + t66 * t20 + t70 * t21) * pkin(4) + t79; t56 + t57 - t82 * pkin(8) + (m(6) * (t18 * t70 + t19 * t66) + (t66 * t39 - t70 * t40) * mrSges(6,3)) * pkin(4) + t78; m(6) * (t66 ^ 2 + t70 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t77 - t87; t83; -Ifges(6,3) * t72 + t79; t78; Ifges(6,3) + t77; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
