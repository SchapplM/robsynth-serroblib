% Calculate joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:23
% EndTime: 2019-12-31 17:29:25
% DurationCPUTime: 0.58s
% Computational Cost: add. (662->184), mult. (1564->281), div. (0->0), fcn. (1523->8), ass. (0->80)
t100 = 2 * pkin(7);
t69 = cos(pkin(4));
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t68 = sin(pkin(4));
t72 = sin(qJ(2));
t88 = t68 * t72;
t37 = t69 * t71 + t74 * t88;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t75 = cos(qJ(2));
t87 = t68 * t75;
t19 = -t37 * t70 - t73 * t87;
t99 = t19 / 0.2e1;
t20 = t37 * t73 - t70 * t87;
t98 = t20 / 0.2e1;
t97 = -t70 / 0.2e1;
t96 = t70 / 0.2e1;
t95 = t73 / 0.2e1;
t94 = pkin(1) * t75;
t93 = pkin(7) * t74;
t92 = Ifges(5,4) * t70;
t91 = Ifges(5,4) * t73;
t53 = pkin(6) * t88;
t38 = t69 * t94 - t53;
t90 = t38 * mrSges(3,1);
t39 = t69 * t72 * pkin(1) + pkin(6) * t87;
t89 = t39 * mrSges(3,2);
t86 = t70 * t71;
t85 = t71 * t73;
t28 = t69 * pkin(7) + t39;
t29 = (-pkin(2) * t75 - pkin(7) * t72 - pkin(1)) * t68;
t13 = t74 * t28 + t71 * t29;
t36 = -t69 * t74 + t71 * t88;
t84 = -Ifges(4,5) * t37 + Ifges(4,6) * t36;
t46 = Ifges(5,5) * t70 + Ifges(5,6) * t73;
t83 = Ifges(4,5) * t71 + Ifges(4,6) * t74;
t82 = t70 ^ 2 + t73 ^ 2;
t3 = Ifges(5,5) * t20 + Ifges(5,6) * t19 + Ifges(5,3) * t36;
t81 = Ifges(3,5) * t88 + Ifges(3,6) * t87 + Ifges(3,3) * t69;
t27 = t53 + (-pkin(2) - t94) * t69;
t7 = t36 * pkin(3) - t37 * pkin(8) + t27;
t9 = -pkin(8) * t87 + t13;
t1 = t73 * t7 - t70 * t9;
t2 = t70 * t7 + t73 * t9;
t80 = -t1 * t70 + t2 * t73;
t79 = mrSges(5,1) * t70 + mrSges(5,2) * t73;
t43 = -t74 * pkin(3) - t71 * pkin(8) - pkin(2);
t25 = t73 * t43 - t70 * t93;
t26 = t70 * t43 + t73 * t93;
t78 = -t25 * t70 + t26 * t73;
t12 = -t71 * t28 + t74 * t29;
t30 = Ifges(5,5) * t85 - Ifges(5,6) * t86 - Ifges(5,3) * t74;
t77 = pkin(7) ^ 2;
t67 = t74 ^ 2;
t65 = t71 ^ 2;
t63 = t65 * t77;
t50 = Ifges(4,1) * t71 + Ifges(4,4) * t74;
t49 = Ifges(5,1) * t70 + t91;
t48 = Ifges(4,4) * t71 + Ifges(4,2) * t74;
t47 = Ifges(5,2) * t73 + t92;
t45 = -t74 * mrSges(4,1) + t71 * mrSges(4,2);
t44 = -t73 * mrSges(5,1) + t70 * mrSges(5,2);
t42 = -t74 * mrSges(5,1) - mrSges(5,3) * t85;
t41 = t74 * mrSges(5,2) - mrSges(5,3) * t86;
t40 = t79 * t71;
t32 = -Ifges(5,5) * t74 + (Ifges(5,1) * t73 - t92) * t71;
t31 = -Ifges(5,6) * t74 + (-Ifges(5,2) * t70 + t91) * t71;
t22 = -mrSges(4,1) * t87 - t37 * mrSges(4,3);
t21 = mrSges(4,2) * t87 - t36 * mrSges(4,3);
t16 = t36 * mrSges(4,1) + t37 * mrSges(4,2);
t15 = Ifges(4,1) * t37 - Ifges(4,4) * t36 - Ifges(4,5) * t87;
t14 = Ifges(4,4) * t37 - Ifges(4,2) * t36 - Ifges(4,6) * t87;
t11 = t36 * mrSges(5,1) - t20 * mrSges(5,3);
t10 = -t36 * mrSges(5,2) + t19 * mrSges(5,3);
t8 = pkin(3) * t87 - t12;
t6 = -t19 * mrSges(5,1) + t20 * mrSges(5,2);
t5 = Ifges(5,1) * t20 + Ifges(5,4) * t19 + Ifges(5,5) * t36;
t4 = Ifges(5,4) * t20 + Ifges(5,2) * t19 + Ifges(5,6) * t36;
t17 = [0.2e1 * t1 * t11 + 0.2e1 * t2 * t10 + 0.2e1 * t12 * t22 + 0.2e1 * t13 * t21 + t37 * t15 + 0.2e1 * t27 * t16 + t19 * t4 + t20 * t5 + 0.2e1 * t8 * t6 + Ifges(2,3) + (t3 - t14) * t36 + (t81 - 0.2e1 * t89 + 0.2e1 * t90) * t69 + ((-0.2e1 * t38 * mrSges(3,3) + Ifges(3,5) * t69 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t72) * t68) * t72 + (0.2e1 * t39 * mrSges(3,3) + Ifges(3,6) * t69 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t72 + (Ifges(4,3) + Ifges(3,2)) * t75) * t68 + t84) * t75) * t68 + m(5) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) + m(3) * (t68 ^ 2 * pkin(1) ^ 2 + t38 ^ 2 + t39 ^ 2); t81 + m(4) * (-pkin(2) * t27 + (-t12 * t71 + t13 * t74) * pkin(7)) + (t15 / 0.2e1 - t12 * mrSges(4,3) + t4 * t97 + t5 * t95 + (t6 - t22) * pkin(7)) * t71 + (t14 / 0.2e1 - t3 / 0.2e1 + t13 * mrSges(4,3) + pkin(7) * t21) * t74 + m(5) * (t71 * pkin(7) * t8 + t25 * t1 + t26 * t2) + (-t48 / 0.2e1 + t30 / 0.2e1) * t36 + t1 * t42 + t27 * t45 + t37 * t50 / 0.2e1 + t90 - t89 + t8 * t40 + t2 * t41 + t25 * t11 + t26 * t10 + t31 * t99 + t32 * t98 - pkin(2) * t16 - t83 * t87 / 0.2e1; -0.2e1 * pkin(2) * t45 + 0.2e1 * t25 * t42 + 0.2e1 * t26 * t41 + Ifges(3,3) + (-t30 + t48) * t74 + (t65 + t67) * mrSges(4,3) * t100 + m(5) * (t25 ^ 2 + t26 ^ 2 + t63) + m(4) * (pkin(2) ^ 2 + t67 * t77 + t63) + (t40 * t100 - t70 * t31 + t73 * t32 + t50) * t71; -Ifges(4,3) * t87 + t8 * t44 + t5 * t96 + t4 * t95 + t49 * t98 + t47 * t99 + t36 * t46 / 0.2e1 - t13 * mrSges(4,2) + t12 * mrSges(4,1) + t80 * mrSges(5,3) - t84 + (-m(5) * t8 - t6) * pkin(3) + (m(5) * t80 + t73 * t10 - t70 * t11) * pkin(8); -pkin(3) * t40 + t32 * t96 + t31 * t95 + (m(5) * t78 + t73 * t41 - t70 * t42) * pkin(8) + (t49 * t95 + t47 * t97 + (-m(5) * pkin(3) - mrSges(4,1) + t44) * pkin(7)) * t71 + (-pkin(7) * mrSges(4,2) - t46 / 0.2e1) * t74 + t78 * mrSges(5,3) + t83; Ifges(4,3) + m(5) * (t82 * pkin(8) ^ 2 + pkin(3) ^ 2) - 0.2e1 * pkin(3) * t44 + t70 * t49 + t73 * t47 + 0.2e1 * t82 * pkin(8) * mrSges(5,3); t1 * mrSges(5,1) - t2 * mrSges(5,2) + t3; t25 * mrSges(5,1) - t26 * mrSges(5,2) + t30; -pkin(8) * t79 + t46; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1), t17(2), t17(4), t17(7); t17(2), t17(3), t17(5), t17(8); t17(4), t17(5), t17(6), t17(9); t17(7), t17(8), t17(9), t17(10);];
Mq = res;
