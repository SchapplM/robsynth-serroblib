% Calculate joint inertia matrix for
% S5RRRPR3
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:19
% EndTime: 2019-12-05 18:42:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (933->153), mult. (1714->206), div. (0->0), fcn. (1730->8), ass. (0->65)
t68 = cos(qJ(3));
t96 = t68 ^ 2;
t66 = sin(qJ(2));
t54 = t66 * pkin(1) + pkin(7);
t65 = sin(qJ(3));
t39 = (-qJ(4) - t54) * t65;
t57 = t68 * qJ(4);
t40 = t68 * t54 + t57;
t62 = sin(pkin(9));
t63 = cos(pkin(9));
t17 = t63 * t39 - t62 * t40;
t42 = t62 * t68 + t63 * t65;
t89 = t42 * pkin(8);
t11 = t17 - t89;
t18 = t62 * t39 + t63 * t40;
t41 = -t62 * t65 + t63 * t68;
t38 = t41 * pkin(8);
t12 = t38 + t18;
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t2 = t67 * t11 - t64 * t12;
t3 = t64 * t11 + t67 * t12;
t95 = t2 * mrSges(6,1) - t3 * mrSges(6,2);
t49 = (-qJ(4) - pkin(7)) * t65;
t51 = t68 * pkin(7) + t57;
t25 = t63 * t49 - t62 * t51;
t15 = t25 - t89;
t26 = t62 * t49 + t63 * t51;
t16 = t38 + t26;
t7 = t67 * t15 - t64 * t16;
t8 = t64 * t15 + t67 * t16;
t94 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t22 = t67 * t41 - t64 * t42;
t23 = t64 * t41 + t67 * t42;
t9 = -t22 * mrSges(6,1) + t23 * mrSges(6,2);
t93 = 0.2e1 * t9;
t24 = -t41 * mrSges(5,1) + t42 * mrSges(5,2);
t92 = 0.2e1 * t24;
t91 = pkin(3) * t62;
t69 = cos(qJ(2));
t88 = t69 * pkin(1);
t53 = t63 * pkin(3) + pkin(4);
t33 = t67 * t53 - t64 * t91;
t86 = t33 * mrSges(6,1);
t34 = t64 * t53 + t67 * t91;
t85 = t34 * mrSges(6,2);
t84 = Ifges(6,5) * t23 + Ifges(6,6) * t22;
t83 = t65 ^ 2 + t96;
t82 = 2 * mrSges(5,3);
t81 = 2 * mrSges(6,3);
t80 = t63 * t42 * mrSges(5,3);
t56 = -t68 * pkin(3) - pkin(2);
t79 = t83 * t54;
t78 = -t65 * mrSges(4,1) - t68 * mrSges(4,2);
t77 = 0.2e1 * mrSges(4,3) * t83;
t28 = -t41 * pkin(4) + t56;
t76 = Ifges(5,1) * t42 ^ 2 + Ifges(6,1) * t23 ^ 2 + Ifges(4,2) * t96 + Ifges(3,3) + (Ifges(4,1) * t65 + 0.2e1 * Ifges(4,4) * t68) * t65 + (0.2e1 * Ifges(5,4) * t42 + Ifges(5,2) * t41) * t41 + (0.2e1 * Ifges(6,4) * t23 + Ifges(6,2) * t22) * t22;
t75 = (t69 * mrSges(3,1) - t66 * mrSges(3,2)) * pkin(1);
t74 = t24 + t9;
t73 = Ifges(4,5) * t65 + Ifges(5,5) * t42 + Ifges(4,6) * t68 + t84 + (mrSges(5,3) * t91 + Ifges(5,6)) * t41 + (t22 * t34 - t23 * t33) * mrSges(6,3);
t55 = -pkin(2) - t88;
t50 = -t68 * mrSges(4,1) + t65 * mrSges(4,2);
t48 = t56 - t88;
t27 = t28 - t88;
t1 = [t76 + (-t17 * t42 + t18 * t41) * t82 + (-t2 * t23 + t3 * t22) * t81 + m(3) * (t66 ^ 2 + t69 ^ 2) * pkin(1) ^ 2 + m(4) * (t83 * t54 ^ 2 + t55 ^ 2) + t54 * t77 + m(6) * (t2 ^ 2 + t27 ^ 2 + t3 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t48 ^ 2) + 0.2e1 * t55 * t50 + t48 * t92 + t27 * t93 + Ifges(2,3) + 0.2e1 * t75; t75 + t76 + m(4) * (-pkin(2) * t55 + pkin(7) * t79) + (pkin(7) * t83 + t79) * mrSges(4,3) + ((-t17 - t25) * t42 + (t18 + t26) * t41) * mrSges(5,3) + ((-t2 - t7) * t23 + (t3 + t8) * t22) * mrSges(6,3) + m(6) * (t7 * t2 + t28 * t27 + t8 * t3) + m(5) * (t25 * t17 + t26 * t18 + t56 * t48) + (t27 + t28) * t9 + (t55 - pkin(2)) * t50 + (t56 + t48) * t24; -0.2e1 * pkin(2) * t50 + t56 * t92 + t28 * t93 + (t8 * t22 - t7 * t23) * t81 + (-t25 * t42 + t26 * t41) * t82 + pkin(7) * t77 + m(6) * (t28 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2 + t56 ^ 2) + m(4) * (t83 * pkin(7) ^ 2 + pkin(2) ^ 2) + t76; t73 + (-t80 + m(5) * (t17 * t63 + t18 * t62)) * pkin(3) + m(6) * (t33 * t2 + t34 * t3) + t78 * t54 + t17 * mrSges(5,1) - t18 * mrSges(5,2) + t95; t73 + (-t80 + m(5) * (t25 * t63 + t26 * t62)) * pkin(3) + m(6) * (t33 * t7 + t34 * t8) + t78 * pkin(7) + t25 * mrSges(5,1) - t26 * mrSges(5,2) + t94; 0.2e1 * t86 - 0.2e1 * t85 + Ifges(4,3) + Ifges(5,3) + Ifges(6,3) + m(6) * (t33 ^ 2 + t34 ^ 2) + (0.2e1 * t63 * mrSges(5,1) - 0.2e1 * t62 * mrSges(5,2) + m(5) * (t62 ^ 2 + t63 ^ 2) * pkin(3)) * pkin(3); m(5) * t48 + m(6) * t27 + t74; m(5) * t56 + m(6) * t28 + t74; 0; m(5) + m(6); t84 + t95; t84 + t94; Ifges(6,3) - t85 + t86; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
