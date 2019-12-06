% Calculate joint inertia matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:03
% EndTime: 2019-12-05 16:35:07
% DurationCPUTime: 0.64s
% Computational Cost: add. (544->191), mult. (1263->283), div. (0->0), fcn. (1203->10), ass. (0->85)
t96 = 2 * pkin(7);
t95 = 2 * qJ(4);
t64 = cos(pkin(5));
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t62 = sin(pkin(5));
t67 = sin(qJ(2));
t80 = t62 * t67;
t30 = t64 * t66 + t69 * t80;
t61 = sin(pkin(10));
t63 = cos(pkin(10));
t70 = cos(qJ(2));
t79 = t62 * t70;
t9 = t30 * t61 + t63 * t79;
t94 = t9 ^ 2;
t28 = -t64 * t69 + t66 * t80;
t93 = t28 ^ 2;
t60 = t69 ^ 2;
t92 = m(5) * pkin(3);
t91 = pkin(7) * t69;
t90 = t9 * t61;
t77 = t63 * t66;
t39 = -t69 * mrSges(5,1) - mrSges(5,3) * t77;
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t31 = -t65 * t77 - t68 * t69;
t32 = -t65 * t69 + t68 * t77;
t8 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t89 = -t39 + t8;
t88 = Ifges(5,4) * t61;
t87 = Ifges(5,4) * t63;
t86 = Ifges(5,5) * t69;
t85 = Ifges(5,6) * t69;
t84 = t28 * t66;
t83 = t61 * t65;
t82 = t61 * t66;
t81 = t61 * t68;
t41 = -t69 * pkin(3) - t66 * qJ(4) - pkin(2);
t78 = t63 * t41;
t42 = -t63 * mrSges(5,1) + t61 * mrSges(5,2);
t76 = -mrSges(4,1) + t42;
t33 = mrSges(5,1) * t82 + mrSges(5,2) * t77;
t19 = t61 * t41 + t63 * t91;
t75 = qJ(4) * t63;
t5 = Ifges(6,5) * t32 + Ifges(6,6) * t31 + Ifges(6,3) * t82;
t11 = t30 * t63 - t61 * t79;
t74 = t11 * t63 + t90;
t73 = t30 * t69 + t84;
t20 = Ifges(6,5) * t81 - Ifges(6,6) * t83 - Ifges(6,3) * t63;
t72 = pkin(7) ^ 2;
t71 = qJ(4) ^ 2;
t59 = t66 ^ 2;
t58 = t63 ^ 2;
t57 = t62 ^ 2;
t56 = t61 ^ 2;
t55 = t59 * t72;
t53 = t56 * t71;
t51 = t57 * t70 ^ 2;
t45 = -t69 * mrSges(4,1) + t66 * mrSges(4,2);
t44 = Ifges(5,1) * t61 + t87;
t43 = Ifges(5,2) * t63 + t88;
t40 = -t63 * pkin(4) - t61 * pkin(8) - pkin(3);
t38 = t69 * mrSges(5,2) - mrSges(5,3) * t82;
t37 = -t63 * mrSges(6,1) - mrSges(6,3) * t81;
t36 = t63 * mrSges(6,2) - mrSges(6,3) * t83;
t34 = (mrSges(6,1) * t65 + mrSges(6,2) * t68) * t61;
t25 = (pkin(4) * t61 - pkin(8) * t63 + pkin(7)) * t66;
t24 = -t86 + (Ifges(5,1) * t63 - t88) * t66;
t23 = -t85 + (-Ifges(5,2) * t61 + t87) * t66;
t22 = -Ifges(6,5) * t63 + (Ifges(6,1) * t68 - Ifges(6,4) * t65) * t61;
t21 = -Ifges(6,6) * t63 + (Ifges(6,4) * t68 - Ifges(6,2) * t65) * t61;
t18 = -t61 * t91 + t78;
t17 = t65 * t40 + t68 * t75;
t16 = t68 * t40 - t65 * t75;
t15 = mrSges(6,1) * t82 - t32 * mrSges(6,3);
t14 = -mrSges(6,2) * t82 + t31 * mrSges(6,3);
t13 = -t69 * pkin(8) + t19;
t12 = -t78 + (pkin(7) * t61 + pkin(4)) * t69;
t7 = Ifges(6,1) * t32 + Ifges(6,4) * t31 + Ifges(6,5) * t82;
t6 = Ifges(6,4) * t32 + Ifges(6,2) * t31 + Ifges(6,6) * t82;
t4 = t68 * t13 + t65 * t25;
t3 = -t65 * t13 + t68 * t25;
t2 = t11 * t68 + t28 * t65;
t1 = -t11 * t65 + t28 * t68;
t10 = [m(2) + m(5) * (t11 ^ 2 + t93 + t94) + m(6) * (t1 ^ 2 + t2 ^ 2 + t94) + m(4) * (t30 ^ 2 + t51 + t93) + m(3) * (t57 * t67 ^ 2 + t64 ^ 2 + t51); t1 * t15 + t11 * t38 + t2 * t14 + t28 * t33 + t89 * t9 + t73 * mrSges(4,3) + (-t67 * mrSges(3,2) + (mrSges(3,1) - t45) * t70) * t62 + m(5) * (pkin(7) * t84 + t19 * t11 - t18 * t9) + m(6) * (t3 * t1 + t12 * t9 + t4 * t2) + m(4) * (pkin(2) * t79 + pkin(7) * t73); -0.2e1 * pkin(2) * t45 + 0.2e1 * t12 * t8 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t18 * t39 + 0.2e1 * t19 * t38 + t31 * t6 + t32 * t7 + Ifges(3,3) + (t59 + t60) * mrSges(4,3) * t96 + (Ifges(4,2) + Ifges(5,3)) * t60 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2 + t55) + m(4) * (pkin(2) ^ 2 + t60 * t72 + t55) + (Ifges(4,1) * t66 + 0.2e1 * Ifges(4,4) * t69 + t33 * t96 + (t24 - t86) * t63 + (-t23 + t5 + t85) * t61) * t66; -t30 * mrSges(4,2) + t1 * t37 + t2 * t36 + t9 * t34 + t76 * t28 + t74 * mrSges(5,3) + m(5) * (-pkin(3) * t28 + qJ(4) * t74) + m(6) * (qJ(4) * t90 + t16 * t1 + t17 * t2); t4 * t36 + t3 * t37 + t31 * t21 / 0.2e1 + t32 * t22 / 0.2e1 - pkin(3) * t33 + t12 * t34 + t16 * t15 + t17 * t14 + m(6) * (t16 * t3 + t17 * t4) + (-pkin(7) * mrSges(4,2) + Ifges(4,6)) * t69 + (-t85 / 0.2e1 + t23 / 0.2e1 - t5 / 0.2e1 + t19 * mrSges(5,3)) * t63 + (-t86 / 0.2e1 + t24 / 0.2e1 - t18 * mrSges(5,3) - t65 * t6 / 0.2e1 + t68 * t7 / 0.2e1) * t61 + ((m(5) * t19 + t38) * t63 + (-m(5) * t18 + m(6) * t12 + t89) * t61) * qJ(4) + (Ifges(4,5) + t63 * t44 / 0.2e1 + (t20 / 0.2e1 - t43 / 0.2e1) * t61 + (t76 - t92) * pkin(7)) * t66; -0.2e1 * pkin(3) * t42 + 0.2e1 * t16 * t37 + 0.2e1 * t17 * t36 + Ifges(4,3) + (-t20 + t43) * t63 + (t56 + t58) * mrSges(5,3) * t95 + m(6) * (t16 ^ 2 + t17 ^ 2 + t53) + m(5) * (pkin(3) ^ 2 + t58 * t71 + t53) + (-t65 * t21 + t68 * t22 + t34 * t95 + t44) * t61; m(5) * t28 + m(6) * (t68 * t1 + t65 * t2); m(6) * (t68 * t3 + t65 * t4) + t65 * t14 + t68 * t15 + m(5) * t66 * pkin(7) + t33; m(6) * (t68 * t16 + t65 * t17) + t65 * t36 + t68 * t37 - t92 + t42; m(5) + m(6) * (t65 ^ 2 + t68 ^ 2); t1 * mrSges(6,1) - t2 * mrSges(6,2); t3 * mrSges(6,1) - t4 * mrSges(6,2) + t5; t16 * mrSges(6,1) - t17 * mrSges(6,2) + t20; t68 * mrSges(6,1) - t65 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
