% Calculate joint inertia matrix for
% S5PRRPR6
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:09
% EndTime: 2019-12-05 16:30:12
% DurationCPUTime: 0.56s
% Computational Cost: add. (553->178), mult. (1269->267), div. (0->0), fcn. (1242->10), ass. (0->81)
t94 = 2 * pkin(7);
t65 = cos(pkin(10));
t68 = sin(qJ(3));
t81 = t65 * t68;
t63 = sin(pkin(10));
t84 = t63 * t68;
t32 = mrSges(5,1) * t84 + mrSges(5,2) * t81;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t39 = t70 * t63 + t67 * t65;
t26 = t39 * t68;
t38 = -t67 * t63 + t70 * t65;
t27 = t38 * t68;
t7 = t26 * mrSges(6,1) + t27 * mrSges(6,2);
t93 = t32 + t7;
t66 = cos(pkin(5));
t71 = cos(qJ(3));
t64 = sin(pkin(5));
t69 = sin(qJ(2));
t83 = t64 * t69;
t29 = -t66 * t71 + t68 * t83;
t28 = t29 ^ 2;
t92 = m(5) * pkin(3);
t91 = -t63 / 0.2e1;
t90 = t65 / 0.2e1;
t89 = pkin(7) * t71;
t88 = Ifges(5,4) * t63;
t87 = Ifges(5,4) * t65;
t86 = t29 * t68;
t31 = t66 * t68 + t71 * t83;
t85 = t31 * t71;
t72 = cos(qJ(2));
t82 = t64 * t72;
t80 = pkin(8) + qJ(4);
t79 = -Ifges(6,5) * t27 + Ifges(6,6) * t26;
t45 = -t65 * mrSges(5,1) + t63 * mrSges(5,2);
t78 = t45 - mrSges(4,1);
t43 = -t71 * pkin(3) - t68 * qJ(4) - pkin(2);
t20 = t63 * t43 + t65 * t89;
t77 = t63 ^ 2 + t65 ^ 2;
t9 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t12 = -t31 * t63 - t65 * t82;
t13 = t31 * t65 - t63 * t82;
t76 = -t12 * t63 + t13 * t65;
t37 = t65 * t43;
t19 = -t63 * t89 + t37;
t75 = -t19 * t63 + t20 * t65;
t74 = pkin(7) ^ 2;
t62 = t71 ^ 2;
t61 = t68 ^ 2;
t59 = t64 ^ 2;
t57 = t61 * t74;
t55 = -t65 * pkin(4) - pkin(3);
t53 = t59 * t72 ^ 2;
t49 = -t71 * mrSges(4,1) + t68 * mrSges(4,2);
t48 = Ifges(5,1) * t63 + t87;
t47 = Ifges(5,2) * t65 + t88;
t46 = t80 * t65;
t44 = t80 * t63;
t42 = (pkin(4) * t63 + pkin(7)) * t68;
t41 = -t71 * mrSges(5,1) - mrSges(5,3) * t81;
t40 = t71 * mrSges(5,2) - mrSges(5,3) * t84;
t35 = Ifges(6,5) * t39;
t34 = Ifges(6,6) * t38;
t25 = -Ifges(5,5) * t71 + (Ifges(5,1) * t65 - t88) * t68;
t24 = -Ifges(5,6) * t71 + (-Ifges(5,2) * t63 + t87) * t68;
t18 = -t71 * mrSges(6,1) - t27 * mrSges(6,3);
t17 = t71 * mrSges(6,2) - t26 * mrSges(6,3);
t16 = -t67 * t44 + t70 * t46;
t15 = -t70 * t44 - t67 * t46;
t14 = -pkin(8) * t84 + t20;
t11 = Ifges(6,1) * t39 + Ifges(6,4) * t38;
t10 = Ifges(6,4) * t39 + Ifges(6,2) * t38;
t8 = -pkin(8) * t81 + t37 + (-pkin(7) * t63 - pkin(4)) * t71;
t6 = Ifges(6,1) * t27 - Ifges(6,4) * t26 - Ifges(6,5) * t71;
t5 = Ifges(6,4) * t27 - Ifges(6,2) * t26 - Ifges(6,6) * t71;
t4 = t67 * t12 + t70 * t13;
t3 = t70 * t12 - t67 * t13;
t2 = t70 * t14 + t67 * t8;
t1 = -t67 * t14 + t70 * t8;
t21 = [m(2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t28) + m(6) * (t3 ^ 2 + t4 ^ 2 + t28) + m(4) * (t31 ^ 2 + t28 + t53) + m(3) * (t59 * t69 ^ 2 + t66 ^ 2 + t53); mrSges(4,3) * t85 + t12 * t41 + t13 * t40 + t4 * t17 + t3 * t18 + (-t69 * mrSges(3,2) + (mrSges(3,1) - t49) * t72) * t64 + (t68 * mrSges(4,3) + t93) * t29 + m(5) * (pkin(7) * t86 + t19 * t12 + t20 * t13) + m(6) * (t1 * t3 + t2 * t4 + t42 * t29) + m(4) * (pkin(2) * t82 + (t85 + t86) * pkin(7)); -0.2e1 * pkin(2) * t49 + 0.2e1 * t1 * t18 + 0.2e1 * t2 * t17 + 0.2e1 * t19 * t41 + 0.2e1 * t20 * t40 - t26 * t5 + t27 * t6 + 0.2e1 * t42 * t7 + Ifges(3,3) + (t61 + t62) * mrSges(4,3) * t94 + m(6) * (t1 ^ 2 + t2 ^ 2 + t42 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t57) + m(4) * (pkin(2) ^ 2 + t62 * t74 + t57) + ((Ifges(6,3) + Ifges(5,3) + Ifges(4,2)) * t71 + t79) * t71 + (Ifges(4,1) * t68 + t32 * t94 - t63 * t24 + t65 * t25 + (-Ifges(5,5) * t65 + Ifges(5,6) * t63 + (2 * Ifges(4,4))) * t71) * t68; -t31 * mrSges(4,2) + (-t3 * t39 + t4 * t38) * mrSges(6,3) + t76 * mrSges(5,3) + (t9 + t78) * t29 + m(5) * (-pkin(3) * t29 + t76 * qJ(4)) + m(6) * (t15 * t3 + t16 * t4 + t55 * t29); t63 * t25 / 0.2e1 + t24 * t90 + t55 * t7 - pkin(3) * t32 + t38 * t5 / 0.2e1 + t39 * t6 / 0.2e1 + t42 * t9 + t16 * t17 + t15 * t18 - t26 * t10 / 0.2e1 + t27 * t11 / 0.2e1 + m(6) * (t15 * t1 + t16 * t2 + t55 * t42) + (-pkin(7) * mrSges(4,2) + Ifges(5,5) * t91 - Ifges(5,6) * t65 / 0.2e1 - t35 / 0.2e1 - t34 / 0.2e1 + Ifges(4,6)) * t71 + (-t1 * t39 + t2 * t38) * mrSges(6,3) + t75 * mrSges(5,3) + (m(5) * t75 + t65 * t40 - t63 * t41) * qJ(4) + (t47 * t91 + t48 * t90 + Ifges(4,5) + (t78 - t92) * pkin(7)) * t68; -0.2e1 * pkin(3) * t45 + t38 * t10 + t39 * t11 + t65 * t47 + t63 * t48 + 0.2e1 * t55 * t9 + Ifges(4,3) + m(6) * (t15 ^ 2 + t16 ^ 2 + t55 ^ 2) + m(5) * (t77 * qJ(4) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t15 * t39 + t16 * t38) * mrSges(6,3) + 0.2e1 * t77 * qJ(4) * mrSges(5,3); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t29; m(5) * t68 * pkin(7) + m(6) * t42 + t93; m(6) * t55 + t45 + t9 - t92; m(5) + m(6); t3 * mrSges(6,1) - t4 * mrSges(6,2); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,3) * t71 - t79; t15 * mrSges(6,1) - t16 * mrSges(6,2) + t34 + t35; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
