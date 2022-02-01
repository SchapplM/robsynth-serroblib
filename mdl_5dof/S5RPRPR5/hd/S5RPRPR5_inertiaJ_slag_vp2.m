% Calculate joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:54
% EndTime: 2022-01-23 09:24:55
% DurationCPUTime: 0.45s
% Computational Cost: add. (753->144), mult. (1560->219), div. (0->0), fcn. (1568->8), ass. (0->63)
t76 = 2 * qJ(2);
t75 = m(5) * pkin(3);
t52 = sin(pkin(9));
t54 = cos(pkin(9));
t57 = sin(qJ(3));
t59 = cos(qJ(3));
t38 = t52 * t59 + t54 * t57;
t53 = sin(pkin(8));
t31 = t38 * t53;
t37 = -t52 * t57 + t54 * t59;
t32 = t37 * t53;
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t13 = -t58 * t31 - t56 * t32;
t14 = -t56 * t31 + t58 * t32;
t69 = -Ifges(6,5) * t14 - Ifges(6,6) * t13;
t74 = pkin(3) * t52;
t30 = Ifges(5,5) * t32;
t29 = Ifges(5,6) * t31;
t46 = t54 * pkin(3) + pkin(4);
t33 = t58 * t46 - t56 * t74;
t73 = t33 * mrSges(6,1);
t34 = t56 * t46 + t58 * t74;
t72 = t34 * mrSges(6,2);
t71 = t53 * t57;
t70 = t59 * t53;
t55 = cos(pkin(8));
t42 = -pkin(2) * t55 - t53 * pkin(6) - pkin(1);
t36 = t59 * t42;
t67 = qJ(4) * t53;
t68 = qJ(2) * t57;
t18 = -t59 * t67 + t36 + (-pkin(3) - t68) * t55;
t27 = t59 * t55 * qJ(2) + t57 * t42;
t23 = -t57 * t67 + t27;
t7 = t52 * t18 + t54 * t23;
t41 = pkin(3) * t71 + t53 * qJ(2);
t66 = -Ifges(6,3) - Ifges(5,3) - Ifges(4,3);
t65 = t31 * mrSges(5,1) + t32 * mrSges(5,2);
t64 = -t13 * mrSges(6,1) + t14 * mrSges(6,2);
t19 = t58 * t37 - t56 * t38;
t20 = t56 * t37 + t58 * t38;
t63 = t19 * mrSges(6,1) - t20 * mrSges(6,2);
t6 = t54 * t18 - t52 * t23;
t4 = -t55 * pkin(4) - t32 * pkin(7) + t6;
t5 = -t31 * pkin(7) + t7;
t2 = t58 * t4 - t56 * t5;
t3 = t56 * t4 + t58 * t5;
t62 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t69;
t60 = qJ(2) ^ 2;
t51 = t55 ^ 2;
t50 = t53 ^ 2;
t49 = t53 * mrSges(3,2);
t47 = t50 * t60;
t44 = Ifges(4,5) * t70;
t40 = -t55 * mrSges(4,1) - mrSges(4,3) * t70;
t39 = t55 * mrSges(4,2) - mrSges(4,3) * t71;
t26 = -t55 * t68 + t36;
t25 = -t55 * mrSges(5,1) - t32 * mrSges(5,3);
t24 = t55 * mrSges(5,2) - t31 * mrSges(5,3);
t21 = t31 * pkin(4) + t41;
t9 = -t55 * mrSges(6,1) - t14 * mrSges(6,3);
t8 = t55 * mrSges(6,2) + t13 * mrSges(6,3);
t1 = [Ifges(5,2) * t31 ^ 2 + Ifges(6,2) * t13 ^ 2 - 0.2e1 * pkin(1) * t49 + 0.2e1 * t27 * t39 + 0.2e1 * t26 * t40 + 0.2e1 * t41 * t65 + 0.2e1 * t21 * t64 + 0.2e1 * t7 * t24 + 0.2e1 * t6 * t25 + 0.2e1 * t3 * t8 + 0.2e1 * t2 * t9 + Ifges(2,3) + (Ifges(5,1) * t32 - 0.2e1 * Ifges(5,4) * t31) * t32 + (Ifges(6,1) * t14 + 0.2e1 * Ifges(6,4) * t13) * t14 + (t50 + t51) * mrSges(3,3) * t76 + m(3) * (pkin(1) ^ 2 + t51 * t60 + t47) + m(4) * (t26 ^ 2 + t27 ^ 2 + t47) + m(5) * (t41 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(6) * (t2 ^ 2 + t21 ^ 2 + t3 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * t30 + 0.2e1 * t29 - t44 + (Ifges(3,2) - t66) * t55 + 0.2e1 * t69) * t55 + ((Ifges(3,1) + (mrSges(4,2) * t76 + Ifges(4,1) * t59) * t59 + (mrSges(4,1) * t76 - 0.2e1 * Ifges(4,4) * t59 + Ifges(4,2) * t57) * t57) * t53 + (-t59 * Ifges(4,5) + 0.2e1 * Ifges(4,6) * t57 + (2 * Ifges(3,4))) * t55) * t53; -m(3) * pkin(1) - t55 * mrSges(3,1) + t19 * t9 + t20 * t8 + t38 * t24 + t37 * t25 + t57 * t39 + t59 * t40 + t49 + m(6) * (t19 * t2 + t20 * t3) + m(5) * (t37 * t6 + t38 * t7) + m(4) * (t59 * t26 + t57 * t27); m(3) + m(4) * (t57 ^ 2 + t59 ^ 2) + m(5) * (t37 ^ 2 + t38 ^ 2) + m(6) * (t19 ^ 2 + t20 ^ 2); t44 - Ifges(4,6) * t71 + t30 - t29 + m(6) * (t33 * t2 + t34 * t3) + t34 * t8 + t33 * t9 + t6 * mrSges(5,1) - t7 * mrSges(5,2) - t27 * mrSges(4,2) + t26 * mrSges(4,1) + t66 * t55 + (t52 * t24 + t54 * t25 + m(5) * (t52 * t7 + t54 * t6)) * pkin(3) + t62; t59 * mrSges(4,1) + t37 * mrSges(5,1) - t57 * mrSges(4,2) - t38 * mrSges(5,2) + m(6) * (t33 * t19 + t34 * t20) + (t37 * t54 + t38 * t52) * t75 + t63; 0.2e1 * t73 - 0.2e1 * t72 + m(6) * (t33 ^ 2 + t34 ^ 2) - t66 + (0.2e1 * t54 * mrSges(5,1) - 0.2e1 * t52 * mrSges(5,2) + (t52 ^ 2 + t54 ^ 2) * t75) * pkin(3); m(5) * t41 + m(6) * t21 + t64 + t65; 0; 0; m(5) + m(6); -Ifges(6,3) * t55 + t62; t63; Ifges(6,3) - t72 + t73; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
