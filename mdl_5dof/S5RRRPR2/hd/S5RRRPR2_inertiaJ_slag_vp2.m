% Calculate joint inertia matrix for
% S5RRRPR2
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:23
% EndTime: 2020-01-03 12:07:24
% DurationCPUTime: 0.28s
% Computational Cost: add. (490->102), mult. (872->143), div. (0->0), fcn. (636->8), ass. (0->62)
t50 = cos(qJ(5));
t44 = t50 ^ 2;
t47 = sin(qJ(5));
t30 = -t50 * mrSges(6,1) + t47 * mrSges(6,2);
t46 = cos(pkin(9));
t73 = t46 * pkin(3);
t36 = -pkin(4) - t73;
t25 = t36 * t30;
t45 = sin(pkin(9));
t35 = t45 * pkin(3) + pkin(8);
t43 = t47 ^ 2;
t69 = mrSges(6,3) * t43;
t28 = t35 * t69;
t68 = mrSges(6,3) * t44;
t29 = t35 * t68;
t39 = mrSges(5,1) * t73;
t76 = t25 + t28 + t29 + t39;
t52 = cos(qJ(2));
t38 = t52 * pkin(1) + pkin(2);
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t49 = sin(qJ(2));
t74 = pkin(1) * t49;
t23 = t51 * t38 - t48 * t74;
t20 = pkin(3) + t23;
t24 = t48 * t38 + t51 * t74;
t8 = t46 * t20 - t45 * t24;
t6 = -pkin(4) - t8;
t1 = t6 * t30;
t17 = t23 * mrSges(4,1);
t9 = t45 * t20 + t46 * t24;
t7 = pkin(8) + t9;
t2 = t7 * t69;
t3 = t7 * t68;
t4 = t8 * mrSges(5,1);
t75 = t1 + t17 + t2 + t3 + t4;
t72 = t48 * pkin(2);
t71 = t51 * pkin(2);
t70 = t9 * mrSges(5,2);
t37 = pkin(3) + t71;
t22 = t45 * t37 + t46 * t72;
t67 = t22 * mrSges(5,2);
t66 = t24 * mrSges(4,2);
t65 = t45 * mrSges(5,2);
t64 = Ifges(6,5) * t47 + Ifges(6,6) * t50;
t63 = t43 + t44;
t62 = mrSges(4,2) * t72;
t61 = t63 * t35;
t60 = Ifges(6,2) * t44 + Ifges(4,3) + Ifges(5,3) + (Ifges(6,1) * t47 + 0.2e1 * Ifges(6,4) * t50) * t47;
t59 = Ifges(3,3) + t60;
t58 = -mrSges(6,1) * t47 - mrSges(6,2) * t50;
t21 = t46 * t37 - t45 * t72;
t57 = (t52 * mrSges(3,1) - t49 * mrSges(3,2)) * pkin(1);
t18 = -pkin(4) - t21;
t10 = t18 * t30;
t19 = pkin(8) + t22;
t11 = t19 * t69;
t12 = t19 * t68;
t15 = t21 * mrSges(5,1);
t40 = mrSges(4,1) * t71;
t56 = t10 + t11 + t12 + t15 + t40 + t60;
t5 = [-0.2e1 * t66 - 0.2e1 * t70 + Ifges(2,3) + 0.2e1 * t1 + 0.2e1 * t17 + 0.2e1 * t2 + 0.2e1 * t3 + 0.2e1 * t4 + 0.2e1 * t57 + m(6) * (t63 * t7 ^ 2 + t6 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(3) * (t49 ^ 2 + t52 ^ 2) * pkin(1) ^ 2 + t59; t56 + m(6) * (t63 * t7 * t19 + t18 * t6) + m(5) * (t21 * t8 + t22 * t9) + t57 + (-t9 - t22) * mrSges(5,2) + (-t24 - t72) * mrSges(4,2) + m(4) * (t23 * t51 + t24 * t48) * pkin(2) + Ifges(3,3) + t75; -0.2e1 * t62 - 0.2e1 * t67 + 0.2e1 * t10 + 0.2e1 * t11 + 0.2e1 * t12 + 0.2e1 * t15 + 0.2e1 * t40 + m(5) * (t21 ^ 2 + t22 ^ 2) + m(6) * (t63 * t19 ^ 2 + t18 ^ 2) + m(4) * (t48 ^ 2 + t51 ^ 2) * pkin(2) ^ 2 + t59; m(6) * (t36 * t6 + t7 * t61) - t66 - t70 + (m(5) * (t45 * t9 + t46 * t8) - t65) * pkin(3) + t60 + t75 + t76; t56 + m(6) * (t36 * t18 + t19 * t61) + (-t65 + m(5) * (t21 * t46 + t22 * t45)) * pkin(3) - t67 - t62 + t76; 0.2e1 * t25 + 0.2e1 * t28 + 0.2e1 * t29 + 0.2e1 * t39 + m(6) * (t63 * t35 ^ 2 + t36 ^ 2) + t60 + (-0.2e1 * t65 + m(5) * (t45 ^ 2 + t46 ^ 2) * pkin(3)) * pkin(3); 0; 0; 0; m(6) * t63 + m(5); t58 * t7 + t64; t58 * t19 + t64; t58 * t35 + t64; -t30; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
