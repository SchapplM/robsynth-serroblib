% Calculate joint inertia matrix for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaJ_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:55
% EndTime: 2019-07-18 13:24:56
% DurationCPUTime: 0.38s
% Computational Cost: add. (223->131), mult. (593->191), div. (0->0), fcn. (516->6), ass. (0->61)
t70 = 2 * qJ(2);
t42 = sin(qJ(4));
t44 = cos(qJ(5));
t45 = cos(qJ(4));
t41 = sin(qJ(5));
t60 = Ifges(6,4) * t41;
t10 = -Ifges(6,5) * t45 + (Ifges(6,1) * t44 - t60) * t42;
t69 = t10 / 0.2e1;
t46 = cos(qJ(3));
t54 = t46 * t44;
t43 = sin(qJ(3));
t56 = t43 * t45;
t16 = -t41 * t56 - t54;
t68 = t16 / 0.2e1;
t59 = Ifges(6,4) * t44;
t27 = Ifges(6,1) * t41 + t59;
t67 = t27 / 0.2e1;
t66 = t44 / 0.2e1;
t55 = t46 * t41;
t17 = t44 * t56 - t55;
t65 = -t46 * mrSges(5,1) + t16 * mrSges(6,1) - t17 * mrSges(6,2) - mrSges(5,3) * t56;
t64 = mrSges(5,2) * t45;
t63 = mrSges(6,3) * t42;
t62 = Ifges(5,4) * t42;
t61 = Ifges(5,4) * t45;
t40 = t46 ^ 2;
t47 = qJ(2) ^ 2;
t58 = t40 * t47;
t57 = t42 * t43;
t53 = t44 * mrSges(6,1) - t41 * mrSges(6,2) + mrSges(5,1);
t52 = Ifges(5,5) * t42 + Ifges(5,6) * t45;
t51 = t41 ^ 2 + t44 ^ 2;
t36 = t42 ^ 2;
t39 = t45 ^ 2;
t50 = t36 + t39;
t49 = qJ(2) * t46;
t1 = Ifges(6,5) * t17 + Ifges(6,6) * t16 + Ifges(6,3) * t57;
t12 = (t43 * t44 - t45 * t55) * qJ(2);
t13 = (t41 * t43 + t45 * t54) * qJ(2);
t48 = -t12 * t41 + t13 * t44;
t37 = t43 ^ 2;
t32 = t37 * t47;
t31 = Ifges(5,5) * t56;
t29 = t36 * t58;
t28 = Ifges(5,1) * t42 + t61;
t26 = Ifges(5,2) * t45 + t62;
t25 = Ifges(6,2) * t44 + t60;
t24 = Ifges(6,5) * t41 + Ifges(6,6) * t44;
t21 = -t45 * mrSges(6,1) - t44 * t63;
t20 = t46 * mrSges(5,2) - mrSges(5,3) * t57;
t19 = t45 * mrSges(6,2) - t41 * t63;
t18 = (mrSges(6,1) * t41 + mrSges(6,2) * t44) * t42;
t11 = -Ifges(5,5) * t46 + (Ifges(5,1) * t45 - t62) * t43;
t9 = -Ifges(5,6) * t46 + (-Ifges(5,2) * t42 + t61) * t43;
t8 = -Ifges(6,6) * t45 + (-Ifges(6,2) * t41 + t59) * t42;
t7 = -Ifges(6,3) * t45 + (Ifges(6,5) * t44 - Ifges(6,6) * t41) * t42;
t6 = mrSges(6,1) * t57 - t17 * mrSges(6,3);
t5 = -mrSges(6,2) * t57 + t16 * mrSges(6,3);
t3 = Ifges(6,1) * t17 + Ifges(6,4) * t16 + Ifges(6,5) * t57;
t2 = Ifges(6,4) * t17 + Ifges(6,2) * t16 + Ifges(6,6) * t57;
t4 = [(m(3) * t47) + 0.2e1 * t12 * t6 + 0.2e1 * t13 * t5 + t16 * t2 + t17 * t3 + Ifges(3,2) + Ifges(2,3) + (-t31 + (Ifges(5,3) + Ifges(4,2)) * t46) * t46 + m(6) * (t12 ^ 2 + t13 ^ 2 + t29) + m(5) * (t39 * t58 + t29 + t32) + m(4) * (t32 + t58) + (mrSges(3,3) + (t37 + t40) * mrSges(4,3) + (t20 * t45 - t42 * t65) * t46) * t70 + (t42 * t1 + t45 * t11 - t42 * t9 + (Ifges(4,1) + (mrSges(5,1) * t42 + t64) * t70) * t43 + (Ifges(5,6) * t42 + (2 * Ifges(4,4))) * t46) * t43; -t46 * mrSges(4,1) + t43 * mrSges(4,2) - mrSges(3,1) + t65 * t45 + (m(6) * (-t45 * t49 + t48) + t44 * t5 - t41 * t6 + t20) * t42; m(3) + m(4) + m(5) * t50 + m(6) * (t51 * t36 + t39); t12 * t21 + t13 * t19 + t17 * t69 + t8 * t68 + (-t1 / 0.2e1 + t9 / 0.2e1) * t45 + (t3 * t66 - t41 * t2 / 0.2e1 + t11 / 0.2e1) * t42 + (Ifges(4,5) + t45 * t28 / 0.2e1 + (t7 / 0.2e1 - t26 / 0.2e1) * t42 + (-t45 * mrSges(5,1) + t42 * mrSges(5,2) - mrSges(4,1)) * qJ(2)) * t43 + (-t52 / 0.2e1 + Ifges(4,6) + (t50 * mrSges(5,3) + t42 * t18 - mrSges(4,2)) * qJ(2)) * t46; -t45 * t18 + (t19 * t44 - t21 * t41) * t42; Ifges(4,3) + (t26 - t7) * t45 + (t10 * t44 - t41 * t8 + t28) * t42; t31 + t17 * t67 + t25 * t68 + t41 * t3 / 0.2e1 + t2 * t66 + (-qJ(2) * t64 - Ifges(5,3)) * t46 + t48 * mrSges(6,3) + ((t24 / 0.2e1 - Ifges(5,6)) * t43 - t53 * t49) * t42; t53 * t45 + (t51 * mrSges(6,3) - mrSges(5,2)) * t42; -t45 * t24 / 0.2e1 + (t42 * t67 + t8 / 0.2e1) * t44 + (-t42 * t25 / 0.2e1 + t69) * t41 + t52; t44 * t25 + t41 * t27 + Ifges(5,3); t12 * mrSges(6,1) - t13 * mrSges(6,2) + t1; -t18; t7; t24; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq  = res;
