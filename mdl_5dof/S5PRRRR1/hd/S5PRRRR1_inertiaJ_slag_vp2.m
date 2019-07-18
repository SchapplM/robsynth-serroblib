% Calculate joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaJ_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:21
% EndTime: 2019-07-18 13:28:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (257->96), mult. (660->147), div. (0->0), fcn. (655->8), ass. (0->51)
t43 = sin(qJ(5));
t47 = cos(qJ(5));
t28 = -t47 * mrSges(6,1) + t43 * mrSges(6,2);
t73 = mrSges(5,1) - t28;
t44 = sin(qJ(4));
t48 = cos(qJ(4));
t59 = t43 ^ 2 + t47 ^ 2;
t72 = (t73 * t48 + (t59 * mrSges(6,3) - mrSges(5,2)) * t44) * pkin(2);
t71 = 0.2e1 * pkin(2);
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t27 = t44 * t49 + t48 * t45;
t46 = sin(qJ(2));
t17 = t27 * t46;
t70 = t17 ^ 2;
t38 = t45 ^ 2;
t69 = pkin(2) * t49;
t68 = Ifges(6,4) * t43;
t67 = Ifges(6,4) * t47;
t26 = t44 * t45 - t48 * t49;
t19 = t26 * t46;
t50 = cos(qJ(2));
t12 = -t47 * t19 - t50 * t43;
t66 = t12 * t47;
t65 = t17 * t48;
t64 = t27 * t43;
t63 = t27 * t47;
t51 = pkin(2) ^ 2;
t62 = t44 ^ 2 * t51;
t61 = t43 * mrSges(6,3);
t60 = Ifges(6,5) * t63 + Ifges(6,3) * t26;
t29 = Ifges(6,5) * t43 + Ifges(6,6) * t47;
t41 = t49 ^ 2;
t58 = t38 + t41;
t30 = Ifges(6,2) * t47 + t68;
t31 = Ifges(6,1) * t43 + t67;
t57 = t47 * t30 + t43 * t31 + Ifges(5,3);
t54 = mrSges(6,1) * t43 + mrSges(6,2) * t47;
t3 = Ifges(6,6) * t26 + (-Ifges(6,2) * t43 + t67) * t27;
t4 = Ifges(6,5) * t26 + (Ifges(6,1) * t47 - t68) * t27;
t53 = -t30 * t64 / 0.2e1 + t31 * t63 / 0.2e1 + t43 * t4 / 0.2e1 + t47 * t3 / 0.2e1 + Ifges(5,5) * t27 + (t29 / 0.2e1 - Ifges(5,6)) * t26;
t11 = t43 * t19 - t50 * t47;
t52 = t19 * mrSges(5,2) + mrSges(6,3) * t66 - t11 * t61 - t73 * t17;
t42 = t50 ^ 2;
t39 = t46 ^ 2;
t35 = t48 ^ 2 * t51;
t9 = t26 * mrSges(5,1) + t27 * mrSges(5,2);
t8 = t26 * mrSges(6,1) - mrSges(6,3) * t63;
t7 = -t26 * mrSges(6,2) - t27 * t61;
t5 = t54 * t27;
t1 = [m(2) + m(3) * (t42 + t39) + m(4) * (t58 * t39 + t42) + m(5) * (t19 ^ 2 + t42 + t70) + m(6) * (t11 ^ 2 + t12 ^ 2 + t70); t11 * t8 + t12 * t7 + t17 * t5 + m(6) * (-t11 * t47 - t12 * t43) * t69 + (t17 * t27 + t19 * t26) * mrSges(5,3) + (t58 * mrSges(4,3) - mrSges(3,2)) * t46 + (-t45 * mrSges(4,2) + mrSges(3,1) - t9 + (m(5) * pkin(2) + mrSges(4,1)) * t49) * t50; Ifges(4,1) * t38 + Ifges(3,3) + (Ifges(5,2) * t26 + t60) * t26 + (m(6) * t59 + m(5)) * t51 * t41 + (Ifges(5,1) * t27 - t43 * t3 + t47 * t4 + (-Ifges(6,6) * t43 - (2 * Ifges(5,4))) * t26) * t27 + (0.2e1 * Ifges(4,4) * t45 + Ifges(4,2) * t49 + (-t43 * t7 - t47 * t8 - t9) * t71) * t49; (-t45 * mrSges(4,1) - t49 * mrSges(4,2)) * t46 + (m(5) * (-t19 * t44 - t65) / 0.2e1 + m(6) * (-t65 + (-t11 * t43 + t66) * t44) / 0.2e1) * t71 + t52; Ifges(4,5) * t45 + Ifges(4,6) * t49 + ((-t27 * mrSges(5,3) - t5) * t48 + (-mrSges(5,3) * t26 - t43 * t8 + t47 * t7) * t44) * pkin(2) + t53; Ifges(4,3) + m(6) * (t59 * t62 + t35) + m(5) * (t35 + t62) + 0.2e1 * t72 + t57; t52; t53; t72 + t57; t57; t11 * mrSges(6,1) - t12 * mrSges(6,2); -Ifges(6,6) * t64 + t28 * t69 + t60; -t54 * t44 * pkin(2) + t29; t29; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
