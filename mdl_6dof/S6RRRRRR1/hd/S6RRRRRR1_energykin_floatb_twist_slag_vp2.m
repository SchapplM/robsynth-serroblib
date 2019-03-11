% Calculate kinetic energy for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:46
% EndTime: 2019-03-10 03:27:48
% DurationCPUTime: 1.44s
% Computational Cost: add. (4691->156), mult. (6231->230), div. (0->0), fcn. (5148->12), ass. (0->61)
t68 = sin(qJ(1));
t74 = cos(qJ(1));
t54 = -t68 * V_base(4) + t74 * V_base(5);
t55 = t68 * V_base(5) + t74 * V_base(4);
t42 = -pkin(1) * t54 - pkin(7) * t55 + V_base(3);
t59 = V_base(5) * pkin(6) + V_base(1);
t60 = -V_base(4) * pkin(6) + V_base(2);
t49 = t74 * t59 + t68 * t60;
t62 = V_base(6) + qJD(1);
t45 = pkin(7) * t62 + t49;
t67 = sin(qJ(2));
t73 = cos(qJ(2));
t35 = t73 * t42 - t45 * t67;
t47 = t55 * t73 + t62 * t67;
t53 = qJD(2) - t54;
t32 = pkin(2) * t53 - pkin(8) * t47 + t35;
t36 = t67 * t42 + t73 * t45;
t46 = -t55 * t67 + t62 * t73;
t34 = pkin(8) * t46 + t36;
t66 = sin(qJ(3));
t72 = cos(qJ(3));
t25 = t72 * t32 - t34 * t66;
t38 = t46 * t66 + t47 * t72;
t52 = qJD(3) + t53;
t16 = pkin(3) * t52 - pkin(9) * t38 + t25;
t26 = t66 * t32 + t72 * t34;
t37 = t46 * t72 - t47 * t66;
t20 = pkin(9) * t37 + t26;
t65 = sin(qJ(4));
t71 = cos(qJ(4));
t13 = t65 * t16 + t71 * t20;
t27 = t37 * t71 - t38 * t65;
t11 = pkin(10) * t27 + t13;
t64 = sin(qJ(5));
t70 = cos(qJ(5));
t12 = t71 * t16 - t20 * t65;
t28 = t37 * t65 + t38 * t71;
t51 = qJD(4) + t52;
t9 = pkin(4) * t51 - pkin(10) * t28 + t12;
t6 = t70 * t11 + t64 * t9;
t48 = -t68 * t59 + t60 * t74;
t5 = -t11 * t64 + t70 * t9;
t22 = t27 * t70 - t28 * t64;
t44 = -pkin(1) * t62 - t48;
t39 = -pkin(2) * t46 + t44;
t29 = -pkin(3) * t37 + t39;
t24 = -pkin(4) * t27 + t29;
t75 = V_base(3) ^ 2;
t69 = cos(qJ(6));
t63 = sin(qJ(6));
t50 = qJD(5) + t51;
t23 = t27 * t64 + t28 * t70;
t21 = qJD(6) - t22;
t18 = t23 * t69 + t50 * t63;
t17 = -t23 * t63 + t50 * t69;
t7 = -pkin(5) * t22 - pkin(11) * t23 + t24;
t4 = pkin(11) * t50 + t6;
t3 = -pkin(5) * t50 - t5;
t2 = t4 * t69 + t63 * t7;
t1 = -t4 * t63 + t69 * t7;
t8 = (t39 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(4,5) * t52 + Ifges(4,1) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(2,2) - t48 * mrSges(2,3) + Ifges(2,5) * t62 + Ifges(2,1) * t55 / 0.2e1) * t55 + (-V_base(3) * mrSges(2,1) + t49 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t62 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t24 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t50 + Ifges(6,1) * t23 / 0.2e1) * t23 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t51 / 0.2e1) * t51 + (t48 * mrSges(2,1) - t49 * mrSges(2,2) + Ifges(2,3) * t62 / 0.2e1) * t62 + (-t44 * mrSges(3,1) + t36 * mrSges(3,3) + Ifges(3,4) * t47 + Ifges(3,6) * t53 + Ifges(3,2) * t46 / 0.2e1) * t46 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t21 / 0.2e1) * t21 + (-t39 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t38 + Ifges(4,6) * t52 + Ifges(4,2) * t37 / 0.2e1) * t37 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t50 / 0.2e1) * t50 + (t25 * mrSges(4,1) - t26 * mrSges(4,2) + Ifges(4,3) * t52 / 0.2e1) * t52 + m(2) * (t48 ^ 2 + t49 ^ 2 + t75) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t75) / 0.2e1 + m(4) * (t25 ^ 2 + t26 ^ 2 + t39 ^ 2) / 0.2e1 + m(3) * (t35 ^ 2 + t36 ^ 2 + t44 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t29 ^ 2) / 0.2e1 + (t44 * mrSges(3,2) - t35 * mrSges(3,3) + Ifges(3,5) * t53 + Ifges(3,1) * t47 / 0.2e1) * t47 + m(6) * (t24 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t29 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t28 + Ifges(5,6) * t51 + Ifges(5,2) * t27 / 0.2e1) * t27 + (-t24 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t23 + Ifges(6,6) * t50 + Ifges(6,2) * t22 / 0.2e1) * t22 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t18 + Ifges(7,6) * t21 + Ifges(7,2) * t17 / 0.2e1) * t17 + (t29 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t51 + Ifges(5,1) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t35 * mrSges(3,1) - t36 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t21 + Ifges(7,1) * t18 / 0.2e1) * t18;
T  = t8;
