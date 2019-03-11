% Calculate kinetic energy for
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:53:55
% EndTime: 2019-03-09 06:53:57
% DurationCPUTime: 1.29s
% Computational Cost: add. (3971->156), mult. (5815->228), div. (0->0), fcn. (4792->12), ass. (0->60)
t59 = V_base(5) * pkin(6) + V_base(1);
t60 = -V_base(4) * pkin(6) + V_base(2);
t68 = sin(qJ(1));
t73 = cos(qJ(1));
t50 = -t59 * t68 + t73 * t60;
t55 = t68 * V_base(5) + t73 * V_base(4);
t61 = V_base(6) + qJD(1);
t41 = pkin(1) * t61 - qJ(2) * t55 + t50;
t51 = t73 * t59 + t68 * t60;
t54 = -t68 * V_base(4) + t73 * V_base(5);
t46 = qJ(2) * t54 + t51;
t62 = sin(pkin(11));
t63 = cos(pkin(11));
t37 = t62 * t41 + t63 * t46;
t32 = pkin(7) * t61 + t37;
t48 = t54 * t63 - t62 * t55;
t49 = t54 * t62 + t55 * t63;
t52 = -pkin(1) * t54 + qJD(2) + V_base(3);
t35 = -pkin(2) * t48 - pkin(7) * t49 + t52;
t67 = sin(qJ(3));
t72 = cos(qJ(3));
t25 = -t32 * t67 + t72 * t35;
t40 = t49 * t72 + t61 * t67;
t47 = qJD(3) - t48;
t16 = pkin(3) * t47 - pkin(8) * t40 + t25;
t26 = t72 * t32 + t67 * t35;
t39 = -t49 * t67 + t61 * t72;
t21 = pkin(8) * t39 + t26;
t66 = sin(qJ(4));
t71 = cos(qJ(4));
t13 = t66 * t16 + t71 * t21;
t28 = t39 * t71 - t40 * t66;
t11 = pkin(9) * t28 + t13;
t65 = sin(qJ(5));
t70 = cos(qJ(5));
t12 = t71 * t16 - t21 * t66;
t29 = t39 * t66 + t40 * t71;
t45 = qJD(4) + t47;
t8 = pkin(4) * t45 - pkin(9) * t29 + t12;
t6 = t70 * t11 + t65 * t8;
t36 = t41 * t63 - t62 * t46;
t5 = -t11 * t65 + t70 * t8;
t23 = t28 * t70 - t29 * t65;
t31 = -pkin(2) * t61 - t36;
t27 = -pkin(3) * t39 + t31;
t19 = -pkin(4) * t28 + t27;
t74 = V_base(3) ^ 2;
t69 = cos(qJ(6));
t64 = sin(qJ(6));
t44 = qJD(5) + t45;
t24 = t28 * t65 + t29 * t70;
t22 = qJD(6) - t23;
t18 = t24 * t69 + t44 * t64;
t17 = -t24 * t64 + t44 * t69;
t9 = -pkin(5) * t23 - pkin(10) * t24 + t19;
t4 = pkin(10) * t44 + t6;
t3 = -pkin(5) * t44 - t5;
t2 = t4 * t69 + t64 * t9;
t1 = -t4 * t64 + t69 * t9;
t7 = (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,2) * t54 / 0.2e1) * t54 + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,1) * t55 / 0.2e1) * t55 + (t31 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(4,5) * t47 + Ifges(4,1) * t40 / 0.2e1) * t40 + m(2) * (t50 ^ 2 + t51 ^ 2 + t74) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t74) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t52 ^ 2) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t19 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t44 + Ifges(6,1) * t24 / 0.2e1) * t24 + (t50 * mrSges(2,1) + t36 * mrSges(3,1) - t51 * mrSges(2,2) - t37 * mrSges(3,2) + Ifges(2,5) * t55 + Ifges(3,5) * t49 + Ifges(2,6) * t54 + Ifges(3,6) * t48 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t61) * t61 + (t52 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,1) * t49 / 0.2e1) * t49 + m(4) * (t25 ^ 2 + t26 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(6) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t31 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,6) * t47 + Ifges(4,2) * t39 / 0.2e1) * t39 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t22 / 0.2e1) * t22 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,6) * t45 + Ifges(5,2) * t28 / 0.2e1) * t28 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t44 / 0.2e1) * t44 + (-t19 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t24 + Ifges(6,6) * t44 + Ifges(6,2) * t23 / 0.2e1) * t23 + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t45 + Ifges(5,1) * t29 / 0.2e1) * t29 + (t25 * mrSges(4,1) - t26 * mrSges(4,2) + Ifges(4,3) * t47 / 0.2e1) * t47 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t18 + Ifges(7,6) * t22 + Ifges(7,2) * t17 / 0.2e1) * t17 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t22 + Ifges(7,1) * t18 / 0.2e1) * t18 + (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t45 / 0.2e1) * t45 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t52 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t49 + Ifges(3,2) * t48 / 0.2e1) * t48 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T  = t7;
