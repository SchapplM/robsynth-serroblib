% Calculate kinetic energy for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:40
% EndTime: 2019-03-09 08:49:41
% DurationCPUTime: 1.20s
% Computational Cost: add. (4157->156), mult. (5519->226), div. (0->0), fcn. (4524->12), ass. (0->59)
t69 = sin(qJ(1));
t73 = cos(qJ(1));
t54 = -t69 * V_base(4) + t73 * V_base(5);
t55 = t69 * V_base(5) + t73 * V_base(4);
t44 = -pkin(1) * t54 - pkin(7) * t55 + V_base(3);
t59 = V_base(5) * pkin(6) + V_base(1);
t60 = -V_base(4) * pkin(6) + V_base(2);
t52 = t73 * t59 + t69 * t60;
t61 = V_base(6) + qJD(1);
t48 = pkin(7) * t61 + t52;
t68 = sin(qJ(2));
t72 = cos(qJ(2));
t35 = t72 * t44 - t48 * t68;
t50 = t55 * t72 + t61 * t68;
t53 = qJD(2) - t54;
t29 = pkin(2) * t53 - qJ(3) * t50 + t35;
t36 = t68 * t44 + t72 * t48;
t49 = -t55 * t68 + t61 * t72;
t32 = qJ(3) * t49 + t36;
t63 = sin(pkin(10));
t65 = cos(pkin(10));
t22 = t63 * t29 + t65 * t32;
t20 = qJ(4) * t53 + t22;
t39 = -t65 * t49 + t50 * t63;
t40 = t49 * t63 + t50 * t65;
t51 = -t69 * t59 + t60 * t73;
t47 = -pkin(1) * t61 - t51;
t41 = -pkin(2) * t49 + qJD(3) + t47;
t25 = pkin(3) * t39 - qJ(4) * t40 + t41;
t62 = sin(pkin(11));
t64 = cos(pkin(11));
t14 = t64 * t20 + t62 * t25;
t33 = -t40 * t62 + t53 * t64;
t11 = pkin(8) * t33 + t14;
t67 = sin(qJ(5));
t71 = cos(qJ(5));
t13 = -t20 * t62 + t64 * t25;
t34 = t40 * t64 + t53 * t62;
t9 = pkin(4) * t39 - pkin(8) * t34 + t13;
t6 = t71 * t11 + t67 * t9;
t5 = -t11 * t67 + t71 * t9;
t21 = t29 * t65 - t63 * t32;
t38 = qJD(5) + t39;
t19 = -pkin(3) * t53 + qJD(4) - t21;
t17 = -pkin(4) * t33 + t19;
t74 = V_base(3) ^ 2;
t70 = cos(qJ(6));
t66 = sin(qJ(6));
t37 = qJD(6) + t38;
t27 = t33 * t67 + t34 * t71;
t26 = t33 * t71 - t34 * t67;
t16 = t26 * t66 + t27 * t70;
t15 = t26 * t70 - t27 * t66;
t12 = -pkin(5) * t26 + t17;
t4 = pkin(9) * t26 + t6;
t3 = pkin(5) * t38 - pkin(9) * t27 + t5;
t2 = t3 * t66 + t4 * t70;
t1 = t3 * t70 - t4 * t66;
t7 = (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t37 / 0.2e1) * t37 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t19 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,1) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t19 * mrSges(5,1) + t14 * mrSges(5,3) + Ifges(5,4) * t34 + Ifges(5,2) * t33 / 0.2e1) * t33 + (t12 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t37 + Ifges(7,1) * t16 / 0.2e1) * t16 + (-V_base(3) * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t61 + Ifges(2,2) * t54 / 0.2e1) * t54 + m(2) * (t51 ^ 2 + t52 ^ 2 + t74) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t74) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t41 ^ 2) / 0.2e1 + m(3) * (t35 ^ 2 + t36 ^ 2 + t47 ^ 2) / 0.2e1 + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) / 0.2e1 + (t41 * mrSges(4,1) + t13 * mrSges(5,1) - t14 * mrSges(5,2) - t22 * mrSges(4,3) - Ifges(4,4) * t40 + Ifges(5,5) * t34 - Ifges(4,6) * t53 + Ifges(5,6) * t33 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t39) * t39 + (-t47 * mrSges(3,1) + t36 * mrSges(3,3) + Ifges(3,4) * t50 + Ifges(3,2) * t49 / 0.2e1) * t49 + (t47 * mrSges(3,2) - t35 * mrSges(3,3) + Ifges(3,1) * t50 / 0.2e1) * t50 + (t41 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,1) * t40 / 0.2e1) * t40 + (t17 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t38 + Ifges(6,1) * t27 / 0.2e1) * t27 + (-t17 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t38 + Ifges(6,2) * t26 / 0.2e1) * t26 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t55 / 0.2e1) * t55 + (-t12 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t37 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t51 * mrSges(2,1) - t52 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (t35 * mrSges(3,1) + t21 * mrSges(4,1) - t36 * mrSges(3,2) - t22 * mrSges(4,2) + Ifges(3,5) * t50 + Ifges(4,5) * t40 + Ifges(3,6) * t49 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t53) * t53;
T  = t7;
