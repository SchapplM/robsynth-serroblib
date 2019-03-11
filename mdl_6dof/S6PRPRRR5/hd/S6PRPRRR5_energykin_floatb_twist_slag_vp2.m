% Calculate kinetic energy for
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:37
% EndTime: 2019-03-08 20:40:38
% DurationCPUTime: 1.15s
% Computational Cost: add. (3465->160), mult. (5927->226), div. (0->0), fcn. (4930->12), ass. (0->63)
t78 = pkin(2) + pkin(8);
t60 = sin(pkin(11));
t62 = cos(pkin(11));
t52 = t60 * V_base(5) + t62 * V_base(4);
t77 = pkin(7) * t52;
t67 = sin(qJ(2));
t51 = -t60 * V_base(4) + t62 * V_base(5);
t61 = sin(pkin(6));
t63 = cos(pkin(6));
t72 = t51 * t63 + t61 * V_base(6);
t76 = cos(qJ(2));
t41 = t52 * t76 + t67 * t72;
t47 = -t51 * t61 + t63 * V_base(6) + qJD(2);
t57 = V_base(5) * qJ(1) + V_base(1);
t58 = -V_base(4) * qJ(1) + V_base(2);
t49 = t62 * t57 + t60 * t58;
t38 = pkin(7) * t72 + t49;
t48 = -t57 * t60 + t62 * t58;
t42 = V_base(6) * pkin(1) - t63 * t77 + t48;
t59 = V_base(3) + qJD(1);
t45 = -pkin(1) * t51 - t61 * t77 + t59;
t74 = t63 * t76;
t75 = t61 * t76;
t29 = -t67 * t38 + t42 * t74 + t45 * t75;
t71 = qJD(3) - t29;
t17 = t41 * pkin(3) - t47 * t78 + t71;
t40 = -t51 * t74 + t52 * t67 - t75 * V_base(6);
t31 = -t42 * t61 + t63 * t45;
t73 = -qJ(3) * t41 + t31;
t19 = t40 * t78 + t73;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t13 = t66 * t17 + t70 * t19;
t32 = t40 * t70 - t66 * t47;
t11 = pkin(9) * t32 + t13;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t12 = t70 * t17 - t19 * t66;
t33 = t66 * t40 + t47 * t70;
t39 = qJD(4) + t41;
t9 = pkin(4) * t39 - pkin(9) * t33 + t12;
t6 = t69 * t11 + t65 * t9;
t30 = t76 * t38 + (t42 * t63 + t45 * t61) * t67;
t25 = -t47 * qJ(3) - t30;
t5 = -t11 * t65 + t69 * t9;
t27 = t32 * t69 - t33 * t65;
t22 = -pkin(3) * t40 - t25;
t14 = -pkin(4) * t32 + t22;
t68 = cos(qJ(6));
t64 = sin(qJ(6));
t37 = qJD(5) + t39;
t28 = t32 * t65 + t33 * t69;
t26 = qJD(6) - t27;
t24 = -t47 * pkin(2) + t71;
t23 = pkin(2) * t40 + t73;
t21 = t28 * t68 + t37 * t64;
t20 = -t28 * t64 + t37 * t68;
t7 = -pkin(5) * t27 - pkin(10) * t28 + t14;
t4 = pkin(10) * t37 + t6;
t3 = -pkin(5) * t37 - t5;
t2 = t4 * t68 + t64 * t7;
t1 = -t4 * t64 + t68 * t7;
t8 = m(2) * (t48 ^ 2 + t49 ^ 2 + t59 ^ 2) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t59 * mrSges(2,2) - t48 * mrSges(2,3) + Ifges(2,1) * t52 / 0.2e1) * t52 + (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t39 / 0.2e1) * t39 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t37 / 0.2e1) * t37 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t59 * mrSges(2,1) + t49 * mrSges(2,3) + Ifges(2,4) * t52 + Ifges(2,2) * t51 / 0.2e1) * t51 + (t22 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t39 + Ifges(5,1) * t33 / 0.2e1) * t33 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t37 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t26 + Ifges(7,1) * t21 / 0.2e1) * t21 + (-t22 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,6) * t39 + Ifges(5,2) * t32 / 0.2e1) * t32 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t37 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t26 + Ifges(7,2) * t20 / 0.2e1) * t20 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + t24 * mrSges(4,2) - t25 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t47) * t47 + (t24 * mrSges(4,1) + t31 * mrSges(3,2) - t29 * mrSges(3,3) - t23 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t41 + (-Ifges(4,4) + Ifges(3,5)) * t47) * t41 + (V_base(2) * mrSges(1,1) + t48 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t49 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t52 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t51 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t31 * mrSges(3,1) + t25 * mrSges(4,1) - t23 * mrSges(4,2) - t30 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t40 + (Ifges(4,5) - Ifges(3,6)) * t47 + (-Ifges(3,4) - Ifges(4,6)) * t41) * t40;
T  = t8;
