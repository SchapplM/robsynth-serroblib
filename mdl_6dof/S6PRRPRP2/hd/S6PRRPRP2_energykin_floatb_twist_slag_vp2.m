% Calculate kinetic energy for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:07
% EndTime: 2019-03-08 21:28:08
% DurationCPUTime: 1.16s
% Computational Cost: add. (4419->157), mult. (7627->220), div. (0->0), fcn. (6454->12), ass. (0->58)
t57 = V_base(5) * qJ(1) + V_base(1);
t58 = -V_base(4) * qJ(1) + V_base(2);
t61 = sin(pkin(10));
t64 = cos(pkin(10));
t50 = -t57 * t61 + t64 * t58;
t65 = cos(pkin(6));
t53 = t61 * V_base(5) + t64 * V_base(4);
t77 = pkin(7) * t53;
t45 = V_base(6) * pkin(1) - t65 * t77 + t50;
t52 = -t61 * V_base(4) + t64 * V_base(5);
t59 = V_base(3) + qJD(1);
t62 = sin(pkin(6));
t48 = -pkin(1) * t52 - t62 * t77 + t59;
t78 = t45 * t65 + t48 * t62;
t51 = t64 * t57 + t61 * t58;
t71 = t52 * t65 + t62 * V_base(6);
t41 = t71 * pkin(7) + t51;
t68 = sin(qJ(2));
t70 = cos(qJ(2));
t32 = -t68 * t41 + t78 * t70;
t43 = -t68 * t53 + t71 * t70;
t49 = -t52 * t62 + t65 * V_base(6) + qJD(2);
t27 = -t49 * pkin(2) - t32;
t44 = t53 * t70 + t71 * t68;
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t35 = -t44 * t67 + t49 * t69;
t20 = -t35 * pkin(3) + qJD(4) + t27;
t36 = t44 * t69 + t49 * t67;
t60 = sin(pkin(11));
t63 = cos(pkin(11));
t30 = t35 * t63 - t36 * t60;
t31 = t35 * t60 + t36 * t63;
t12 = -t30 * pkin(4) - t31 * pkin(9) + t20;
t66 = sin(qJ(5));
t76 = cos(qJ(5));
t34 = -t45 * t62 + t65 * t48;
t25 = -pkin(2) * t43 - pkin(8) * t44 + t34;
t33 = t70 * t41 + t78 * t68;
t28 = pkin(8) * t49 + t33;
t18 = t69 * t25 - t28 * t67;
t42 = qJD(3) - t43;
t14 = pkin(3) * t42 - qJ(4) * t36 + t18;
t19 = t67 * t25 + t69 * t28;
t17 = qJ(4) * t35 + t19;
t10 = t60 * t14 + t63 * t17;
t8 = pkin(9) * t42 + t10;
t4 = t66 * t12 + t76 * t8;
t9 = t14 * t63 - t60 * t17;
t7 = -pkin(4) * t42 - t9;
t3 = t76 * t12 - t66 * t8;
t29 = qJD(5) - t30;
t22 = t76 * t31 + t66 * t42;
t21 = t31 * t66 - t76 * t42;
t5 = pkin(5) * t21 - qJ(6) * t22 + t7;
t2 = qJ(6) * t29 + t4;
t1 = -t29 * pkin(5) + qJD(6) - t3;
t6 = m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t27 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t20 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t50 ^ 2 + t51 ^ 2 + t59 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t59 * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,1) * t53 / 0.2e1) * t53 + (t32 * mrSges(3,1) - t33 * mrSges(3,2) + Ifges(3,3) * t49 / 0.2e1) * t49 + (t27 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t36 / 0.2e1) * t36 + (t20 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,1) * t31 / 0.2e1) * t31 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t59 * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t53 + Ifges(2,2) * t52 / 0.2e1) * t52 + (t34 * mrSges(3,2) - t32 * mrSges(3,3) + Ifges(3,5) * t49 + Ifges(3,1) * t44 / 0.2e1) * t44 + (-t27 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t36 + Ifges(4,2) * t35 / 0.2e1) * t35 + (-t20 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,2) * t30 / 0.2e1) * t30 + (-t34 * mrSges(3,1) + t33 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,6) * t49 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t29) * t29 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t29) * t22 + (V_base(2) * mrSges(1,1) + t50 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t51 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t53 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t52 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t18 * mrSges(4,1) + t9 * mrSges(5,1) - t19 * mrSges(4,2) - t10 * mrSges(5,2) + Ifges(4,5) * t36 + Ifges(5,5) * t31 + Ifges(4,6) * t35 + Ifges(5,6) * t30 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t42) * t42 + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t29 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21;
T  = t6;
