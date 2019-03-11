% Calculate kinetic energy for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:03
% EndTime: 2019-03-09 00:01:04
% DurationCPUTime: 1.26s
% Computational Cost: add. (4455->157), mult. (7627->222), div. (0->0), fcn. (6454->12), ass. (0->59)
t59 = qJ(1) * V_base(5) + V_base(1);
t60 = -qJ(1) * V_base(4) + V_base(2);
t62 = sin(pkin(11));
t64 = cos(pkin(11));
t51 = -t59 * t62 + t60 * t64;
t65 = cos(pkin(6));
t55 = t62 * V_base(5) + t64 * V_base(4);
t78 = pkin(7) * t55;
t46 = V_base(6) * pkin(1) - t65 * t78 + t51;
t54 = -t62 * V_base(4) + t64 * V_base(5);
t61 = V_base(3) + qJD(1);
t63 = sin(pkin(6));
t49 = -pkin(1) * t54 - t63 * t78 + t61;
t79 = t46 * t65 + t49 * t63;
t52 = t59 * t64 + t60 * t62;
t73 = t54 * t65 + t63 * V_base(6);
t42 = pkin(7) * t73 + t52;
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t32 = -t42 * t69 + t72 * t79;
t44 = -t55 * t69 + t72 * t73;
t50 = -t54 * t63 + t65 * V_base(6) + qJD(2);
t27 = -t50 * pkin(2) - t32;
t45 = t55 * t72 + t69 * t73;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t35 = -t45 * t68 + t50 * t71;
t20 = -t35 * pkin(3) + t27;
t36 = t45 * t71 + t50 * t68;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t30 = t35 * t70 - t36 * t67;
t31 = t35 * t67 + t36 * t70;
t12 = -t30 * pkin(4) - t31 * pkin(10) + t20;
t66 = sin(qJ(5));
t77 = cos(qJ(5));
t34 = -t46 * t63 + t49 * t65;
t25 = -pkin(2) * t44 - pkin(8) * t45 + t34;
t33 = t72 * t42 + t69 * t79;
t28 = pkin(8) * t50 + t33;
t18 = t25 * t71 - t28 * t68;
t43 = qJD(3) - t44;
t14 = pkin(3) * t43 - pkin(9) * t36 + t18;
t19 = t25 * t68 + t28 * t71;
t17 = pkin(9) * t35 + t19;
t10 = t14 * t67 + t17 * t70;
t41 = qJD(4) + t43;
t8 = pkin(10) * t41 + t10;
t4 = t12 * t66 + t77 * t8;
t9 = t14 * t70 - t17 * t67;
t7 = -pkin(4) * t41 - t9;
t3 = t12 * t77 - t66 * t8;
t29 = qJD(5) - t30;
t22 = t31 * t77 + t41 * t66;
t21 = t31 * t66 - t41 * t77;
t5 = pkin(5) * t21 - qJ(6) * t22 + t7;
t2 = qJ(6) * t29 + t4;
t1 = -pkin(5) * t29 + qJD(6) - t3;
t6 = m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t27 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t20 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t51 ^ 2 + t52 ^ 2 + t61 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t61 * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,1) * t55 / 0.2e1) * t55 + (t32 * mrSges(3,1) - t33 * mrSges(3,2) + Ifges(3,3) * t50 / 0.2e1) * t50 + (t18 * mrSges(4,1) - t19 * mrSges(4,2) + Ifges(4,3) * t43 / 0.2e1) * t43 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t61 * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t34 * mrSges(3,2) - t32 * mrSges(3,3) + Ifges(3,5) * t50 + Ifges(3,1) * t45 / 0.2e1) * t45 + (t27 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,5) * t43 + Ifges(4,1) * t36 / 0.2e1) * t36 + (t20 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t31 / 0.2e1) * t31 + (-t34 * mrSges(3,1) + t33 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,6) * t50 + Ifges(3,2) * t44 / 0.2e1) * t44 + (-t27 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t36 + Ifges(4,6) * t43 + Ifges(4,2) * t35 / 0.2e1) * t35 + (-t20 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,6) * t41 + Ifges(5,2) * t30 / 0.2e1) * t30 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t29) * t29 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t29) * t22 + (V_base(2) * mrSges(1,1) + t51 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t52 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t55 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t54 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t29 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21;
T  = t6;
