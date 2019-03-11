% Calculate kinetic energy for
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:22
% EndTime: 2019-03-10 01:34:24
% DurationCPUTime: 1.36s
% Computational Cost: add. (5127->157), mult. (7649->223), div. (0->0), fcn. (6476->12), ass. (0->60)
t59 = V_base(5) * pkin(7) + V_base(1);
t60 = -V_base(4) * pkin(7) + V_base(2);
t68 = sin(qJ(1));
t73 = cos(qJ(1));
t51 = -t59 * t68 + t73 * t60;
t61 = V_base(6) + qJD(1);
t63 = cos(pkin(6));
t55 = t68 * V_base(5) + t73 * V_base(4);
t79 = pkin(8) * t55;
t46 = pkin(1) * t61 - t63 * t79 + t51;
t54 = -t68 * V_base(4) + t73 * V_base(5);
t62 = sin(pkin(6));
t49 = -pkin(1) * t54 - t62 * t79 + V_base(3);
t80 = t46 * t63 + t49 * t62;
t52 = t73 * t59 + t68 * t60;
t75 = t54 * t63 + t61 * t62;
t43 = t75 * pkin(8) + t52;
t67 = sin(qJ(2));
t72 = cos(qJ(2));
t33 = -t67 * t43 + t80 * t72;
t44 = -t67 * t55 + t75 * t72;
t50 = -t54 * t62 + t61 * t63 + qJD(2);
t28 = -pkin(2) * t50 - t33;
t45 = t55 * t72 + t75 * t67;
t66 = sin(qJ(3));
t71 = cos(qJ(3));
t36 = -t45 * t66 + t50 * t71;
t21 = -pkin(3) * t36 + t28;
t37 = t45 * t71 + t50 * t66;
t65 = sin(qJ(4));
t70 = cos(qJ(4));
t31 = t36 * t70 - t37 * t65;
t32 = t36 * t65 + t37 * t70;
t13 = -pkin(4) * t31 - pkin(11) * t32 + t21;
t64 = sin(qJ(5));
t69 = cos(qJ(5));
t35 = -t46 * t62 + t63 * t49;
t26 = -pkin(2) * t44 - pkin(9) * t45 + t35;
t34 = t72 * t43 + t80 * t67;
t29 = pkin(9) * t50 + t34;
t19 = t71 * t26 - t29 * t66;
t42 = qJD(3) - t44;
t15 = pkin(3) * t42 - pkin(10) * t37 + t19;
t20 = t66 * t26 + t71 * t29;
t18 = pkin(10) * t36 + t20;
t10 = t65 * t15 + t70 * t18;
t41 = qJD(4) + t42;
t8 = pkin(11) * t41 + t10;
t4 = t64 * t13 + t69 * t8;
t3 = t69 * t13 - t64 * t8;
t9 = t15 * t70 - t65 * t18;
t7 = -pkin(4) * t41 - t9;
t74 = V_base(3) ^ 2;
t30 = qJD(5) - t31;
t23 = t32 * t69 + t41 * t64;
t22 = -t32 * t64 + t41 * t69;
t5 = -pkin(5) * t22 + qJD(6) + t7;
t2 = qJ(6) * t22 + t4;
t1 = pkin(5) * t30 - qJ(6) * t23 + t3;
t6 = m(2) * (t51 ^ 2 + t52 ^ 2 + t74) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t74) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t20 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t21 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (t33 * mrSges(3,1) - t34 * mrSges(3,2) + Ifges(3,3) * t50 / 0.2e1) * t50 + (-t35 * mrSges(3,1) + t34 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,6) * t50 + Ifges(3,2) * t44 / 0.2e1) * t44 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (-t28 * mrSges(4,1) + t20 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,6) * t42 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t51 * mrSges(2,1) - t52 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (-t21 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t32 + Ifges(5,6) * t41 + Ifges(5,2) * t31 / 0.2e1) * t31 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t28 * mrSges(4,2) - t19 * mrSges(4,3) + Ifges(4,5) * t42 + Ifges(4,1) * t37 / 0.2e1) * t37 + (t35 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,5) * t50 + Ifges(3,1) * t45 / 0.2e1) * t45 + (t19 * mrSges(4,1) - t20 * mrSges(4,2) + Ifges(4,3) * t42 / 0.2e1) * t42 + (-V_base(3) * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t61 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t21 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t32 / 0.2e1) * t32 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t30) * t30 + (V_base(3) * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t55 / 0.2e1) * t55 + (-t7 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t22 + (Ifges(6,6) + Ifges(7,6)) * t30 + (Ifges(6,4) + Ifges(7,4)) * t23) * t22 + (t7 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t23 + (Ifges(6,5) + Ifges(7,5)) * t30) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T  = t6;
