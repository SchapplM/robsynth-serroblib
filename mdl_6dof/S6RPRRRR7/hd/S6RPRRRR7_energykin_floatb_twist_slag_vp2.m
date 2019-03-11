% Calculate kinetic energy for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:52
% EndTime: 2019-03-09 07:16:53
% DurationCPUTime: 1.10s
% Computational Cost: add. (2587->153), mult. (3277->213), div. (0->0), fcn. (2464->10), ass. (0->57)
t69 = pkin(1) + pkin(7);
t60 = sin(qJ(1));
t65 = cos(qJ(1));
t47 = t60 * V_base(5) + t65 * V_base(4);
t55 = V_base(6) + qJD(1);
t51 = V_base(5) * pkin(6) + V_base(1);
t52 = -V_base(4) * pkin(6) + V_base(2);
t41 = -t60 * t51 + t52 * t65;
t67 = qJD(2) - t41;
t32 = pkin(2) * t47 - t55 * t69 + t67;
t46 = t60 * V_base(4) - t65 * V_base(5);
t68 = -qJ(2) * t47 + V_base(3);
t34 = t46 * t69 + t68;
t59 = sin(qJ(3));
t64 = cos(qJ(3));
t25 = t64 * t32 - t34 * t59;
t40 = t46 * t59 + t55 * t64;
t45 = qJD(3) + t47;
t18 = pkin(3) * t45 - pkin(8) * t40 + t25;
t26 = t59 * t32 + t64 * t34;
t39 = t46 * t64 - t55 * t59;
t20 = pkin(8) * t39 + t26;
t58 = sin(qJ(4));
t63 = cos(qJ(4));
t13 = t58 * t18 + t63 * t20;
t28 = t39 * t63 - t40 * t58;
t11 = pkin(9) * t28 + t13;
t57 = sin(qJ(5));
t62 = cos(qJ(5));
t12 = t63 * t18 - t20 * t58;
t29 = t39 * t58 + t40 * t63;
t44 = qJD(4) + t45;
t9 = pkin(4) * t44 - pkin(9) * t29 + t12;
t6 = t62 * t11 + t57 * t9;
t42 = t65 * t51 + t60 * t52;
t38 = -t55 * qJ(2) - t42;
t5 = -t11 * t57 + t62 * t9;
t23 = t28 * t62 - t29 * t57;
t35 = -pkin(2) * t46 - t38;
t27 = -pkin(3) * t39 + t35;
t21 = -pkin(4) * t28 + t27;
t66 = V_base(3) ^ 2;
t61 = cos(qJ(6));
t56 = sin(qJ(6));
t43 = qJD(5) + t44;
t37 = -pkin(1) * t55 + t67;
t36 = pkin(1) * t46 + t68;
t24 = t28 * t57 + t29 * t62;
t22 = qJD(6) - t23;
t17 = t24 * t61 + t43 * t56;
t16 = -t24 * t56 + t43 * t61;
t7 = -pkin(5) * t23 - pkin(10) * t24 + t21;
t4 = pkin(10) * t43 + t6;
t3 = -pkin(5) * t43 - t5;
t2 = t4 * t61 + t56 * t7;
t1 = -t4 * t56 + t61 * t7;
t8 = (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t44 / 0.2e1) * t44 + (-t21 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t24 + Ifges(6,6) * t43 + Ifges(6,2) * t23 / 0.2e1) * t23 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t22 / 0.2e1) * t22 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t22 + Ifges(7,2) * t16 / 0.2e1) * t16 + (-t35 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,6) * t45 + Ifges(4,2) * t39 / 0.2e1) * t39 + (t35 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(4,5) * t45 + Ifges(4,1) * t40 / 0.2e1) * t40 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t22 + Ifges(7,1) * t17 / 0.2e1) * t17 + m(2) * (t41 ^ 2 + t42 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(4) * (t25 ^ 2 + t26 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(6) * (t21 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(3) * mrSges(2,1) + t38 * mrSges(3,1) - t36 * mrSges(3,2) - t42 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t46 + (Ifges(3,5) - Ifges(2,6)) * t55 + (-Ifges(2,4) - Ifges(3,6)) * t47) * t46 + (t37 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t41 * mrSges(2,3) - t36 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t47 + (-Ifges(3,4) + Ifges(2,5)) * t55) * t47 + (t25 * mrSges(4,1) - t26 * mrSges(4,2) + Ifges(4,3) * t45 / 0.2e1) * t45 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t21 * mrSges(6,2) - t5 * mrSges(6,3) + t43 * Ifges(6,5) + Ifges(6,1) * t24 / 0.2e1) * t24 + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t44 + Ifges(5,1) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t43 / 0.2e1) * t43 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,6) * t44 + Ifges(5,2) * t28 / 0.2e1) * t28 + (t41 * mrSges(2,1) - t42 * mrSges(2,2) + t37 * mrSges(3,2) - t38 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t55) * t55;
T  = t8;
