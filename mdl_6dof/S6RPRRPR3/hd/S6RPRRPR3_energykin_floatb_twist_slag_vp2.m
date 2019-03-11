% Calculate kinetic energy for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:20
% EndTime: 2019-03-09 05:04:21
% DurationCPUTime: 1.19s
% Computational Cost: add. (2555->153), mult. (3675->211), div. (0->0), fcn. (2916->10), ass. (0->56)
t72 = -pkin(4) - pkin(5);
t71 = cos(qJ(4));
t55 = pkin(6) * V_base(5) + V_base(1);
t56 = -pkin(6) * V_base(4) + V_base(2);
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t46 = -t55 * t63 + t56 * t66;
t50 = t63 * V_base(5) + t66 * V_base(4);
t57 = V_base(6) + qJD(1);
t38 = pkin(1) * t57 - qJ(2) * t50 + t46;
t47 = t55 * t66 + t56 * t63;
t49 = -t63 * V_base(4) + t66 * V_base(5);
t42 = qJ(2) * t49 + t47;
t58 = sin(pkin(10));
t59 = cos(pkin(10));
t31 = t38 * t58 + t42 * t59;
t24 = pkin(7) * t57 + t31;
t44 = t49 * t59 - t50 * t58;
t45 = t49 * t58 + t50 * t59;
t48 = -pkin(1) * t49 + qJD(2) + V_base(3);
t27 = -pkin(2) * t44 - pkin(7) * t45 + t48;
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t17 = t24 * t65 + t27 * t62;
t43 = qJD(3) - t44;
t14 = pkin(8) * t43 + t17;
t30 = t38 * t59 - t42 * t58;
t23 = -pkin(2) * t57 - t30;
t36 = -t45 * t62 + t57 * t65;
t37 = t45 * t65 + t57 * t62;
t18 = -pkin(3) * t36 - pkin(8) * t37 + t23;
t61 = sin(qJ(4));
t9 = t14 * t71 + t18 * t61;
t16 = -t24 * t62 + t27 * t65;
t34 = qJD(4) - t36;
t7 = qJ(5) * t34 + t9;
t70 = pkin(3) * t43 + t16;
t8 = -t14 * t61 + t18 * t71;
t69 = qJD(5) - t8;
t29 = t37 * t71 + t43 * t61;
t68 = qJ(5) * t29 + t70;
t67 = V_base(3) ^ 2;
t64 = cos(qJ(6));
t60 = sin(qJ(6));
t33 = qJD(6) - t34;
t28 = t37 * t61 - t43 * t71;
t20 = t28 * t60 + t29 * t64;
t19 = t28 * t64 - t29 * t60;
t10 = pkin(4) * t28 - t68;
t6 = -pkin(4) * t34 + t69;
t5 = t28 * t72 + t68;
t4 = pkin(9) * t28 + t7;
t3 = -pkin(9) * t29 + t34 * t72 + t69;
t2 = t3 * t60 + t4 * t64;
t1 = t3 * t64 - t4 * t60;
t11 = m(2) * (t46 ^ 2 + t47 ^ 2 + t67) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t67) / 0.2e1 + m(3) * (t30 ^ 2 + t31 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t23 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t46 * mrSges(2,3) + Ifges(2,1) * t50 / 0.2e1) * t50 + (t48 * mrSges(3,2) - t30 * mrSges(3,3) + Ifges(3,1) * t45 / 0.2e1) * t45 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t43 / 0.2e1) * t43 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t33 / 0.2e1) * t33 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t47 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,2) * t49 / 0.2e1) * t49 + (-t48 * mrSges(3,1) + t31 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,2) * t44 / 0.2e1) * t44 + (t23 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t43 + Ifges(4,1) * t37 / 0.2e1) * t37 + (t5 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t33 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t23 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,6) * t43 + Ifges(4,2) * t36 / 0.2e1) * t36 + (-t5 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t33 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t8 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) + t7 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t34) * t34 + (-t70 * mrSges(5,2) + t6 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t29 + (Ifges(6,4) + Ifges(5,5)) * t34) * t29 + (t46 * mrSges(2,1) + t30 * mrSges(3,1) - t47 * mrSges(2,2) - t31 * mrSges(3,2) + Ifges(2,5) * t50 + Ifges(3,5) * t45 + Ifges(2,6) * t49 + Ifges(3,6) * t44 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t57) * t57 + (-t70 * mrSges(5,1) + t10 * mrSges(6,1) - t7 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t28 + (-Ifges(5,6) + Ifges(6,6)) * t34 + (-Ifges(5,4) + Ifges(6,5)) * t29) * t28;
T  = t11;
