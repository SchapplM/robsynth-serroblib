% Calculate kinetic energy for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:03
% EndTime: 2019-03-08 21:08:04
% DurationCPUTime: 0.93s
% Computational Cost: add. (2509->158), mult. (4275->209), div. (0->0), fcn. (3494->10), ass. (0->57)
t76 = -pkin(4) - pkin(9);
t59 = sin(pkin(10));
t61 = cos(pkin(10));
t51 = t59 * V_base(5) + t61 * V_base(4);
t75 = pkin(7) * t51;
t74 = cos(qJ(2));
t56 = V_base(5) * qJ(1) + V_base(1);
t57 = -V_base(4) * qJ(1) + V_base(2);
t47 = -t59 * t56 + t61 * t57;
t73 = cos(pkin(6));
t39 = V_base(6) * pkin(1) - t73 * t75 + t47;
t50 = -t59 * V_base(4) + t61 * V_base(5);
t58 = V_base(3) + qJD(1);
t60 = sin(pkin(6));
t43 = -pkin(1) * t50 - t60 * t75 + t58;
t24 = -t39 * t60 + t73 * t43;
t64 = sin(qJ(2));
t70 = t73 * t74;
t72 = t60 * t74;
t37 = t50 * t70 - t51 * t64 + V_base(6) * t72;
t67 = t73 * t50 + t60 * V_base(6);
t38 = t74 * t51 + t67 * t64;
t15 = -pkin(2) * t37 - pkin(8) * t38 + t24;
t48 = t61 * t56 + t59 * t57;
t35 = t67 * pkin(7) + t48;
t21 = t74 * t35 + (t39 * t73 + t43 * t60) * t64;
t46 = -t60 * t50 + t73 * V_base(6) + qJD(2);
t19 = pkin(8) * t46 + t21;
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t12 = t63 * t15 + t66 * t19;
t36 = qJD(3) - t37;
t10 = t36 * qJ(4) + t12;
t20 = -t64 * t35 + t39 * t70 + t43 * t72;
t11 = t15 * t66 - t63 * t19;
t18 = -t46 * pkin(2) - t20;
t71 = qJD(4) - t11;
t29 = t38 * t66 + t63 * t46;
t28 = t38 * t63 - t66 * t46;
t7 = -qJ(5) * t28 - t10;
t13 = t28 * pkin(3) - t29 * qJ(4) + t18;
t69 = qJD(5) - t13;
t68 = -t29 * qJ(5) + t71;
t65 = cos(qJ(6));
t62 = sin(qJ(6));
t27 = qJD(6) + t29;
t23 = t28 * t65 - t36 * t62;
t22 = -t28 * t62 - t36 * t65;
t9 = -t36 * pkin(3) + t71;
t8 = -pkin(4) * t28 + t69;
t6 = pkin(5) * t36 - t7;
t5 = (-pkin(3) - pkin(4)) * t36 + t68;
t4 = pkin(5) * t29 + t76 * t28 + t69;
t3 = (-pkin(3) + t76) * t36 + t68;
t2 = t3 * t65 + t4 * t62;
t1 = -t3 * t62 + t4 * t65;
t14 = m(2) * (t47 ^ 2 + t48 ^ 2 + t58 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t24 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t18 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t13 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t58 * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,1) * t51 / 0.2e1) * t51 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t58 * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,2) * t50 / 0.2e1) * t50 + (t24 * mrSges(3,2) - t20 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t6 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t23 / 0.2e1) * t23 + (-t24 * mrSges(3,1) + t21 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t46 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t6 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t23 + Ifges(7,6) * t27 + Ifges(7,2) * t22 / 0.2e1) * t22 + (V_base(2) * mrSges(1,1) + t47 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t48 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t51 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t50 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t11 * mrSges(4,1) - t9 * mrSges(5,1) - t7 * mrSges(6,1) - t12 * mrSges(4,2) + t5 * mrSges(6,2) + t10 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t36) * t36 + (t8 * mrSges(6,1) + t18 * mrSges(4,2) + t9 * mrSges(5,2) - t11 * mrSges(4,3) - t13 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t29 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t36) * t29 + (t18 * mrSges(4,1) + t13 * mrSges(5,1) - t10 * mrSges(5,2) + t8 * mrSges(6,2) - t12 * mrSges(4,3) - t7 * mrSges(6,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t28 + (-Ifges(6,5) - Ifges(4,6) + Ifges(5,6)) * t36 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t29) * t28;
T  = t14;
