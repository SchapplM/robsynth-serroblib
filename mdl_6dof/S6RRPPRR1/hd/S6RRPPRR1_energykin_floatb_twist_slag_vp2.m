% Calculate kinetic energy for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:42
% EndTime: 2019-03-09 08:45:43
% DurationCPUTime: 1.01s
% Computational Cost: add. (2579->153), mult. (3391->211), div. (0->0), fcn. (2660->10), ass. (0->54)
t65 = sin(qJ(1));
t69 = cos(qJ(1));
t50 = -t65 * V_base(4) + t69 * V_base(5);
t51 = t65 * V_base(5) + t69 * V_base(4);
t37 = -pkin(1) * t50 - pkin(7) * t51 + V_base(3);
t56 = V_base(5) * pkin(6) + V_base(1);
t57 = -V_base(4) * pkin(6) + V_base(2);
t46 = t69 * t56 + t65 * t57;
t60 = V_base(6) + qJD(1);
t41 = pkin(7) * t60 + t46;
t64 = sin(qJ(2));
t68 = cos(qJ(2));
t28 = t68 * t37 - t41 * t64;
t44 = t51 * t68 + t60 * t64;
t49 = qJD(2) - t50;
t24 = pkin(2) * t49 - qJ(3) * t44 + t28;
t29 = t64 * t37 + t68 * t41;
t43 = -t51 * t64 + t60 * t68;
t27 = qJ(3) * t43 + t29;
t61 = sin(pkin(10));
t72 = cos(pkin(10));
t16 = t61 * t24 + t72 * t27;
t14 = t49 * qJ(4) + t16;
t32 = -t72 * t43 + t44 * t61;
t11 = pkin(8) * t32 + t14;
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t33 = t61 * t43 + t72 * t44;
t15 = t72 * t24 - t61 * t27;
t71 = qJD(4) - t15;
t9 = -t33 * pkin(8) + (-pkin(3) - pkin(4)) * t49 + t71;
t6 = t67 * t11 + t63 * t9;
t45 = -t65 * t56 + t69 * t57;
t40 = -t60 * pkin(1) - t45;
t5 = -t11 * t63 + t67 * t9;
t34 = -t43 * pkin(2) + qJD(3) + t40;
t21 = t32 * t67 - t33 * t63;
t17 = t32 * pkin(3) - t33 * qJ(4) + t34;
t12 = -pkin(4) * t32 - t17;
t70 = V_base(3) ^ 2;
t66 = cos(qJ(6));
t62 = sin(qJ(6));
t48 = qJD(5) - t49;
t22 = t32 * t63 + t33 * t67;
t20 = qJD(6) - t21;
t19 = t22 * t66 + t48 * t62;
t18 = -t22 * t62 + t48 * t66;
t13 = -t49 * pkin(3) + t71;
t7 = -pkin(5) * t21 - pkin(9) * t22 + t12;
t4 = pkin(9) * t48 + t6;
t3 = -pkin(5) * t48 - t5;
t2 = t4 * t66 + t62 * t7;
t1 = -t4 * t62 + t66 * t7;
t8 = m(2) * (t45 ^ 2 + t46 ^ 2 + t70) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t70) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t34 ^ 2) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t40 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t17 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t45 * mrSges(2,1) - t46 * mrSges(2,2) + Ifges(2,3) * t60 / 0.2e1) * t60 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t48 / 0.2e1) * t48 + (t40 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,1) * t44 / 0.2e1) * t44 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,5) * t60 + Ifges(2,1) * t51 / 0.2e1) * t51 + (-t40 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t12 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t48 + Ifges(6,1) * t22 / 0.2e1) * t22 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t20 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t46 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,6) * t60 + Ifges(2,2) * t50 / 0.2e1) * t50 + (-t12 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t22 + Ifges(6,6) * t48 + Ifges(6,2) * t21 / 0.2e1) * t21 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t20 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t34 * mrSges(4,2) + t13 * mrSges(5,2) - t15 * mrSges(4,3) - t17 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t33) * t33 + (t34 * mrSges(4,1) + t17 * mrSges(5,1) - t14 * mrSges(5,2) - t16 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t32 + (-Ifges(4,4) + Ifges(5,5)) * t33) * t32 + (t28 * mrSges(3,1) + t15 * mrSges(4,1) - t13 * mrSges(5,1) - t29 * mrSges(3,2) - t16 * mrSges(4,2) + t14 * mrSges(5,3) + Ifges(3,5) * t44 + Ifges(3,6) * t43 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t49 + (Ifges(5,4) + Ifges(4,5)) * t33 + (-Ifges(4,6) + Ifges(5,6)) * t32) * t49;
T  = t8;
