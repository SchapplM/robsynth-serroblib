% Calculate kinetic energy for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:17
% EndTime: 2019-03-09 04:57:18
% DurationCPUTime: 1.19s
% Computational Cost: add. (3953->156), mult. (5815->226), div. (0->0), fcn. (4792->12), ass. (0->59)
t58 = V_base(5) * pkin(6) + V_base(1);
t59 = -V_base(4) * pkin(6) + V_base(2);
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t49 = -t58 * t68 + t72 * t59;
t54 = t68 * V_base(5) + t72 * V_base(4);
t60 = V_base(6) + qJD(1);
t41 = pkin(1) * t60 - qJ(2) * t54 + t49;
t50 = t72 * t58 + t68 * t59;
t53 = -t68 * V_base(4) + t72 * V_base(5);
t45 = qJ(2) * t53 + t50;
t62 = sin(pkin(10));
t64 = cos(pkin(10));
t37 = t62 * t41 + t64 * t45;
t32 = pkin(7) * t60 + t37;
t47 = t53 * t64 - t62 * t54;
t48 = t53 * t62 + t54 * t64;
t51 = -pkin(1) * t53 + qJD(2) + V_base(3);
t35 = -pkin(2) * t47 - pkin(7) * t48 + t51;
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t25 = -t32 * t67 + t71 * t35;
t40 = t48 * t71 + t60 * t67;
t46 = qJD(3) - t47;
t16 = pkin(3) * t46 - pkin(8) * t40 + t25;
t26 = t71 * t32 + t67 * t35;
t39 = -t48 * t67 + t60 * t71;
t21 = pkin(8) * t39 + t26;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t13 = t66 * t16 + t70 * t21;
t28 = t39 * t70 - t40 * t66;
t11 = qJ(5) * t28 + t13;
t61 = sin(pkin(11));
t63 = cos(pkin(11));
t12 = t70 * t16 - t21 * t66;
t29 = t39 * t66 + t40 * t70;
t44 = qJD(4) + t46;
t9 = pkin(4) * t44 - qJ(5) * t29 + t12;
t6 = t63 * t11 + t61 * t9;
t36 = t41 * t64 - t62 * t45;
t5 = -t11 * t61 + t63 * t9;
t23 = t28 * t63 - t29 * t61;
t31 = -pkin(2) * t60 - t36;
t27 = -pkin(3) * t39 + t31;
t19 = -pkin(4) * t28 + qJD(5) + t27;
t73 = V_base(3) ^ 2;
t69 = cos(qJ(6));
t65 = sin(qJ(6));
t24 = t28 * t61 + t29 * t63;
t22 = qJD(6) - t23;
t18 = t24 * t69 + t44 * t65;
t17 = -t24 * t65 + t44 * t69;
t8 = -pkin(5) * t23 - pkin(9) * t24 + t19;
t4 = pkin(9) * t44 + t6;
t3 = -pkin(5) * t44 - t5;
t2 = t4 * t69 + t65 * t8;
t1 = -t4 * t65 + t69 * t8;
t7 = (t31 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t40 / 0.2e1) * t40 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t73) / 0.2e1 + m(2) * (t49 ^ 2 + t50 ^ 2 + t73) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t19 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t24 / 0.2e1) * t24 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,2) * t28 / 0.2e1) * t28 + (t49 * mrSges(2,1) + t36 * mrSges(3,1) - t50 * mrSges(2,2) - t37 * mrSges(3,2) + Ifges(2,5) * t54 + Ifges(3,5) * t48 + Ifges(2,6) * t53 + Ifges(3,6) * t47 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t60) * t60 + (t51 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,1) * t48 / 0.2e1) * t48 + (-t19 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t24 + Ifges(6,2) * t23 / 0.2e1) * t23 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t22 / 0.2e1) * t22 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t22 + Ifges(7,1) * t18 / 0.2e1) * t18 + m(3) * (t36 ^ 2 + t37 ^ 2 + t51 ^ 2) / 0.2e1 + m(4) * (t25 ^ 2 + t26 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(6) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t13 * mrSges(5,2) - t6 * mrSges(6,2) + Ifges(5,5) * t29 + Ifges(6,5) * t24 + Ifges(5,6) * t28 + Ifges(6,6) * t23 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t44) * t44 + (-t31 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,6) * t46 + Ifges(4,2) * t39 / 0.2e1) * t39 + (-V_base(3) * mrSges(2,1) + t50 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,2) * t53 / 0.2e1) * t53 + (-t51 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t48 + Ifges(3,2) * t47 / 0.2e1) * t47 + (t25 * mrSges(4,1) - t26 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t18 + Ifges(7,6) * t22 + Ifges(7,2) * t17 / 0.2e1) * t17 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(2,2) - t49 * mrSges(2,3) + Ifges(2,1) * t54 / 0.2e1) * t54;
T  = t7;
