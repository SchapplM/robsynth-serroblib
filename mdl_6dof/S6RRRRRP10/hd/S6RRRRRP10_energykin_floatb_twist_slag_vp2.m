% Calculate kinetic energy for
% S6RRRRRP10
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:01
% EndTime: 2019-03-10 02:16:03
% DurationCPUTime: 1.34s
% Computational Cost: add. (5107->157), mult. (7627->223), div. (0->0), fcn. (6450->12), ass. (0->60)
t59 = V_base(5) * pkin(7) + V_base(1);
t60 = -V_base(4) * pkin(7) + V_base(2);
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t52 = -t59 * t68 + t72 * t60;
t61 = V_base(6) + qJD(1);
t63 = cos(pkin(6));
t55 = t68 * V_base(5) + t72 * V_base(4);
t80 = pkin(8) * t55;
t47 = pkin(1) * t61 - t63 * t80 + t52;
t54 = -t68 * V_base(4) + t72 * V_base(5);
t62 = sin(pkin(6));
t50 = -pkin(1) * t54 - t62 * t80 + V_base(3);
t81 = t47 * t63 + t50 * t62;
t53 = t72 * t59 + t68 * t60;
t74 = t54 * t63 + t61 * t62;
t44 = t74 * pkin(8) + t53;
t67 = sin(qJ(2));
t71 = cos(qJ(2));
t29 = -t67 * t44 + t81 * t71;
t45 = -t55 * t67 + t74 * t71;
t64 = sin(qJ(5));
t34 = -t47 * t62 + t63 * t50;
t46 = t55 * t71 + t74 * t67;
t24 = -pkin(2) * t45 - pkin(9) * t46 + t34;
t30 = t71 * t44 + t81 * t67;
t51 = -t54 * t62 + t61 * t63 + qJD(2);
t28 = pkin(9) * t51 + t30;
t66 = sin(qJ(3));
t70 = cos(qJ(3));
t17 = t66 * t24 + t70 * t28;
t43 = qJD(3) - t45;
t15 = pkin(10) * t43 + t17;
t27 = -pkin(2) * t51 - t29;
t37 = -t66 * t46 + t51 * t70;
t38 = t46 * t70 + t51 * t66;
t20 = -pkin(3) * t37 - pkin(10) * t38 + t27;
t65 = sin(qJ(4));
t69 = cos(qJ(4));
t10 = -t15 * t65 + t69 * t20;
t33 = t38 * t69 + t43 * t65;
t36 = qJD(4) - t37;
t7 = pkin(4) * t36 - pkin(11) * t33 + t10;
t79 = cos(qJ(5));
t11 = t69 * t15 + t65 * t20;
t32 = -t38 * t65 + t43 * t69;
t9 = pkin(11) * t32 + t11;
t4 = t64 * t7 + t79 * t9;
t16 = t24 * t70 - t66 * t28;
t14 = -pkin(3) * t43 - t16;
t3 = -t64 * t9 + t79 * t7;
t12 = -pkin(4) * t32 + t14;
t73 = V_base(3) ^ 2;
t35 = qJD(5) + t36;
t22 = t64 * t32 + t79 * t33;
t21 = -t79 * t32 + t33 * t64;
t5 = pkin(5) * t21 - qJ(6) * t22 + t12;
t2 = qJ(6) * t35 + t4;
t1 = -t35 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t52 ^ 2 + t53 ^ 2 + t73) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t73) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t27 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (-t27 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t38 + Ifges(4,6) * t43 + Ifges(4,2) * t37 / 0.2e1) * t37 + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t35 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21 + (t52 * mrSges(2,1) - t53 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (-t14 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,6) * t36 + Ifges(5,2) * t32 / 0.2e1) * t32 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t35) * t22 + (V_base(3) * mrSges(2,2) - t52 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t55 / 0.2e1) * t55 + (t34 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t51 + Ifges(3,1) * t46 / 0.2e1) * t46 + (t27 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t43 + Ifges(4,1) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t14 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,5) * t36 + Ifges(5,1) * t33 / 0.2e1) * t33 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t43 / 0.2e1) * t43 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t53 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t61 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t51 / 0.2e1) * t51 + (t10 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,3) * t36 / 0.2e1) * t36 + (-t34 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t46 + Ifges(3,6) * t51 + Ifges(3,2) * t45 / 0.2e1) * t45;
T  = t6;
