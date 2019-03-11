% Calculate kinetic energy for
% S6PRRRRP4
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:01
% EndTime: 2019-03-09 00:13:02
% DurationCPUTime: 1.26s
% Computational Cost: add. (4457->157), mult. (7627->222), div. (0->0), fcn. (6450->12), ass. (0->59)
t59 = qJ(1) * V_base(5) + V_base(1);
t60 = -qJ(1) * V_base(4) + V_base(2);
t62 = sin(pkin(11));
t64 = cos(pkin(11));
t52 = -t59 * t62 + t60 * t64;
t65 = cos(pkin(6));
t55 = t62 * V_base(5) + t64 * V_base(4);
t79 = pkin(7) * t55;
t47 = V_base(6) * pkin(1) - t65 * t79 + t52;
t54 = -t62 * V_base(4) + t64 * V_base(5);
t61 = V_base(3) + qJD(1);
t63 = sin(pkin(6));
t50 = -pkin(1) * t54 - t63 * t79 + t61;
t80 = t47 * t65 + t50 * t63;
t53 = t59 * t64 + t60 * t62;
t73 = t54 * t65 + t63 * V_base(6);
t43 = pkin(7) * t73 + t53;
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t29 = -t43 * t69 + t72 * t80;
t45 = -t69 * t55 + t72 * t73;
t66 = sin(qJ(5));
t34 = -t47 * t63 + t50 * t65;
t46 = t55 * t72 + t69 * t73;
t24 = -pkin(2) * t45 - pkin(8) * t46 + t34;
t30 = t72 * t43 + t69 * t80;
t51 = -t54 * t63 + t65 * V_base(6) + qJD(2);
t28 = pkin(8) * t51 + t30;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t17 = t24 * t68 + t28 * t71;
t44 = qJD(3) - t45;
t15 = pkin(9) * t44 + t17;
t27 = -t51 * pkin(2) - t29;
t37 = -t46 * t68 + t51 * t71;
t38 = t46 * t71 + t51 * t68;
t20 = -t37 * pkin(3) - t38 * pkin(9) + t27;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t10 = -t15 * t67 + t20 * t70;
t33 = t38 * t70 + t44 * t67;
t36 = qJD(4) - t37;
t7 = pkin(4) * t36 - pkin(10) * t33 + t10;
t78 = cos(qJ(5));
t11 = t15 * t70 + t20 * t67;
t32 = -t38 * t67 + t44 * t70;
t9 = pkin(10) * t32 + t11;
t4 = t66 * t7 + t78 * t9;
t16 = t24 * t71 - t28 * t68;
t14 = -pkin(3) * t44 - t16;
t3 = -t66 * t9 + t7 * t78;
t12 = -pkin(4) * t32 + t14;
t35 = qJD(5) + t36;
t22 = t32 * t66 + t33 * t78;
t21 = -t32 * t78 + t33 * t66;
t5 = pkin(5) * t21 - qJ(6) * t22 + t12;
t2 = qJ(6) * t35 + t4;
t1 = -pkin(5) * t35 + qJD(6) - t3;
t6 = m(3) * (t29 ^ 2 + t30 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t27 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t52 ^ 2 + t53 ^ 2 + t61 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t61 * mrSges(2,2) - t52 * mrSges(2,3) + Ifges(2,1) * t55 / 0.2e1) * t55 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t51 / 0.2e1) * t51 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t44 / 0.2e1) * t44 + (t10 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,3) * t36 / 0.2e1) * t36 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t61 * mrSges(2,1) + t53 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t34 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t51 + Ifges(3,1) * t46 / 0.2e1) * t46 + (t27 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t44 + Ifges(4,1) * t38 / 0.2e1) * t38 + (t14 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,5) * t36 + Ifges(5,1) * t33 / 0.2e1) * t33 + (-t34 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t46 + Ifges(3,6) * t51 + Ifges(3,2) * t45 / 0.2e1) * t45 + (-t27 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t38 + Ifges(4,6) * t44 + Ifges(4,2) * t37 / 0.2e1) * t37 + (-t14 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,6) * t36 + Ifges(5,2) * t32 / 0.2e1) * t32 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t35) * t35 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t35) * t22 + (V_base(2) * mrSges(1,1) + t52 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t53 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t55 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t54 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t35 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21;
T  = t6;
