% Calculate kinetic energy for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:06
% EndTime: 2019-03-09 01:39:07
% DurationCPUTime: 1.07s
% Computational Cost: add. (3689->156), mult. (5477->224), div. (0->0), fcn. (4496->12), ass. (0->58)
t59 = V_base(5) * pkin(6) + V_base(1);
t60 = -V_base(4) * pkin(6) + V_base(2);
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t50 = -t59 * t70 + t72 * t60;
t55 = t70 * V_base(5) + t72 * V_base(4);
t61 = V_base(6) + qJD(1);
t43 = pkin(1) * t61 - qJ(2) * t55 + t50;
t51 = t72 * t59 + t70 * t60;
t54 = -t70 * V_base(4) + t72 * V_base(5);
t46 = qJ(2) * t54 + t51;
t64 = sin(pkin(9));
t67 = cos(pkin(9));
t38 = t64 * t43 + t67 * t46;
t33 = qJ(3) * t61 + t38;
t48 = -t67 * t54 + t55 * t64;
t49 = t54 * t64 + t55 * t67;
t52 = -pkin(1) * t54 + qJD(2) + V_base(3);
t36 = pkin(2) * t48 - qJ(3) * t49 + t52;
t63 = sin(pkin(10));
t66 = cos(pkin(10));
t23 = -t33 * t63 + t66 * t36;
t41 = t49 * t66 + t61 * t63;
t19 = pkin(3) * t48 - pkin(7) * t41 + t23;
t24 = t66 * t33 + t63 * t36;
t40 = -t49 * t63 + t61 * t66;
t22 = pkin(7) * t40 + t24;
t69 = sin(qJ(4));
t74 = cos(qJ(4));
t12 = t69 * t19 + t74 * t22;
t47 = qJD(4) + t48;
t10 = qJ(5) * t47 + t12;
t37 = t43 * t67 - t64 * t46;
t32 = -pkin(2) * t61 + qJD(3) - t37;
t27 = -pkin(3) * t40 + t32;
t29 = -t74 * t40 + t41 * t69;
t30 = t69 * t40 + t41 * t74;
t15 = pkin(4) * t29 - qJ(5) * t30 + t27;
t62 = sin(pkin(11));
t65 = cos(pkin(11));
t6 = t65 * t10 + t62 * t15;
t5 = -t10 * t62 + t65 * t15;
t11 = t19 * t74 - t69 * t22;
t9 = -t47 * pkin(4) + qJD(5) - t11;
t73 = V_base(3) ^ 2;
t71 = cos(qJ(6));
t68 = sin(qJ(6));
t28 = qJD(6) + t29;
t26 = t30 * t65 + t47 * t62;
t25 = -t30 * t62 + t47 * t65;
t17 = t25 * t68 + t26 * t71;
t16 = t25 * t71 - t26 * t68;
t7 = -t25 * pkin(5) + t9;
t4 = pkin(8) * t25 + t6;
t3 = pkin(5) * t29 - pkin(8) * t26 + t5;
t2 = t3 * t68 + t4 * t71;
t1 = t3 * t71 - t4 * t68;
t8 = (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,2) * t25 / 0.2e1) * t25 + (t32 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,1) * t41 / 0.2e1) * t41 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t17 / 0.2e1) * t17 + m(2) * (t50 ^ 2 + t51 ^ 2 + t73) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t73) / 0.2e1 + m(3) * (t37 ^ 2 + t38 ^ 2 + t52 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t52 * mrSges(3,2) - t37 * mrSges(3,3) + Ifges(3,1) * t49 / 0.2e1) * t49 + (t27 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t47 + Ifges(5,1) * t30 / 0.2e1) * t30 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(4) * (t23 ^ 2 + t24 ^ 2 + t32 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t27 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (t52 * mrSges(3,1) + t23 * mrSges(4,1) - t24 * mrSges(4,2) - t38 * mrSges(3,3) - Ifges(3,4) * t49 + Ifges(4,5) * t41 - Ifges(3,6) * t61 + Ifges(4,6) * t40 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t48) * t48 + (t50 * mrSges(2,1) + t37 * mrSges(3,1) - t51 * mrSges(2,2) - t38 * mrSges(3,2) + Ifges(2,5) * t55 + Ifges(3,5) * t49 + Ifges(2,6) * t54 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t61) * t61 + (t27 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t30 + Ifges(6,5) * t26 - Ifges(5,6) * t47 + Ifges(6,6) * t25 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t29) * t29 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t28 + Ifges(7,2) * t16 / 0.2e1) * t16 + (-t32 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t41 + Ifges(4,2) * t40 / 0.2e1) * t40 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,1) * t55 / 0.2e1) * t55;
T  = t8;
