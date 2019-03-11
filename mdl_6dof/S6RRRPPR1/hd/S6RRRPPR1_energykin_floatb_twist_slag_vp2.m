% Calculate kinetic energy for
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:07
% EndTime: 2019-03-09 15:21:08
% DurationCPUTime: 1.21s
% Computational Cost: add. (4265->156), mult. (5649->226), div. (0->0), fcn. (4648->12), ass. (0->59)
t69 = sin(qJ(1));
t73 = cos(qJ(1));
t54 = -t69 * V_base(4) + t73 * V_base(5);
t55 = t69 * V_base(5) + t73 * V_base(4);
t44 = -pkin(1) * t54 - pkin(7) * t55 + V_base(3);
t59 = V_base(5) * pkin(6) + V_base(1);
t60 = -V_base(4) * pkin(6) + V_base(2);
t51 = t73 * t59 + t69 * t60;
t62 = V_base(6) + qJD(1);
t47 = pkin(7) * t62 + t51;
t68 = sin(qJ(2));
t72 = cos(qJ(2));
t36 = t72 * t44 - t47 * t68;
t49 = t55 * t72 + t62 * t68;
t53 = qJD(2) - t54;
t33 = pkin(2) * t53 - pkin(8) * t49 + t36;
t37 = t68 * t44 + t72 * t47;
t48 = -t55 * t68 + t62 * t72;
t35 = pkin(8) * t48 + t37;
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t23 = t71 * t33 - t35 * t67;
t40 = t48 * t67 + t49 * t71;
t52 = qJD(3) + t53;
t19 = pkin(3) * t52 - qJ(4) * t40 + t23;
t24 = t67 * t33 + t71 * t35;
t39 = t48 * t71 - t49 * t67;
t22 = qJ(4) * t39 + t24;
t64 = sin(pkin(10));
t75 = cos(pkin(10));
t12 = t64 * t19 + t75 * t22;
t10 = qJ(5) * t52 + t12;
t28 = -t75 * t39 + t40 * t64;
t29 = t64 * t39 + t75 * t40;
t50 = -t69 * t59 + t60 * t73;
t46 = -pkin(1) * t62 - t50;
t41 = -pkin(2) * t48 + t46;
t30 = -pkin(3) * t39 + qJD(4) + t41;
t15 = pkin(4) * t28 - qJ(5) * t29 + t30;
t63 = sin(pkin(11));
t65 = cos(pkin(11));
t6 = t65 * t10 + t63 * t15;
t5 = -t10 * t63 + t65 * t15;
t11 = t75 * t19 - t64 * t22;
t9 = -t52 * pkin(4) + qJD(5) - t11;
t74 = V_base(3) ^ 2;
t70 = cos(qJ(6));
t66 = sin(qJ(6));
t27 = qJD(6) + t28;
t26 = t29 * t65 + t52 * t63;
t25 = -t29 * t63 + t52 * t65;
t17 = t25 * t66 + t26 * t70;
t16 = t25 * t70 - t26 * t66;
t7 = -t25 * pkin(5) + t9;
t4 = pkin(9) * t25 + t6;
t3 = pkin(5) * t28 - pkin(9) * t26 + t5;
t2 = t3 * t66 + t4 * t70;
t1 = t3 * t70 - t4 * t66;
t8 = (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t17 / 0.2e1) * t17 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t62 + Ifges(2,2) * t54 / 0.2e1) * t54 + m(2) * (t50 ^ 2 + t51 ^ 2 + t74) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t74) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t41 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t46 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t30 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (-t46 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t49 + Ifges(3,6) * t53 + Ifges(3,2) * t48 / 0.2e1) * t48 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t41 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,1) * t40 / 0.2e1) * t40 + (t30 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t29 + Ifges(6,5) * t26 - Ifges(5,6) * t52 + Ifges(6,6) * t25 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t28) * t28 + (t50 * mrSges(2,1) - t51 * mrSges(2,2) + Ifges(2,3) * t62 / 0.2e1) * t62 + (t46 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,5) * t53 + Ifges(3,1) * t49 / 0.2e1) * t49 + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,5) * t62 + Ifges(2,1) * t55 / 0.2e1) * t55 + (-t41 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,2) * t39 / 0.2e1) * t39 + (t30 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t29 / 0.2e1) * t29 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t27 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t23 * mrSges(4,1) + t11 * mrSges(5,1) - t24 * mrSges(4,2) - t12 * mrSges(5,2) + Ifges(4,5) * t40 + Ifges(5,5) * t29 + Ifges(4,6) * t39 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t52) * t52;
T  = t8;
