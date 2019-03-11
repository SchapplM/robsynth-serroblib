% Calculate kinetic energy for
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:18
% EndTime: 2019-03-09 21:57:19
% DurationCPUTime: 1.32s
% Computational Cost: add. (4319->156), mult. (5649->228), div. (0->0), fcn. (4648->12), ass. (0->60)
t70 = sin(qJ(1));
t74 = cos(qJ(1));
t55 = -t70 * V_base(4) + t74 * V_base(5);
t56 = t70 * V_base(5) + t74 * V_base(4);
t44 = -pkin(1) * t55 - pkin(7) * t56 + V_base(3);
t60 = V_base(5) * pkin(6) + V_base(1);
t61 = -V_base(4) * pkin(6) + V_base(2);
t51 = t74 * t60 + t70 * t61;
t63 = V_base(6) + qJD(1);
t47 = pkin(7) * t63 + t51;
t69 = sin(qJ(2));
t73 = cos(qJ(2));
t36 = t73 * t44 - t47 * t69;
t49 = t56 * t73 + t63 * t69;
t54 = qJD(2) - t55;
t33 = pkin(2) * t54 - pkin(8) * t49 + t36;
t37 = t69 * t44 + t73 * t47;
t48 = -t56 * t69 + t63 * t73;
t35 = pkin(8) * t48 + t37;
t68 = sin(qJ(3));
t72 = cos(qJ(3));
t23 = t72 * t33 - t35 * t68;
t40 = t48 * t68 + t49 * t72;
t53 = qJD(3) + t54;
t19 = pkin(3) * t53 - pkin(9) * t40 + t23;
t24 = t68 * t33 + t72 * t35;
t39 = t48 * t72 - t49 * t68;
t22 = pkin(9) * t39 + t24;
t67 = sin(qJ(4));
t76 = cos(qJ(4));
t12 = t67 * t19 + t76 * t22;
t52 = qJD(4) + t53;
t10 = qJ(5) * t52 + t12;
t28 = -t76 * t39 + t40 * t67;
t29 = t67 * t39 + t76 * t40;
t50 = -t70 * t60 + t61 * t74;
t46 = -pkin(1) * t63 - t50;
t41 = -pkin(2) * t48 + t46;
t30 = -pkin(3) * t39 + t41;
t15 = pkin(4) * t28 - qJ(5) * t29 + t30;
t64 = sin(pkin(11));
t65 = cos(pkin(11));
t6 = t65 * t10 + t64 * t15;
t5 = -t10 * t64 + t65 * t15;
t11 = t76 * t19 - t67 * t22;
t9 = -t52 * pkin(4) + qJD(5) - t11;
t75 = V_base(3) ^ 2;
t71 = cos(qJ(6));
t66 = sin(qJ(6));
t27 = qJD(6) + t28;
t26 = t29 * t65 + t52 * t64;
t25 = -t29 * t64 + t52 * t65;
t17 = t25 * t66 + t26 * t71;
t16 = t25 * t71 - t26 * t66;
t7 = -t25 * pkin(5) + t9;
t4 = pkin(10) * t25 + t6;
t3 = pkin(5) * t28 - pkin(10) * t26 + t5;
t2 = t3 * t66 + t4 * t71;
t1 = t3 * t71 - t4 * t66;
t8 = (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t30 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t52 + Ifges(5,1) * t29 / 0.2e1) * t29 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t56 + Ifges(2,6) * t63 + Ifges(2,2) * t55 / 0.2e1) * t55 + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + Ifges(4,3) * t53 / 0.2e1) * t53 + (t50 * mrSges(2,1) - t51 * mrSges(2,2) + Ifges(2,3) * t63 / 0.2e1) * t63 + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,5) * t63 + Ifges(2,1) * t56 / 0.2e1) * t56 + m(2) * (t50 ^ 2 + t51 ^ 2 + t75) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t75) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t41 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t46 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(5) * (t11 ^ 2 + t12 ^ 2 + t30 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (-t46 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t49 + Ifges(3,6) * t54 + Ifges(3,2) * t48 / 0.2e1) * t48 + (-t41 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,6) * t53 + Ifges(4,2) * t39 / 0.2e1) * t39 + (t41 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,5) * t53 + Ifges(4,1) * t40 / 0.2e1) * t40 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t52 / 0.2e1) * t52 + (t30 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t29 + Ifges(6,5) * t26 - Ifges(5,6) * t52 + Ifges(6,6) * t25 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t28) * t28 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t17 / 0.2e1) * t17 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t26 / 0.2e1) * t26 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t27 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t54 / 0.2e1) * t54 + (t46 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,5) * t54 + Ifges(3,1) * t49 / 0.2e1) * t49 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T  = t8;
