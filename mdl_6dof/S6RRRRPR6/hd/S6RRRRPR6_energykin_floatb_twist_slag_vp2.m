% Calculate kinetic energy for
% S6RRRRPR6
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:32
% EndTime: 2019-03-09 22:18:34
% DurationCPUTime: 1.33s
% Computational Cost: add. (4699->156), mult. (6061->228), div. (0->0), fcn. (4972->12), ass. (0->60)
t69 = sin(qJ(1));
t74 = cos(qJ(1));
t55 = -t69 * V_base(4) + t74 * V_base(5);
t56 = t69 * V_base(5) + t74 * V_base(4);
t41 = -pkin(1) * t55 - pkin(7) * t56 + V_base(3);
t60 = V_base(5) * pkin(6) + V_base(1);
t61 = -V_base(4) * pkin(6) + V_base(2);
t52 = t74 * t60 + t69 * t61;
t62 = V_base(6) + qJD(1);
t47 = pkin(7) * t62 + t52;
t68 = sin(qJ(2));
t73 = cos(qJ(2));
t37 = t68 * t41 + t73 * t47;
t54 = qJD(2) - t55;
t34 = pkin(8) * t54 + t37;
t51 = -t69 * t60 + t61 * t74;
t46 = -pkin(1) * t62 - t51;
t49 = -t68 * t56 + t62 * t73;
t50 = t56 * t73 + t62 * t68;
t35 = -pkin(2) * t49 - pkin(8) * t50 + t46;
t67 = sin(qJ(3));
t72 = cos(qJ(3));
t25 = -t34 * t67 + t72 * t35;
t39 = t50 * t72 + t54 * t67;
t48 = qJD(3) - t49;
t19 = pkin(3) * t48 - pkin(9) * t39 + t25;
t26 = t72 * t34 + t67 * t35;
t38 = -t50 * t67 + t54 * t72;
t22 = pkin(9) * t38 + t26;
t66 = sin(qJ(4));
t71 = cos(qJ(4));
t13 = t66 * t19 + t71 * t22;
t28 = t38 * t71 - t39 * t66;
t11 = qJ(5) * t28 + t13;
t63 = sin(pkin(11));
t64 = cos(pkin(11));
t12 = t71 * t19 - t22 * t66;
t29 = t38 * t66 + t39 * t71;
t45 = qJD(4) + t48;
t9 = pkin(4) * t45 - qJ(5) * t29 + t12;
t6 = t64 * t11 + t63 * t9;
t5 = -t11 * t63 + t64 * t9;
t36 = t41 * t73 - t68 * t47;
t33 = -pkin(2) * t54 - t36;
t27 = -pkin(3) * t38 + t33;
t20 = -pkin(4) * t28 + qJD(5) + t27;
t75 = V_base(3) ^ 2;
t70 = cos(qJ(6));
t65 = sin(qJ(6));
t44 = qJD(6) + t45;
t24 = t28 * t63 + t29 * t64;
t23 = t28 * t64 - t29 * t63;
t16 = t23 * t65 + t24 * t70;
t15 = t23 * t70 - t24 * t65;
t14 = -pkin(5) * t23 + t20;
t4 = pkin(10) * t23 + t6;
t3 = pkin(5) * t45 - pkin(10) * t24 + t5;
t2 = t3 * t65 + t4 * t70;
t1 = t3 * t70 - t4 * t65;
t7 = (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t20 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t24 + Ifges(6,2) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t56 + Ifges(2,6) * t62 + Ifges(2,2) * t55 / 0.2e1) * t55 + (V_base(3) * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,5) * t62 + Ifges(2,1) * t56 / 0.2e1) * t56 + (-t46 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t50 + Ifges(3,6) * t54 + Ifges(3,2) * t49 / 0.2e1) * t49 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t13 * mrSges(5,2) - t6 * mrSges(6,2) + Ifges(5,5) * t29 + Ifges(6,5) * t24 + Ifges(5,6) * t28 + Ifges(6,6) * t23 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t45) * t45 + m(2) * (t51 ^ 2 + t52 ^ 2 + t75) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t75) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t25 ^ 2 + t26 ^ 2 + t33 ^ 2) / 0.2e1 + m(6) * (t20 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) / 0.2e1 + (-t33 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t39 + Ifges(4,6) * t48 + Ifges(4,2) * t38 / 0.2e1) * t38 + (t25 * mrSges(4,1) - t26 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t46 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,5) * t54 + Ifges(3,1) * t50 / 0.2e1) * t50 + (t14 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t44 + Ifges(7,1) * t16 / 0.2e1) * t16 + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t29 / 0.2e1) * t29 + (-t14 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t44 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t54 / 0.2e1) * t54 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t20 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t24 / 0.2e1) * t24 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t44 / 0.2e1) * t44 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t51 * mrSges(2,1) - t52 * mrSges(2,2) + Ifges(2,3) * t62 / 0.2e1) * t62 + (t33 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t39 / 0.2e1) * t39 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,2) * t28 / 0.2e1) * t28;
T  = t7;
