% Calculate kinetic energy for
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:55
% EndTime: 2019-03-09 06:56:57
% DurationCPUTime: 1.30s
% Computational Cost: add. (3779->156), mult. (5477->228), div. (0->0), fcn. (4496->12), ass. (0->60)
t61 = V_base(5) * pkin(6) + V_base(1);
t62 = -V_base(4) * pkin(6) + V_base(2);
t70 = sin(qJ(1));
t75 = cos(qJ(1));
t52 = -t61 * t70 + t75 * t62;
t57 = t70 * V_base(5) + t75 * V_base(4);
t63 = V_base(6) + qJD(1);
t44 = pkin(1) * t63 - qJ(2) * t57 + t52;
t53 = t75 * t61 + t70 * t62;
t56 = -t70 * V_base(4) + t75 * V_base(5);
t48 = qJ(2) * t56 + t53;
t64 = sin(pkin(11));
t65 = cos(pkin(11));
t39 = t64 * t44 + t65 * t48;
t34 = pkin(7) * t63 + t39;
t50 = t56 * t65 - t64 * t57;
t51 = t56 * t64 + t57 * t65;
t54 = -pkin(1) * t56 + qJD(2) + V_base(3);
t37 = -pkin(2) * t50 - pkin(7) * t51 + t54;
t69 = sin(qJ(3));
t74 = cos(qJ(3));
t23 = -t34 * t69 + t74 * t37;
t43 = t51 * t74 + t63 * t69;
t49 = qJD(3) - t50;
t19 = pkin(3) * t49 - pkin(8) * t43 + t23;
t24 = t74 * t34 + t69 * t37;
t42 = -t51 * t69 + t63 * t74;
t22 = pkin(8) * t42 + t24;
t68 = sin(qJ(4));
t73 = cos(qJ(4));
t12 = t68 * t19 + t73 * t22;
t47 = qJD(4) + t49;
t10 = pkin(9) * t47 + t12;
t38 = t44 * t65 - t64 * t48;
t33 = -pkin(2) * t63 - t38;
t27 = -pkin(3) * t42 + t33;
t30 = t42 * t73 - t68 * t43;
t31 = t42 * t68 + t43 * t73;
t15 = -pkin(4) * t30 - pkin(9) * t31 + t27;
t67 = sin(qJ(5));
t72 = cos(qJ(5));
t6 = t72 * t10 + t67 * t15;
t5 = -t10 * t67 + t72 * t15;
t11 = t19 * t73 - t68 * t22;
t29 = qJD(5) - t30;
t9 = -pkin(4) * t47 - t11;
t76 = V_base(3) ^ 2;
t71 = cos(qJ(6));
t66 = sin(qJ(6));
t28 = qJD(6) + t29;
t26 = t31 * t72 + t47 * t67;
t25 = -t31 * t67 + t47 * t72;
t17 = t25 * t66 + t26 * t71;
t16 = t25 * t71 - t26 * t66;
t7 = -pkin(5) * t25 + t9;
t4 = pkin(10) * t25 + t6;
t3 = pkin(5) * t29 - pkin(10) * t26 + t5;
t2 = t3 * t66 + t4 * t71;
t1 = t3 * t71 - t4 * t66;
t8 = (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(2) * (t52 ^ 2 + t53 ^ 2 + t76) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t76) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t54 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t33 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t27 ^ 2) / 0.2e1 + (-t54 * mrSges(3,1) + t39 * mrSges(3,3) + Ifges(3,4) * t51 + Ifges(3,2) * t50 / 0.2e1) * t50 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + Ifges(4,3) * t49 / 0.2e1) * t49 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t33 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,5) * t49 + Ifges(4,1) * t43 / 0.2e1) * t43 + (t54 * mrSges(3,2) - t38 * mrSges(3,3) + Ifges(3,1) * t51 / 0.2e1) * t51 + (t52 * mrSges(2,1) + t38 * mrSges(3,1) - t53 * mrSges(2,2) - t39 * mrSges(3,2) + Ifges(2,5) * t57 + Ifges(3,5) * t51 + Ifges(2,6) * t56 + Ifges(3,6) * t50 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * t63) * t63 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + (-t33 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t49 + Ifges(4,2) * t42 / 0.2e1) * t42 + (V_base(3) * mrSges(2,2) - t52 * mrSges(2,3) + Ifges(2,1) * t57 / 0.2e1) * t57 + (-t27 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,6) * t47 + Ifges(5,2) * t30 / 0.2e1) * t30 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t29 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t28 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t17 / 0.2e1) * t17 + (t27 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t47 + Ifges(5,1) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(2,1) + t53 * mrSges(2,3) + Ifges(2,4) * t57 + Ifges(2,2) * t56 / 0.2e1) * t56 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t8;
