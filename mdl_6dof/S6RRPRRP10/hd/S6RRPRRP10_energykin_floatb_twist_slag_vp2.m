% Calculate kinetic energy for
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:25
% EndTime: 2019-03-09 12:36:27
% DurationCPUTime: 1.27s
% Computational Cost: add. (5057->157), mult. (7627->221), div. (0->0), fcn. (6454->12), ass. (0->59)
t69 = sin(qJ(1));
t72 = cos(qJ(1));
t53 = -t69 * V_base(4) + t72 * V_base(5);
t61 = V_base(6) + qJD(1);
t63 = sin(pkin(6));
t65 = cos(pkin(6));
t74 = t53 * t65 + t61 * t63;
t59 = pkin(7) * V_base(5) + V_base(1);
t60 = -V_base(4) * pkin(7) + V_base(2);
t50 = -t59 * t69 + t60 * t72;
t54 = t69 * V_base(5) + t72 * V_base(4);
t81 = pkin(8) * t54;
t45 = pkin(1) * t61 - t65 * t81 + t50;
t48 = -pkin(1) * t53 - t63 * t81 + V_base(3);
t82 = t45 * t65 + t48 * t63;
t51 = t59 * t72 + t60 * t69;
t42 = pkin(8) * t74 + t51;
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t32 = -t42 * t68 + t71 * t82;
t49 = -t53 * t63 + t61 * t65 + qJD(2);
t27 = -pkin(2) * t49 + qJD(3) - t32;
t44 = t54 * t71 + t68 * t74;
t62 = sin(pkin(11));
t64 = cos(pkin(11));
t35 = -t44 * t62 + t49 * t64;
t20 = -pkin(3) * t35 + t27;
t36 = t44 * t64 + t49 * t62;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t30 = t35 * t70 - t36 * t67;
t31 = t35 * t67 + t36 * t70;
t12 = -pkin(4) * t30 - pkin(10) * t31 + t20;
t66 = sin(qJ(5));
t34 = -t45 * t63 + t48 * t65;
t43 = t54 * t68 - t71 * t74;
t25 = pkin(2) * t43 - qJ(3) * t44 + t34;
t33 = t71 * t42 + t68 * t82;
t29 = qJ(3) * t49 + t33;
t18 = t64 * t25 - t29 * t62;
t14 = pkin(3) * t43 - pkin(9) * t36 + t18;
t19 = t62 * t25 + t64 * t29;
t17 = pkin(9) * t35 + t19;
t10 = t14 * t67 + t17 * t70;
t41 = qJD(4) + t43;
t8 = pkin(10) * t41 + t10;
t80 = cos(qJ(5));
t4 = t12 * t66 + t8 * t80;
t9 = t14 * t70 - t17 * t67;
t7 = -pkin(4) * t41 - t9;
t3 = t12 * t80 - t66 * t8;
t73 = V_base(3) ^ 2;
t28 = qJD(5) - t30;
t22 = t31 * t80 + t41 * t66;
t21 = t31 * t66 - t41 * t80;
t5 = pkin(5) * t21 - qJ(6) * t22 + t7;
t2 = qJ(6) * t28 + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t20 ^ 2 + t9 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t27 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t73) / 0.2e1 + m(2) * (t50 ^ 2 + t51 ^ 2 + t73) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t50 * mrSges(2,1) - t51 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (t32 * mrSges(3,1) - t33 * mrSges(3,2) + Ifges(3,3) * t49 / 0.2e1) * t49 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (t27 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t36 / 0.2e1) * t36 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t54 / 0.2e1) * t54 + (t34 * mrSges(3,2) - t32 * mrSges(3,3) + Ifges(3,5) * t49 + Ifges(3,1) * t44 / 0.2e1) * t44 + (-t27 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t36 + Ifges(4,2) * t35 / 0.2e1) * t35 + (t20 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t31 / 0.2e1) * t31 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,6) * t61 + Ifges(2,2) * t53 / 0.2e1) * t53 + (-t20 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,6) * t41 + Ifges(5,2) * t30 / 0.2e1) * t30 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t28) * t28 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t28) * t22 + (t34 * mrSges(3,1) + t18 * mrSges(4,1) - t19 * mrSges(4,2) - t33 * mrSges(3,3) - Ifges(3,4) * t44 + Ifges(4,5) * t36 - Ifges(3,6) * t49 + Ifges(4,6) * t35 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t43) * t43 + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t28 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21;
T  = t6;
