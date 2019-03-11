% Calculate kinetic energy for
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:35
% EndTime: 2019-03-08 20:04:37
% DurationCPUTime: 1.19s
% Computational Cost: add. (4409->157), mult. (7649->220), div. (0->0), fcn. (6476->12), ass. (0->58)
t63 = sin(pkin(10));
t66 = cos(pkin(10));
t53 = -t63 * V_base(4) + t66 * V_base(5);
t64 = sin(pkin(6));
t67 = cos(pkin(6));
t74 = t53 * t67 + t64 * V_base(6);
t59 = V_base(5) * qJ(1) + V_base(1);
t60 = -V_base(4) * qJ(1) + V_base(2);
t50 = -t59 * t63 + t66 * t60;
t54 = t63 * V_base(5) + t66 * V_base(4);
t80 = pkin(7) * t54;
t45 = V_base(6) * pkin(1) - t67 * t80 + t50;
t61 = V_base(3) + qJD(1);
t48 = -pkin(1) * t53 - t64 * t80 + t61;
t81 = t45 * t67 + t48 * t64;
t51 = t66 * t59 + t63 * t60;
t41 = t74 * pkin(7) + t51;
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t33 = -t70 * t41 + t81 * t73;
t49 = -t53 * t64 + t67 * V_base(6) + qJD(2);
t28 = -t49 * pkin(2) + qJD(3) - t33;
t44 = t54 * t73 + t74 * t70;
t62 = sin(pkin(11));
t65 = cos(pkin(11));
t36 = -t44 * t62 + t49 * t65;
t21 = -t36 * pkin(3) + t28;
t37 = t44 * t65 + t49 * t62;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t31 = t36 * t72 - t37 * t69;
t32 = t36 * t69 + t37 * t72;
t13 = -t31 * pkin(4) - t32 * pkin(9) + t21;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t35 = -t45 * t64 + t67 * t48;
t43 = t54 * t70 - t74 * t73;
t26 = pkin(2) * t43 - qJ(3) * t44 + t35;
t34 = t73 * t41 + t81 * t70;
t29 = qJ(3) * t49 + t34;
t19 = t65 * t26 - t29 * t62;
t15 = pkin(3) * t43 - pkin(8) * t37 + t19;
t20 = t62 * t26 + t65 * t29;
t18 = pkin(8) * t36 + t20;
t10 = t69 * t15 + t72 * t18;
t42 = qJD(4) + t43;
t8 = pkin(9) * t42 + t10;
t4 = t68 * t13 + t71 * t8;
t3 = t71 * t13 - t68 * t8;
t9 = t15 * t72 - t69 * t18;
t7 = -pkin(4) * t42 - t9;
t30 = qJD(5) - t31;
t23 = t32 * t71 + t42 * t68;
t22 = -t32 * t68 + t42 * t71;
t5 = -pkin(5) * t22 + qJD(6) + t7;
t2 = qJ(6) * t22 + t4;
t1 = pkin(5) * t30 - qJ(6) * t23 + t3;
t6 = m(3) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t20 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t21 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t50 ^ 2 + t51 ^ 2 + t61 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t61 * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,1) * t54 / 0.2e1) * t54 + (t33 * mrSges(3,1) - t34 * mrSges(3,2) + Ifges(3,3) * t49 / 0.2e1) * t49 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t42 / 0.2e1) * t42 + (t28 * mrSges(4,2) - t19 * mrSges(4,3) + Ifges(4,1) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t61 * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,2) * t53 / 0.2e1) * t53 + (t35 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,5) * t49 + Ifges(3,1) * t44 / 0.2e1) * t44 + (-t28 * mrSges(4,1) + t20 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t21 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t42 + Ifges(5,1) * t32 / 0.2e1) * t32 + (-t21 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t32 + Ifges(5,6) * t42 + Ifges(5,2) * t31 / 0.2e1) * t31 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t30) * t30 + (t7 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t23 + (Ifges(6,5) + Ifges(7,5)) * t30) * t23 + (V_base(2) * mrSges(1,1) + t50 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t51 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t54 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t53 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t35 * mrSges(3,1) + t19 * mrSges(4,1) - t20 * mrSges(4,2) - t34 * mrSges(3,3) - Ifges(3,4) * t44 + Ifges(4,5) * t37 - Ifges(3,6) * t49 + Ifges(4,6) * t36 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t43) * t43 + (-t7 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t22 + (Ifges(6,6) + Ifges(7,6)) * t30 + (Ifges(6,4) + Ifges(7,4)) * t23) * t22;
T  = t6;
