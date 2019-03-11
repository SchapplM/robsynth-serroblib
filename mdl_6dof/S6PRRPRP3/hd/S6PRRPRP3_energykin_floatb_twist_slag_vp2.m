% Calculate kinetic energy for
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:33:55
% EndTime: 2019-03-08 21:33:57
% DurationCPUTime: 1.16s
% Computational Cost: add. (4433->157), mult. (7627->220), div. (0->0), fcn. (6450->12), ass. (0->58)
t58 = V_base(5) * qJ(1) + V_base(1);
t59 = -V_base(4) * qJ(1) + V_base(2);
t62 = sin(pkin(10));
t65 = cos(pkin(10));
t51 = -t58 * t62 + t65 * t59;
t66 = cos(pkin(6));
t54 = t62 * V_base(5) + t65 * V_base(4);
t78 = pkin(7) * t54;
t45 = V_base(6) * pkin(1) - t66 * t78 + t51;
t53 = -t62 * V_base(4) + t65 * V_base(5);
t60 = V_base(3) + qJD(1);
t63 = sin(pkin(6));
t48 = -pkin(1) * t53 - t63 * t78 + t60;
t79 = t45 * t66 + t48 * t63;
t52 = t65 * t58 + t62 * t59;
t72 = t53 * t66 + t63 * V_base(6);
t41 = t72 * pkin(7) + t52;
t69 = sin(qJ(2));
t71 = cos(qJ(2));
t29 = -t69 * t41 + t79 * t71;
t43 = -t69 * t54 + t72 * t71;
t67 = sin(qJ(5));
t34 = -t45 * t63 + t66 * t48;
t44 = t54 * t71 + t72 * t69;
t24 = -pkin(2) * t43 - pkin(8) * t44 + t34;
t30 = t71 * t41 + t79 * t69;
t50 = -t53 * t63 + t66 * V_base(6) + qJD(2);
t28 = pkin(8) * t50 + t30;
t68 = sin(qJ(3));
t70 = cos(qJ(3));
t17 = t68 * t24 + t70 * t28;
t42 = qJD(3) - t43;
t15 = qJ(4) * t42 + t17;
t27 = -t50 * pkin(2) - t29;
t36 = t44 * t68 - t70 * t50;
t37 = t44 * t70 + t50 * t68;
t20 = t36 * pkin(3) - t37 * qJ(4) + t27;
t61 = sin(pkin(11));
t64 = cos(pkin(11));
t10 = -t15 * t61 + t64 * t20;
t33 = t37 * t64 + t42 * t61;
t7 = pkin(4) * t36 - pkin(9) * t33 + t10;
t77 = cos(qJ(5));
t11 = t64 * t15 + t61 * t20;
t32 = -t37 * t61 + t42 * t64;
t9 = pkin(9) * t32 + t11;
t4 = t67 * t7 + t77 * t9;
t16 = t24 * t70 - t68 * t28;
t3 = -t67 * t9 + t77 * t7;
t14 = -pkin(3) * t42 + qJD(4) - t16;
t12 = -pkin(4) * t32 + t14;
t35 = qJD(5) + t36;
t22 = t67 * t32 + t77 * t33;
t21 = -t77 * t32 + t33 * t67;
t5 = pkin(5) * t21 - qJ(6) * t22 + t12;
t2 = qJ(6) * t35 + t4;
t1 = -t35 * pkin(5) + qJD(6) - t3;
t6 = m(3) * (t29 ^ 2 + t30 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t27 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t51 ^ 2 + t52 ^ 2 + t60 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t60 * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,1) * t54 / 0.2e1) * t54 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t50 / 0.2e1) * t50 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t42 / 0.2e1) * t42 + (t14 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t33 / 0.2e1) * t33 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t60 * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,2) * t53 / 0.2e1) * t53 + (t34 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t50 + Ifges(3,1) * t44 / 0.2e1) * t44 + (t27 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t42 + Ifges(4,1) * t37 / 0.2e1) * t37 + (-t14 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,2) * t32 / 0.2e1) * t32 + (-t34 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,6) * t50 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t35) * t22 + (V_base(2) * mrSges(1,1) + t51 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t52 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t54 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t53 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t27 * mrSges(4,1) + t10 * mrSges(5,1) - t11 * mrSges(5,2) - t17 * mrSges(4,3) - Ifges(4,4) * t37 + Ifges(5,5) * t33 - Ifges(4,6) * t42 + Ifges(5,6) * t32 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t36) * t36 + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t35 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21;
T  = t6;
