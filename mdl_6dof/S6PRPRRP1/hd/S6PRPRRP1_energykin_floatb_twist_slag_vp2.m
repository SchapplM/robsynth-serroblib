% Calculate kinetic energy for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:07
% EndTime: 2019-03-08 19:55:08
% DurationCPUTime: 1.18s
% Computational Cost: add. (4303->157), mult. (7751->220), div. (0->0), fcn. (6578->12), ass. (0->58)
t58 = V_base(5) * qJ(1) + V_base(1);
t59 = -V_base(4) * qJ(1) + V_base(2);
t62 = sin(pkin(10));
t65 = cos(pkin(10));
t51 = -t58 * t62 + t65 * t59;
t66 = cos(pkin(6));
t54 = t62 * V_base(5) + t65 * V_base(4);
t76 = pkin(7) * t54;
t45 = V_base(6) * pkin(1) - t66 * t76 + t51;
t53 = -t62 * V_base(4) + t65 * V_base(5);
t60 = V_base(3) + qJD(1);
t63 = sin(pkin(6));
t49 = -pkin(1) * t53 - t63 * t76 + t60;
t77 = t45 * t66 + t49 * t63;
t52 = t65 * t58 + t62 * t59;
t73 = t53 * t66 + t63 * V_base(6);
t42 = t73 * pkin(7) + t52;
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t29 = -t42 * t69 + t77 * t72;
t44 = t54 * t72 + t73 * t69;
t50 = -t53 * t63 + t66 * V_base(6) + qJD(2);
t25 = pkin(2) * t50 - qJ(3) * t44 + t29;
t30 = t72 * t42 + t77 * t69;
t43 = -t69 * t54 + t73 * t72;
t28 = qJ(3) * t43 + t30;
t61 = sin(pkin(11));
t64 = cos(pkin(11));
t18 = t25 * t64 - t61 * t28;
t16 = -pkin(3) * t50 - t18;
t37 = t43 * t61 + t44 * t64;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t33 = -t37 * t68 + t50 * t71;
t34 = t37 * t71 + t50 * t68;
t13 = -pkin(4) * t33 - pkin(9) * t34 + t16;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t19 = t61 * t25 + t64 * t28;
t17 = pkin(8) * t50 + t19;
t38 = -t45 * t63 + t66 * t49;
t31 = -pkin(2) * t43 + qJD(3) + t38;
t36 = t43 * t64 - t44 * t61;
t21 = -pkin(3) * t36 - pkin(8) * t37 + t31;
t10 = t71 * t17 + t68 * t21;
t35 = qJD(4) - t36;
t8 = pkin(9) * t35 + t10;
t4 = t67 * t13 + t70 * t8;
t3 = t70 * t13 - t67 * t8;
t9 = -t68 * t17 + t21 * t71;
t7 = -pkin(4) * t35 - t9;
t32 = qJD(5) - t33;
t23 = t34 * t70 + t35 * t67;
t22 = -t34 * t67 + t35 * t70;
t5 = -pkin(5) * t22 + qJD(6) + t7;
t2 = qJ(6) * t22 + t4;
t1 = pkin(5) * t32 - qJ(6) * t23 + t3;
t6 = m(3) * (t29 ^ 2 + t30 ^ 2 + t38 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t16 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t51 ^ 2 + t52 ^ 2 + t60 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t60 * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,1) * t54 / 0.2e1) * t54 + (t38 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,1) * t44 / 0.2e1) * t44 + (t31 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t37 / 0.2e1) * t37 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t60 * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,2) * t53 / 0.2e1) * t53 + (-t38 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,2) * t43 / 0.2e1) * t43 + (-t31 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t16 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t35 + Ifges(5,1) * t34 / 0.2e1) * t34 + (-t16 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t34 + Ifges(5,6) * t35 + Ifges(5,2) * t33 / 0.2e1) * t33 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t32) * t32 + (t7 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t23 + (Ifges(6,5) + Ifges(7,5)) * t32) * t23 + (V_base(2) * mrSges(1,1) + t51 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t52 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t54 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t53 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t29 * mrSges(3,1) + t18 * mrSges(4,1) - t30 * mrSges(3,2) - t19 * mrSges(4,2) + Ifges(3,5) * t44 + Ifges(4,5) * t37 + Ifges(3,6) * t43 + Ifges(4,6) * t36 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t50) * t50 + (-t7 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t22 + (Ifges(6,6) + Ifges(7,6)) * t32 + (Ifges(6,4) + Ifges(7,4)) * t23) * t22;
T  = t6;
