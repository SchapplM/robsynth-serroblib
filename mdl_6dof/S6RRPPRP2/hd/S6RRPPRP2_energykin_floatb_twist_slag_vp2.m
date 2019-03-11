% Calculate kinetic energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:33
% EndTime: 2019-03-09 08:29:34
% DurationCPUTime: 0.83s
% Computational Cost: add. (2037->149), mult. (2665->194), div. (0->0), fcn. (2048->8), ass. (0->49)
t62 = pkin(3) + pkin(8);
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t44 = t54 * V_base(5) + t57 * V_base(4);
t50 = V_base(6) + qJD(1);
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t37 = -t44 * t53 + t50 * t56;
t38 = t44 * t56 + t50 * t53;
t51 = sin(pkin(9));
t61 = cos(pkin(9));
t27 = -t37 * t61 + t38 * t51;
t28 = t37 * t51 + t38 * t61;
t48 = pkin(6) * V_base(5) + V_base(1);
t49 = -pkin(6) * V_base(4) + V_base(2);
t39 = -t48 * t54 + t49 * t57;
t35 = -pkin(1) * t50 - t39;
t29 = -pkin(2) * t37 + qJD(3) + t35;
t59 = -qJ(4) * t28 + t29;
t11 = t27 * t62 + t59;
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t43 = -t54 * V_base(4) + t57 * V_base(5);
t42 = qJD(2) - t43;
t32 = -pkin(1) * t43 - pkin(7) * t44 + V_base(3);
t40 = t48 * t57 + t49 * t54;
t36 = pkin(7) * t50 + t40;
t24 = t32 * t56 - t36 * t53;
t18 = pkin(2) * t42 - qJ(3) * t38 + t24;
t25 = t32 * t53 + t36 * t56;
t21 = qJ(3) * t37 + t25;
t14 = t18 * t61 - t21 * t51;
t60 = qJD(4) - t14;
t8 = pkin(4) * t28 - t42 * t62 + t60;
t4 = t11 * t55 + t52 * t8;
t15 = t18 * t51 + t21 * t61;
t13 = -qJ(4) * t42 - t15;
t3 = -t11 * t52 + t55 * t8;
t9 = -pkin(4) * t27 - t13;
t58 = V_base(3) ^ 2;
t26 = qJD(5) + t28;
t23 = t27 * t52 + t42 * t55;
t22 = t27 * t55 - t42 * t52;
t16 = pkin(3) * t27 + t59;
t12 = -pkin(3) * t42 + t60;
t5 = -pkin(5) * t22 + qJD(6) + t9;
t2 = qJ(6) * t22 + t4;
t1 = pkin(5) * t26 - qJ(6) * t23 + t3;
t6 = m(2) * (t39 ^ 2 + t40 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t58) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t35 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t50 / 0.2e1) * t50 + (t35 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,1) * t38 / 0.2e1) * t38 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,5) * t50 + Ifges(2,1) * t44 / 0.2e1) * t44 + (-t35 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,2) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t50 + Ifges(2,2) * t43 / 0.2e1) * t43 + (t12 * mrSges(5,1) + t29 * mrSges(4,2) - t14 * mrSges(4,3) - t16 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t28) * t28 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t26) * t26 + (t29 * mrSges(4,1) + t13 * mrSges(5,1) - t16 * mrSges(5,2) - t15 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t27 + (-Ifges(4,4) - Ifges(5,6)) * t28) * t27 + (t9 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t23 + (Ifges(6,5) + Ifges(7,5)) * t26) * t23 + (-t9 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t22 + (Ifges(6,6) + Ifges(7,6)) * t26 + (Ifges(6,4) + Ifges(7,4)) * t23) * t22 + (t24 * mrSges(3,1) + t14 * mrSges(4,1) - t25 * mrSges(3,2) - t15 * mrSges(4,2) + t12 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(3,5) * t38 + Ifges(3,6) * t37 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t42 + (-Ifges(5,4) + Ifges(4,5)) * t28 + (Ifges(5,5) - Ifges(4,6)) * t27) * t42;
T  = t6;
