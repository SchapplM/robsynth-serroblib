% Calculate kinetic energy for
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:09
% EndTime: 2019-03-09 01:37:10
% DurationCPUTime: 0.75s
% Computational Cost: add. (1363->152), mult. (1781->194), div. (0->0), fcn. (1192->8), ass. (0->52)
t52 = sin(qJ(1));
t61 = cos(qJ(1));
t39 = t52 * V_base(5) + t61 * V_base(4);
t47 = V_base(6) + qJD(1);
t43 = V_base(5) * pkin(6) + V_base(1);
t44 = -V_base(4) * pkin(6) + V_base(2);
t35 = -t52 * t43 + t44 * t61;
t57 = qJD(2) - t35;
t56 = (-pkin(1) - qJ(3)) * t47 + t57;
t60 = pkin(2) + qJ(4);
t20 = t39 * t60 + t56;
t38 = t52 * V_base(4) - t61 * V_base(5);
t36 = t61 * t43 + t52 * t44;
t34 = -t47 * qJ(2) - t36;
t59 = qJD(3) - t34;
t21 = pkin(3) * t47 - t38 * t60 + t59;
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t12 = t49 * t20 + t48 * t21;
t10 = pkin(7) * t47 + t12;
t37 = t38 * qJ(3);
t58 = pkin(1) * t38 + V_base(3);
t22 = qJD(4) + t37 + (-pkin(3) - qJ(2)) * t39 + t58;
t30 = -t38 * t48 + t39 * t49;
t31 = t38 * t49 + t39 * t48;
t14 = -pkin(4) * t30 - pkin(7) * t31 + t22;
t51 = sin(qJ(5));
t54 = cos(qJ(5));
t6 = t54 * t10 + t51 * t14;
t11 = -t48 * t20 + t21 * t49;
t5 = -t10 * t51 + t14 * t54;
t26 = -t31 * t51 + t47 * t54;
t9 = -pkin(4) * t47 - t11;
t32 = -qJ(2) * t39 + t58;
t55 = V_base(3) ^ 2;
t53 = cos(qJ(6));
t50 = sin(qJ(6));
t33 = -t47 * pkin(1) + t57;
t29 = qJD(5) - t30;
t28 = -pkin(2) * t38 + t59;
t27 = t31 * t54 + t47 * t51;
t25 = -t32 - t37;
t24 = qJD(6) - t26;
t23 = t39 * pkin(2) + t56;
t16 = t27 * t53 + t29 * t50;
t15 = -t27 * t50 + t29 * t53;
t7 = -pkin(5) * t26 - pkin(8) * t27 + t9;
t4 = pkin(8) * t29 + t6;
t3 = -pkin(5) * t29 - t5;
t2 = t4 * t53 + t50 * t7;
t1 = -t4 * t50 + t53 * t7;
t8 = m(2) * (t35 ^ 2 + t36 ^ 2 + t55) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t55) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t25 ^ 2 + t28 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t22 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t31 / 0.2e1) * t31 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t22 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,2) * t30 / 0.2e1) * t30 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t29 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t24 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t33 * mrSges(3,1) + t25 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t23 * mrSges(4,2) - t35 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t39) * t39 + (V_base(3) * mrSges(2,1) + t34 * mrSges(3,1) - t32 * mrSges(3,2) + t28 * mrSges(4,2) - t36 * mrSges(2,3) - t25 * mrSges(4,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t38 + (-Ifges(2,4) + Ifges(4,5) - Ifges(3,6)) * t39) * t38 + (t35 * mrSges(2,1) + t28 * mrSges(4,1) + t11 * mrSges(5,1) - t36 * mrSges(2,2) + t33 * mrSges(3,2) - t12 * mrSges(5,2) - t34 * mrSges(3,3) - t23 * mrSges(4,3) + Ifges(5,5) * t31 + Ifges(5,6) * t30 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t47 + (-Ifges(3,4) + Ifges(2,5) - Ifges(4,6)) * t39 + (-Ifges(4,4) + Ifges(3,5) - Ifges(2,6)) * t38) * t47;
T  = t8;
