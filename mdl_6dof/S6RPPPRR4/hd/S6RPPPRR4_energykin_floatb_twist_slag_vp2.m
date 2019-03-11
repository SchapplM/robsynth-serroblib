% Calculate kinetic energy for
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:14
% EndTime: 2019-03-09 01:35:14
% DurationCPUTime: 0.80s
% Computational Cost: add. (1421->150), mult. (1919->194), div. (0->0), fcn. (1340->8), ass. (0->50)
t63 = pkin(3) + pkin(7);
t54 = sin(qJ(1));
t62 = cos(qJ(1));
t40 = t54 * V_base(4) - t62 * V_base(5);
t41 = t54 * V_base(5) + t62 * V_base(4);
t51 = sin(pkin(9));
t61 = cos(pkin(9));
t30 = -t61 * t40 + t41 * t51;
t32 = t40 * pkin(1) - t41 * qJ(2) + V_base(3);
t22 = -pkin(2) * t40 + qJD(3) - t32;
t31 = t51 * t40 + t61 * t41;
t58 = -qJ(4) * t31 + t22;
t10 = t63 * t30 + t58;
t50 = V_base(6) + qJD(1);
t45 = V_base(5) * pkin(6) + V_base(1);
t46 = -V_base(4) * pkin(6) + V_base(2);
t35 = -t54 * t45 + t62 * t46;
t60 = qJD(2) - t35;
t21 = -t41 * qJ(3) + (-pkin(1) - pkin(2)) * t50 + t60;
t36 = t62 * t45 + t54 * t46;
t34 = t50 * qJ(2) + t36;
t28 = qJ(3) * t40 + t34;
t16 = t61 * t21 - t51 * t28;
t59 = qJD(4) - t16;
t11 = t31 * pkin(4) + t63 * t50 + t59;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t6 = t56 * t10 + t53 * t11;
t17 = t51 * t21 + t61 * t28;
t15 = t50 * qJ(4) - t17;
t5 = -t10 * t53 + t11 * t56;
t26 = t30 * t56 + t50 * t53;
t12 = -pkin(4) * t30 - t15;
t57 = V_base(3) ^ 2;
t55 = cos(qJ(6));
t52 = sin(qJ(6));
t33 = -t50 * pkin(1) + t60;
t29 = qJD(5) + t31;
t27 = t30 * t53 - t50 * t56;
t23 = qJD(6) - t26;
t19 = t27 * t55 + t29 * t52;
t18 = -t27 * t52 + t29 * t55;
t14 = t50 * pkin(3) + t59;
t13 = pkin(3) * t30 + t58;
t7 = -pkin(5) * t26 - pkin(8) * t27 + t12;
t4 = pkin(8) * t29 + t6;
t3 = -pkin(5) * t29 - t5;
t2 = t4 * t55 + t52 * t7;
t1 = -t4 * t52 + t55 * t7;
t8 = m(2) * (t35 ^ 2 + t36 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t22 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t12 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t23 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t12 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t29 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t23 + Ifges(7,2) * t18 / 0.2e1) * t18 + (V_base(3) * mrSges(2,2) + t33 * mrSges(3,2) - t35 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t41) * t41 + (t14 * mrSges(5,1) + t22 * mrSges(4,2) - t16 * mrSges(4,3) - t13 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t31) * t31 + (V_base(3) * mrSges(2,1) + t32 * mrSges(3,1) - t34 * mrSges(3,2) - t36 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t40 + (-Ifges(2,4) + Ifges(3,5)) * t41) * t40 + (t22 * mrSges(4,1) + t15 * mrSges(5,1) - t13 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t30 + (-Ifges(4,4) - Ifges(5,6)) * t31) * t30 + (t35 * mrSges(2,1) - t33 * mrSges(3,1) - t16 * mrSges(4,1) - t36 * mrSges(2,2) + t17 * mrSges(4,2) - t14 * mrSges(5,2) + t34 * mrSges(3,3) + t15 * mrSges(5,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t50 + (Ifges(3,4) + Ifges(2,5)) * t41 + (-Ifges(2,6) + Ifges(3,6)) * t40 + (Ifges(5,4) - Ifges(4,5)) * t31 + (-Ifges(5,5) + Ifges(4,6)) * t30) * t50;
T  = t8;
