% Calculate kinetic energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:53
% EndTime: 2019-12-05 15:02:54
% DurationCPUTime: 0.75s
% Computational Cost: add. (1143->129), mult. (1774->175), div. (0->0), fcn. (1316->8), ass. (0->45)
t55 = pkin(3) + pkin(6);
t54 = cos(qJ(3));
t46 = sin(pkin(7));
t48 = cos(pkin(7));
t36 = t46 * V_base(4) - t48 * V_base(5);
t37 = t46 * V_base(5) + t48 * V_base(4);
t44 = V_base(3) + qJD(1);
t26 = pkin(1) * t36 - qJ(2) * t37 + t44;
t41 = qJ(1) * V_base(5) + V_base(1);
t42 = -qJ(1) * V_base(4) + V_base(2);
t34 = t48 * t41 + t46 * t42;
t30 = V_base(6) * qJ(2) + t34;
t45 = sin(pkin(8));
t47 = cos(pkin(8));
t18 = t47 * t26 - t30 * t45;
t32 = t37 * t47 + t45 * V_base(6);
t12 = pkin(2) * t36 - pkin(5) * t32 + t18;
t19 = t45 * t26 + t47 * t30;
t31 = -t37 * t45 + t47 * V_base(6);
t15 = pkin(5) * t31 + t19;
t50 = sin(qJ(3));
t9 = t50 * t12 + t54 * t15;
t33 = -t46 * t41 + t42 * t48;
t35 = qJD(3) + t36;
t6 = -qJ(4) * t35 - t9;
t8 = t12 * t54 - t50 * t15;
t53 = qJD(4) - t8;
t22 = t50 * t31 + t32 * t54;
t28 = -V_base(6) * pkin(1) + qJD(2) - t33;
t23 = -pkin(2) * t31 + t28;
t52 = -qJ(4) * t22 + t23;
t51 = cos(qJ(5));
t49 = sin(qJ(5));
t21 = -t31 * t54 + t32 * t50;
t20 = qJD(5) + t22;
t17 = t21 * t49 + t35 * t51;
t16 = t21 * t51 - t35 * t49;
t10 = pkin(3) * t21 + t52;
t7 = t21 * t55 + t52;
t5 = -t35 * pkin(3) + t53;
t4 = -pkin(4) * t21 - t6;
t3 = t22 * pkin(4) - t35 * t55 + t53;
t2 = t3 * t49 + t51 * t7;
t1 = t3 * t51 - t49 * t7;
t11 = m(2) * (t33 ^ 2 + t34 ^ 2 + t44 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t44 * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,1) * t37 / 0.2e1) * t37 + (t28 * mrSges(3,2) - t18 * mrSges(3,3) + Ifges(3,1) * t32 / 0.2e1) * t32 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t28 * mrSges(3,1) + t19 * mrSges(3,3) + Ifges(3,4) * t32 + Ifges(3,2) * t31 / 0.2e1) * t31 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t20 + Ifges(6,1) * t17 / 0.2e1) * t17 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t20 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t8 * mrSges(4,1) - t9 * mrSges(4,2) + t5 * mrSges(5,2) - t6 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t35) * t35 + (t5 * mrSges(5,1) + t23 * mrSges(4,2) - t8 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t22 + (-Ifges(5,4) + Ifges(4,5)) * t35) * t22 + (V_base(2) * mrSges(1,1) + t33 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t34 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t37 + Ifges(1,6) * V_base(5) + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t44 * mrSges(2,1) + t18 * mrSges(3,1) - t19 * mrSges(3,2) - t34 * mrSges(2,3) - Ifges(2,4) * t37 + Ifges(3,5) * t32 - Ifges(2,6) * V_base(6) + Ifges(3,6) * t31 + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t36) * t36 + (t23 * mrSges(4,1) + t6 * mrSges(5,1) - t10 * mrSges(5,2) - t9 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t21 + (Ifges(5,5) - Ifges(4,6)) * t35 + (-Ifges(4,4) - Ifges(5,6)) * t22) * t21;
T = t11;
