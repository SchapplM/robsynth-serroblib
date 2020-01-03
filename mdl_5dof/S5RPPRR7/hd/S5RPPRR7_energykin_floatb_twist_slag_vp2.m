% Calculate kinetic energy for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:27
% EndTime: 2019-12-31 17:59:28
% DurationCPUTime: 0.70s
% Computational Cost: add. (1179->129), mult. (1702->176), div. (0->0), fcn. (1244->8), ass. (0->46)
t55 = pkin(2) + pkin(6);
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t35 = -t47 * V_base(4) + t50 * V_base(5);
t36 = t47 * V_base(5) + t50 * V_base(4);
t44 = sin(pkin(8));
t54 = cos(pkin(8));
t29 = -t35 * t54 + t36 * t44;
t30 = t44 * t35 + t36 * t54;
t33 = -pkin(1) * t35 + qJD(2) + V_base(3);
t52 = -qJ(3) * t30 + t33;
t12 = t29 * t55 + t52;
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t43 = V_base(6) + qJD(1);
t40 = pkin(5) * V_base(5) + V_base(1);
t41 = -pkin(5) * V_base(4) + V_base(2);
t31 = -t40 * t47 + t50 * t41;
t24 = pkin(1) * t43 - qJ(2) * t36 + t31;
t32 = t50 * t40 + t47 * t41;
t27 = qJ(2) * t35 + t32;
t18 = t24 * t54 - t44 * t27;
t53 = qJD(3) - t18;
t9 = t30 * pkin(3) - t43 * t55 + t53;
t6 = t49 * t12 + t46 * t9;
t19 = t44 * t24 + t54 * t27;
t14 = -t43 * qJ(3) - t19;
t5 = -t12 * t46 + t49 * t9;
t22 = t29 * t49 - t43 * t46;
t10 = -pkin(3) * t29 - t14;
t51 = V_base(3) ^ 2;
t48 = cos(qJ(5));
t45 = sin(qJ(5));
t28 = qJD(4) + t30;
t23 = t29 * t46 + t43 * t49;
t20 = qJD(5) - t22;
t17 = t23 * t48 + t28 * t45;
t16 = -t23 * t45 + t28 * t48;
t15 = pkin(2) * t29 + t52;
t13 = -t43 * pkin(2) + t53;
t7 = -pkin(4) * t22 - pkin(7) * t23 + t10;
t4 = pkin(7) * t28 + t6;
t3 = -pkin(4) * t28 - t5;
t2 = t4 * t48 + t45 * t7;
t1 = -t4 * t45 + t48 * t7;
t8 = m(2) * (t31 ^ 2 + t32 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t33 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t31 * mrSges(2,3) + Ifges(2,1) * t36 / 0.2e1) * t36 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t32 * mrSges(2,3) + Ifges(2,4) * t36 + Ifges(2,2) * t35 / 0.2e1) * t35 + (t10 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t20 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t10 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t28 + Ifges(5,2) * t22 / 0.2e1) * t22 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t20 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t13 * mrSges(4,1) + t33 * mrSges(3,2) - t18 * mrSges(3,3) - t15 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t30) * t30 + (t33 * mrSges(3,1) + t14 * mrSges(4,1) - t15 * mrSges(4,2) - t19 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t29 + (-Ifges(3,4) - Ifges(4,6)) * t30) * t29 + (t31 * mrSges(2,1) + t18 * mrSges(3,1) - t32 * mrSges(2,2) - t19 * mrSges(3,2) + t13 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(2,5) * t36 + Ifges(2,6) * t35 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t43 + (-Ifges(4,4) + Ifges(3,5)) * t30 + (Ifges(4,5) - Ifges(3,6)) * t29) * t43;
T = t8;
