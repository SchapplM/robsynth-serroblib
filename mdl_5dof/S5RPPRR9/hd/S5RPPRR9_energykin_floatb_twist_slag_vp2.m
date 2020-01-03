% Calculate kinetic energy for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:15
% EndTime: 2019-12-31 18:02:16
% DurationCPUTime: 0.71s
% Computational Cost: add. (1115->129), mult. (1500->176), div. (0->0), fcn. (1036->8), ass. (0->44)
t50 = sin(qJ(1));
t55 = cos(qJ(1));
t37 = t50 * V_base(5) + t55 * V_base(4);
t45 = V_base(6) + qJD(1);
t41 = pkin(5) * V_base(5) + V_base(1);
t42 = -pkin(5) * V_base(4) + V_base(2);
t32 = -t50 * t41 + t42 * t55;
t54 = qJD(2) - t32;
t18 = -t37 * qJ(3) + (-pkin(1) - pkin(2)) * t45 + t54;
t33 = t55 * t41 + t50 * t42;
t31 = t45 * qJ(2) + t33;
t36 = t50 * V_base(4) - t55 * V_base(5);
t25 = qJ(3) * t36 + t31;
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t14 = t46 * t18 + t47 * t25;
t12 = -pkin(6) * t45 + t14;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t29 = t36 * pkin(1) - t37 * qJ(2) + V_base(3);
t19 = -pkin(2) * t36 + qJD(3) - t29;
t27 = t36 * t47 - t37 * t46;
t28 = t36 * t46 + t37 * t47;
t9 = -pkin(3) * t27 - pkin(6) * t28 + t19;
t6 = t52 * t12 + t49 * t9;
t13 = t18 * t47 - t46 * t25;
t5 = -t12 * t49 + t52 * t9;
t23 = -t28 * t49 - t45 * t52;
t11 = pkin(3) * t45 - t13;
t53 = V_base(3) ^ 2;
t51 = cos(qJ(5));
t48 = sin(qJ(5));
t30 = -t45 * pkin(1) + t54;
t26 = qJD(4) - t27;
t24 = t28 * t52 - t45 * t49;
t20 = qJD(5) - t23;
t16 = t24 * t51 + t26 * t48;
t15 = -t24 * t48 + t26 * t51;
t7 = -pkin(4) * t23 - pkin(7) * t24 + t11;
t4 = pkin(7) * t26 + t6;
t3 = -pkin(4) * t26 - t5;
t2 = t4 * t51 + t48 * t7;
t1 = -t4 * t48 + t51 * t7;
t8 = m(2) * (t32 ^ 2 + t33 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t19 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,1) * t28 / 0.2e1) * t28 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t26 / 0.2e1) * t26 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t19 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t28 + Ifges(4,2) * t27 / 0.2e1) * t27 + (t11 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t26 + Ifges(5,1) * t24 / 0.2e1) * t24 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t20 + Ifges(6,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t11 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,6) * t26 + Ifges(5,2) * t23 / 0.2e1) * t23 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t20 + Ifges(6,2) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(2,2) + t30 * mrSges(3,2) - t32 * mrSges(2,3) - t29 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t37) * t37 + (V_base(3) * mrSges(2,1) + t29 * mrSges(3,1) - t31 * mrSges(3,2) - t33 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t36 + (-Ifges(2,4) + Ifges(3,5)) * t37) * t36 + (t32 * mrSges(2,1) - t30 * mrSges(3,1) - t13 * mrSges(4,1) - t33 * mrSges(2,2) + t14 * mrSges(4,2) + t31 * mrSges(3,3) - Ifges(4,5) * t28 - Ifges(4,6) * t27 + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t45 + (Ifges(3,4) + Ifges(2,5)) * t37 + (-Ifges(2,6) + Ifges(3,6)) * t36) * t45;
T = t8;
