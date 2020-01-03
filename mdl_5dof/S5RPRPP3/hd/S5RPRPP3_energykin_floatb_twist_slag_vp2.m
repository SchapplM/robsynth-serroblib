% Calculate kinetic energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:16
% DurationCPUTime: 0.58s
% Computational Cost: add. (1099->125), mult. (1468->159), div. (0->0), fcn. (1056->6), ass. (0->39)
t43 = sin(qJ(1));
t44 = cos(qJ(1));
t31 = t43 * V_base(4) - t44 * V_base(5);
t32 = t43 * V_base(5) + t44 * V_base(4);
t21 = pkin(1) * t31 - qJ(2) * t32 + V_base(3);
t36 = V_base(5) * pkin(5) + V_base(1);
t37 = -V_base(4) * pkin(5) + V_base(2);
t29 = t44 * t36 + t43 * t37;
t39 = V_base(6) + qJD(1);
t25 = qJ(2) * t39 + t29;
t40 = sin(pkin(7));
t41 = cos(pkin(7));
t14 = t41 * t21 - t25 * t40;
t27 = t32 * t41 + t39 * t40;
t10 = pkin(2) * t31 - pkin(6) * t27 + t14;
t15 = t40 * t21 + t41 * t25;
t26 = -t32 * t40 + t39 * t41;
t13 = pkin(6) * t26 + t15;
t42 = sin(qJ(3));
t49 = cos(qJ(3));
t7 = t42 * t10 + t49 * t13;
t48 = pkin(3) + qJ(5);
t28 = -t43 * t36 + t37 * t44;
t30 = qJD(3) + t31;
t5 = -qJ(4) * t30 - t7;
t6 = t49 * t10 - t42 * t13;
t47 = qJD(4) - t6;
t23 = -pkin(1) * t39 + qJD(2) - t28;
t18 = -pkin(2) * t26 + t23;
t17 = t42 * t26 + t49 * t27;
t46 = -qJ(4) * t17 + t18;
t45 = V_base(3) ^ 2;
t16 = -t49 * t26 + t27 * t42;
t8 = pkin(3) * t16 + t46;
t4 = -t30 * pkin(3) + t47;
t3 = t48 * t16 + t46;
t2 = -pkin(4) * t16 + qJD(5) - t5;
t1 = t17 * pkin(4) - t48 * t30 + t47;
t9 = m(2) * (t28 ^ 2 + t29 ^ 2 + t45) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t45) / 0.2e1 + m(4) * (t18 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(3) * (t14 ^ 2 + t15 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t28 * mrSges(2,1) - t29 * mrSges(2,2) + Ifges(2,3) * t39 / 0.2e1) * t39 + (t23 * mrSges(3,2) - t14 * mrSges(3,3) + Ifges(3,1) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t28 * mrSges(2,3) + Ifges(2,5) * t39 + Ifges(2,1) * t32 / 0.2e1) * t32 + (-t23 * mrSges(3,1) + t15 * mrSges(3,3) + Ifges(3,4) * t27 + Ifges(3,2) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (V_base(3) * mrSges(2,1) + t14 * mrSges(3,1) - t15 * mrSges(3,2) - t29 * mrSges(2,3) - Ifges(2,4) * t32 + Ifges(3,5) * t27 - Ifges(2,6) * t39 + Ifges(3,6) * t26 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t31) * t31 + (t6 * mrSges(4,1) - t7 * mrSges(4,2) + t4 * mrSges(5,2) + t2 * mrSges(6,2) - t5 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t30) * t30 + (t4 * mrSges(5,1) + t1 * mrSges(6,1) + t18 * mrSges(4,2) - t3 * mrSges(6,2) - t6 * mrSges(4,3) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t17 + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t30) * t17 + (t18 * mrSges(4,1) + t5 * mrSges(5,1) - t2 * mrSges(6,1) - t8 * mrSges(5,2) - t7 * mrSges(4,3) + t3 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t16 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t30 + (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t17) * t16;
T = t9;
