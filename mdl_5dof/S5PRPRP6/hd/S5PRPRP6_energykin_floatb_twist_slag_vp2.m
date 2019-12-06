% Calculate kinetic energy for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:16
% EndTime: 2019-12-05 15:40:16
% DurationCPUTime: 0.60s
% Computational Cost: add. (861->125), mult. (1238->160), div. (0->0), fcn. (848->6), ass. (0->39)
t49 = pkin(2) + pkin(6);
t41 = sin(pkin(7));
t42 = cos(pkin(7));
t33 = t41 * V_base(5) + t42 * V_base(4);
t44 = sin(qJ(2));
t48 = cos(qJ(2));
t26 = t33 * t44 - t48 * V_base(6);
t37 = V_base(5) * qJ(1) + V_base(1);
t38 = -V_base(4) * qJ(1) + V_base(2);
t28 = -t41 * t37 + t38 * t42;
t22 = -V_base(6) * pkin(1) - t28;
t27 = t48 * t33 + t44 * V_base(6);
t45 = -qJ(3) * t27 + t22;
t10 = t49 * t26 + t45;
t43 = sin(qJ(4));
t47 = cos(qJ(4));
t32 = -t41 * V_base(4) + t42 * V_base(5);
t31 = qJD(2) - t32;
t40 = V_base(3) + qJD(1);
t19 = -pkin(1) * t32 - pkin(5) * t33 + t40;
t29 = t42 * t37 + t41 * t38;
t23 = V_base(6) * pkin(5) + t29;
t14 = t48 * t19 - t44 * t23;
t46 = qJD(3) - t14;
t7 = t27 * pkin(3) - t49 * t31 + t46;
t4 = t47 * t10 + t43 * t7;
t15 = t44 * t19 + t48 * t23;
t12 = -t31 * qJ(3) - t15;
t8 = -pkin(3) * t26 - t12;
t3 = -t43 * t10 + t47 * t7;
t25 = qJD(4) + t27;
t17 = t43 * t26 + t47 * t31;
t16 = -t47 * t26 + t31 * t43;
t13 = pkin(2) * t26 + t45;
t11 = -t31 * pkin(2) + t46;
t5 = pkin(4) * t16 - qJ(5) * t17 + t8;
t2 = qJ(5) * t25 + t4;
t1 = -t25 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t28 ^ 2 + t29 ^ 2 + t40 ^ 2) / 0.2e1 + m(3) * (t14 ^ 2 + t15 ^ 2 + t22 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t13 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t40 * mrSges(2,2) - t28 * mrSges(2,3) + Ifges(2,1) * t33 / 0.2e1) * t33 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t40 * mrSges(2,1) + t29 * mrSges(2,3) + Ifges(2,4) * t33 + Ifges(2,2) * t32 / 0.2e1) * t32 + (t14 * mrSges(3,1) - t15 * mrSges(3,2) + t11 * mrSges(4,2) - t12 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t31) * t31 + (t3 * mrSges(5,1) - t1 * mrSges(6,1) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t25) * t25 + (t11 * mrSges(4,1) + t22 * mrSges(3,2) - t14 * mrSges(3,3) - t13 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t27 + (-Ifges(4,4) + Ifges(3,5)) * t31) * t27 + (t8 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t17 + (Ifges(6,4) + Ifges(5,5)) * t25) * t17 + (V_base(2) * mrSges(1,1) + t28 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t29 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t33 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t32 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t22 * mrSges(3,1) + t12 * mrSges(4,1) - t13 * mrSges(4,2) - t15 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t26 + (Ifges(4,5) - Ifges(3,6)) * t31 + (-Ifges(3,4) - Ifges(4,6)) * t27) * t26 + (t8 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t16 + (-Ifges(5,6) + Ifges(6,6)) * t25 + (-Ifges(5,4) + Ifges(6,5)) * t17) * t16;
T = t6;
