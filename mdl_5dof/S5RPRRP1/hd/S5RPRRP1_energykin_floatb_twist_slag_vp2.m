% Calculate kinetic energy for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:10
% DurationCPUTime: 0.64s
% Computational Cost: add. (989->125), mult. (1254->161), div. (0->0), fcn. (836->6), ass. (0->40)
t50 = pkin(1) + pkin(6);
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t32 = t43 * V_base(5) + t46 * V_base(4);
t40 = V_base(6) + qJD(1);
t36 = V_base(5) * pkin(5) + V_base(1);
t37 = -V_base(4) * pkin(5) + V_base(2);
t27 = -t43 * t36 + t37 * t46;
t48 = qJD(2) - t27;
t18 = pkin(2) * t32 - t50 * t40 + t48;
t31 = t43 * V_base(4) - t46 * V_base(5);
t49 = -qJ(2) * t32 + V_base(3);
t20 = t50 * t31 + t49;
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t11 = t45 * t18 - t20 * t42;
t26 = t31 * t42 + t40 * t45;
t30 = qJD(3) + t32;
t7 = pkin(3) * t30 - pkin(7) * t26 + t11;
t12 = t42 * t18 + t45 * t20;
t25 = t31 * t45 - t40 * t42;
t9 = pkin(7) * t25 + t12;
t4 = t41 * t7 + t44 * t9;
t28 = t46 * t36 + t43 * t37;
t24 = -t40 * qJ(2) - t28;
t3 = -t41 * t9 + t44 * t7;
t21 = -pkin(2) * t31 - t24;
t13 = -pkin(3) * t25 + t21;
t47 = V_base(3) ^ 2;
t29 = qJD(4) + t30;
t23 = -pkin(1) * t40 + t48;
t22 = pkin(1) * t31 + t49;
t15 = t25 * t41 + t26 * t44;
t14 = t25 * t44 - t26 * t41;
t10 = -pkin(4) * t14 + qJD(5) + t13;
t2 = qJ(5) * t14 + t4;
t1 = pkin(4) * t29 - qJ(5) * t15 + t3;
t5 = m(2) * (t27 ^ 2 + t28 ^ 2 + t47) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t47) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t21 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t30 / 0.2e1) * t30 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t21 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t30 + Ifges(4,1) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t21 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t26 + Ifges(4,6) * t30 + Ifges(4,2) * t25 / 0.2e1) * t25 + (t27 * mrSges(2,1) - t28 * mrSges(2,2) + t23 * mrSges(3,2) - t24 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t40) * t40 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29) * t29 + (t23 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t27 * mrSges(2,3) - t22 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t32 + (-Ifges(3,4) + Ifges(2,5)) * t40) * t32 + (t13 * mrSges(5,2) + t10 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t15 + (Ifges(5,5) + Ifges(6,5)) * t29) * t15 + (V_base(3) * mrSges(2,1) + t24 * mrSges(3,1) - t22 * mrSges(3,2) - t28 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t31 + (Ifges(3,5) - Ifges(2,6)) * t40 + (-Ifges(2,4) - Ifges(3,6)) * t32) * t31 + (-t13 * mrSges(5,1) - t10 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t14 + (Ifges(5,6) + Ifges(6,6)) * t29 + (Ifges(5,4) + Ifges(6,4)) * t15) * t14;
T = t5;
