% Calculate kinetic energy for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:56
% EndTime: 2019-12-05 17:37:57
% DurationCPUTime: 0.58s
% Computational Cost: add. (757->127), mult. (968->161), div. (0->0), fcn. (588->6), ass. (0->40)
t45 = sin(qJ(1));
t51 = cos(qJ(1));
t30 = t45 * V_base(4) - t51 * V_base(5);
t42 = V_base(6) + qJD(1);
t36 = V_base(5) * pkin(5) + V_base(1);
t37 = -V_base(4) * pkin(5) + V_base(2);
t25 = t51 * t36 + t45 * t37;
t21 = -t42 * qJ(2) - t25;
t50 = qJD(3) - t21;
t11 = -pkin(6) * t42 + (-pkin(2) - pkin(3)) * t30 + t50;
t27 = t30 * qJ(3);
t31 = t45 * V_base(5) + t51 * V_base(4);
t49 = pkin(1) * t30 + V_base(3);
t13 = t27 + (pkin(6) - qJ(2)) * t31 + t49;
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t6 = t44 * t11 + t47 * t13;
t24 = -t45 * t36 + t51 * t37;
t5 = t47 * t11 - t13 * t44;
t20 = -t42 * pkin(1) + qJD(2) - t24;
t16 = t31 * pkin(2) - t42 * qJ(3) + t20;
t29 = qJD(4) - t30;
t19 = -qJ(2) * t31 + t49;
t10 = -pkin(3) * t31 - t16;
t48 = V_base(3) ^ 2;
t46 = cos(qJ(5));
t43 = sin(qJ(5));
t26 = qJD(5) + t29;
t23 = t31 * t44 + t42 * t47;
t22 = t31 * t47 - t42 * t44;
t18 = -pkin(2) * t30 + t50;
t17 = t19 + t27;
t15 = t22 * t43 + t23 * t46;
t14 = t22 * t46 - t23 * t43;
t7 = -pkin(4) * t22 + t10;
t4 = pkin(7) * t22 + t6;
t3 = pkin(4) * t29 - pkin(7) * t23 + t5;
t2 = t3 * t43 + t4 * t46;
t1 = t3 * t46 - t4 * t43;
t8 = m(2) * (t24 ^ 2 + t25 ^ 2 + t48) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t48) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(3) * (t19 ^ 2 + t20 ^ 2 + t21 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t29 / 0.2e1) * t29 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t10 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t29 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t26 + Ifges(6,1) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t10 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t29 + Ifges(5,2) * t22 / 0.2e1) * t22 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t26 + Ifges(6,2) * t14 / 0.2e1) * t14 + (t24 * mrSges(2,1) - t25 * mrSges(2,2) + t20 * mrSges(3,2) + t18 * mrSges(4,2) - t21 * mrSges(3,3) - t16 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t42) * t42 + (t20 * mrSges(3,1) + t16 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t17 * mrSges(4,2) - t24 * mrSges(2,3) - t19 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t31 + (-Ifges(3,4) + Ifges(2,5) + Ifges(4,5)) * t42) * t31 + (V_base(3) * mrSges(2,1) + t21 * mrSges(3,1) - t18 * mrSges(4,1) - t19 * mrSges(3,2) - t25 * mrSges(2,3) + t17 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t30 + (Ifges(4,4) + Ifges(3,5) - Ifges(2,6)) * t42 + (-Ifges(2,4) - Ifges(3,6) + Ifges(4,6)) * t31) * t30;
T = t8;
