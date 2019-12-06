% Calculate kinetic energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:52
% EndTime: 2019-12-05 18:52:53
% DurationCPUTime: 0.92s
% Computational Cost: add. (1133->118), mult. (1612->179), div. (0->0), fcn. (1360->10), ass. (0->43)
t40 = sin(qJ(1));
t45 = cos(qJ(1));
t31 = -t40 * V_base(1) + t45 * V_base(2);
t35 = V_base(6) + qJD(1);
t26 = pkin(1) * t35 + t31;
t32 = t40 * V_base(2) + t45 * V_base(1);
t39 = sin(qJ(2));
t44 = cos(qJ(2));
t19 = t26 * t39 + t32 * t44;
t29 = -t40 * V_base(4) + t45 * V_base(5);
t27 = -pkin(1) * t29 + V_base(3);
t38 = sin(qJ(3));
t43 = cos(qJ(3));
t14 = t19 * t43 + t27 * t38;
t37 = sin(qJ(4));
t42 = cos(qJ(4));
t13 = -t19 * t38 + t43 * t27;
t30 = t40 * V_base(5) + t45 * V_base(4);
t22 = t29 * t44 - t39 * t30;
t21 = qJD(3) - t22;
t47 = pkin(2) * t21 + t13;
t3 = t14 * t37 - t42 * t47;
t49 = t3 ^ 2;
t17 = -t44 * t26 + t32 * t39;
t48 = t17 ^ 2;
t23 = t29 * t39 + t30 * t44;
t33 = qJD(2) + t35;
t15 = -t23 * t38 + t33 * t43;
t16 = t23 * t43 + t33 * t38;
t9 = t15 * t42 - t16 * t37;
t46 = V_base(3) ^ 2;
t41 = cos(qJ(5));
t36 = sin(qJ(5));
t20 = qJD(4) + t21;
t12 = -pkin(2) * t15 + t17;
t10 = t15 * t37 + t16 * t42;
t8 = qJD(5) - t9;
t7 = t10 * t41 + t20 * t36;
t6 = -t10 * t36 + t20 * t41;
t5 = t42 * t14 + t37 * t47;
t2 = t12 * t36 + t41 * t5;
t1 = t12 * t41 - t36 * t5;
t4 = m(2) * (t31 ^ 2 + t32 ^ 2 + t46) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t46) / 0.2e1 + m(3) * (t19 ^ 2 + t27 ^ 2 + t48) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t48) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t49) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t49) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t12 * mrSges(5,1) + t5 * mrSges(5,3) + Ifges(5,2) * t9 / 0.2e1) * t9 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t8 / 0.2e1) * t8 + (t31 * mrSges(2,1) - t32 * mrSges(2,2) + Ifges(2,3) * t35 / 0.2e1) * t35 + (-t17 * mrSges(3,1) - t19 * mrSges(3,2) + Ifges(3,3) * t33 / 0.2e1) * t33 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t8 + Ifges(6,1) * t7 / 0.2e1) * t7 + (V_base(3) * mrSges(2,2) - t31 * mrSges(2,3) + Ifges(2,5) * t35 + Ifges(2,1) * t30 / 0.2e1) * t30 + (t27 * mrSges(3,2) + t17 * mrSges(3,3) + Ifges(3,5) * t33 + Ifges(3,1) * t23 / 0.2e1) * t23 + (-t3 * mrSges(5,1) - t5 * mrSges(5,2) + Ifges(5,6) * t9 + Ifges(5,3) * t20 / 0.2e1) * t20 + (t17 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t21 + Ifges(4,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t7 + Ifges(6,6) * t8 + Ifges(6,2) * t6 / 0.2e1) * t6 + (-V_base(3) * mrSges(2,1) + t32 * mrSges(2,3) + Ifges(2,4) * t30 + Ifges(2,6) * t35 + Ifges(2,2) * t29 / 0.2e1) * t29 + (-t27 * mrSges(3,1) + t19 * mrSges(3,3) + Ifges(3,4) * t23 + Ifges(3,6) * t33 + Ifges(3,2) * t22 / 0.2e1) * t22 + (-t17 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t16 + Ifges(4,6) * t21 + Ifges(4,2) * t15 / 0.2e1) * t15 + (t12 * mrSges(5,2) + t3 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,5) * t20 + Ifges(5,1) * t10 / 0.2e1) * t10;
T = t4;
