% Calculate kinetic energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:20
% EndTime: 2019-07-18 13:27:21
% DurationCPUTime: 0.46s
% Computational Cost: add. (421->96), mult. (599->135), div. (0->0), fcn. (372->6), ass. (0->30)
t25 = -V_base(2) + qJD(1);
t20 = -V_base(6) * qJ(1) - V_base(1);
t21 = -V_base(4) * qJ(1) + V_base(3);
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t13 = t31 * t20 - t21 * t28;
t23 = -V_base(5) + qJD(2);
t12 = pkin(1) * t23 + t13;
t14 = t28 * t20 + t21 * t31;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t6 = t30 * t12 - t14 * t27;
t22 = qJD(3) + t23;
t16 = -t28 * V_base(6) - t31 * V_base(4);
t15 = -pkin(1) * t16 + t25;
t29 = cos(qJ(4));
t26 = sin(qJ(4));
t24 = t25 ^ 2;
t19 = qJD(4) + t22;
t17 = -t28 * V_base(4) + t31 * V_base(6);
t10 = t16 * t27 + t17 * t30;
t9 = t16 * t30 - t17 * t27;
t8 = -pkin(2) * t9 + t15;
t7 = t12 * t27 + t14 * t30;
t5 = pkin(2) * t22 + t6;
t4 = t10 * t29 + t26 * t9;
t3 = -t10 * t26 + t29 * t9;
t2 = t26 * t5 + t29 * t7;
t1 = -t26 * t7 + t29 * t5;
t11 = m(2) * (t20 ^ 2 + t21 ^ 2 + t24) / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t24) / 0.2e1 + m(4) * (t15 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-t15 * mrSges(4,1) + t7 * mrSges(4,3) + Ifges(4,2) * t9 / 0.2e1) * t9 + (t8 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t4 / 0.2e1) * t4 + (t13 * mrSges(3,1) - t14 * mrSges(3,2) + Ifges(3,3) * t23 / 0.2e1) * t23 + (-t8 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t4 + Ifges(5,2) * t3 / 0.2e1) * t3 + (t6 * mrSges(4,1) - t7 * mrSges(4,2) + Ifges(4,6) * t9 + Ifges(4,3) * t22 / 0.2e1) * t22 + (t25 * mrSges(3,2) - t13 * mrSges(3,3) + Ifges(3,5) * t23 + Ifges(3,1) * t17 / 0.2e1) * t17 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t4 + Ifges(5,6) * t3 + Ifges(5,3) * t19 / 0.2e1) * t19 + (-t25 * mrSges(3,1) + t14 * mrSges(3,3) + Ifges(3,4) * t17 + Ifges(3,6) * t23 + Ifges(3,2) * t16 / 0.2e1) * t16 + (t15 * mrSges(4,2) - t6 * mrSges(4,3) + Ifges(4,4) * t9 + Ifges(4,5) * t22 + Ifges(4,1) * t10 / 0.2e1) * t10 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + t25 * mrSges(2,2) - t20 * mrSges(2,3) + (Ifges(1,3) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) - t20 * mrSges(2,1) + t21 * mrSges(2,2) + V_base(1) * mrSges(1,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5) + (-Ifges(2,5) + Ifges(1,6)) * V_base(6)) * V_base(5) + (t25 * mrSges(2,1) + V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) - t21 * mrSges(2,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(4) + (-Ifges(2,4) + Ifges(1,5)) * V_base(6) + (Ifges(1,4) + Ifges(2,6)) * V_base(5)) * V_base(4);
T  = t11;
