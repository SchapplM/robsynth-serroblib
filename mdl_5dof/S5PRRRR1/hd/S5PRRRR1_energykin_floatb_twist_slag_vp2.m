% Calculate kinetic energy for
% S5PRRRR1
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:17
% EndTime: 2019-07-18 13:28:17
% DurationCPUTime: 0.68s
% Computational Cost: add. (749->116), mult. (954->165), div. (0->0), fcn. (636->8), ass. (0->38)
t32 = V_base(3) + qJD(1);
t27 = -V_base(5) * pkin(1) + t32;
t28 = V_base(5) * qJ(1) + V_base(1);
t36 = sin(qJ(2));
t40 = cos(qJ(2));
t19 = t36 * t27 + t28 * t40;
t29 = -V_base(4) * qJ(1) + V_base(2);
t26 = -V_base(6) * pkin(1) - t29;
t35 = sin(qJ(3));
t39 = cos(qJ(3));
t14 = t19 * t39 + t26 * t35;
t34 = sin(qJ(4));
t38 = cos(qJ(4));
t13 = -t19 * t35 + t39 * t26;
t23 = -t36 * V_base(4) + t40 * V_base(6);
t22 = qJD(3) - t23;
t41 = pkin(2) * t22 + t13;
t3 = t14 * t34 - t38 * t41;
t43 = t3 ^ 2;
t17 = -t40 * t27 + t28 * t36;
t42 = t17 ^ 2;
t24 = t36 * V_base(6) + t40 * V_base(4);
t31 = -V_base(5) + qJD(2);
t15 = -t24 * t35 + t31 * t39;
t16 = t24 * t39 + t31 * t35;
t10 = t15 * t38 - t16 * t34;
t37 = cos(qJ(5));
t33 = sin(qJ(5));
t20 = qJD(4) + t22;
t12 = -pkin(2) * t15 + t17;
t11 = t15 * t34 + t16 * t38;
t9 = qJD(5) - t10;
t7 = t11 * t37 + t20 * t33;
t6 = -t11 * t33 + t20 * t37;
t5 = t38 * t14 + t34 * t41;
t2 = t12 * t33 + t37 * t5;
t1 = t12 * t37 - t33 * t5;
t4 = m(2) * (t28 ^ 2 + t29 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t19 ^ 2 + t26 ^ 2 + t42) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t42) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t43) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t43) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t9 / 0.2e1) * t9 + (-t17 * mrSges(3,1) - t19 * mrSges(3,2) + Ifges(3,3) * t31 / 0.2e1) * t31 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t22 / 0.2e1) * t22 + (-t3 * mrSges(5,1) - t5 * mrSges(5,2) + Ifges(5,3) * t20 / 0.2e1) * t20 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t9 + Ifges(6,1) * t7 / 0.2e1) * t7 + (t26 * mrSges(3,2) + t17 * mrSges(3,3) + Ifges(3,5) * t31 + Ifges(3,1) * t24 / 0.2e1) * t24 + (t17 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t22 + Ifges(4,1) * t16 / 0.2e1) * t16 + (t12 * mrSges(5,2) + t3 * mrSges(5,3) + Ifges(5,5) * t20 + Ifges(5,1) * t11 / 0.2e1) * t11 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t7 + Ifges(6,6) * t9 + Ifges(6,2) * t6 / 0.2e1) * t6 + (-t26 * mrSges(3,1) + t19 * mrSges(3,3) + Ifges(3,4) * t24 + Ifges(3,6) * t31 + Ifges(3,2) * t23 / 0.2e1) * t23 + (-t17 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t16 + Ifges(4,6) * t22 + Ifges(4,2) * t15 / 0.2e1) * t15 + (-t12 * mrSges(5,1) + t5 * mrSges(5,3) + Ifges(5,4) * t11 + Ifges(5,6) * t20 + Ifges(5,2) * t10 / 0.2e1) * t10 + (V_base(2) * mrSges(1,1) + t29 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t28 * mrSges(2,2) + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) - t32 * mrSges(2,1) + V_base(1) * mrSges(1,3) + t28 * mrSges(2,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(5) + (Ifges(1,6) + Ifges(2,6)) * V_base(6)) * V_base(5) + (V_base(3) * mrSges(1,2) + t32 * mrSges(2,2) - V_base(2) * mrSges(1,3) - t29 * mrSges(2,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(4) + (Ifges(1,5) + Ifges(2,5)) * V_base(6) + (Ifges(1,4) + Ifges(2,4)) * V_base(5)) * V_base(4);
T  = t4;
