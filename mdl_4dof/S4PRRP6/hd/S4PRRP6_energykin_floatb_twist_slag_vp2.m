% Calculate kinetic energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:10
% EndTime: 2019-12-31 16:30:10
% DurationCPUTime: 0.47s
% Computational Cost: add. (641->104), mult. (931->142), div. (0->0), fcn. (632->6), ass. (0->33)
t33 = V_base(5) * qJ(1) + V_base(1);
t34 = -V_base(4) * qJ(1) + V_base(2);
t36 = sin(pkin(6));
t37 = cos(pkin(6));
t24 = -t36 * t33 + t34 * t37;
t19 = -V_base(6) * pkin(1) - t24;
t29 = t36 * V_base(5) + t37 * V_base(4);
t39 = sin(qJ(2));
t40 = cos(qJ(2));
t22 = -t39 * t29 + t40 * V_base(6);
t23 = t29 * t40 + t39 * V_base(6);
t10 = -pkin(2) * t22 - pkin(5) * t23 + t19;
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t28 = -t36 * V_base(4) + t37 * V_base(5);
t35 = V_base(3) + qJD(1);
t16 = -pkin(1) * t28 - pkin(4) * t29 + t35;
t25 = t37 * t33 + t36 * t34;
t20 = V_base(6) * pkin(4) + t25;
t12 = t39 * t16 + t40 * t20;
t27 = qJD(2) - t28;
t8 = pkin(5) * t27 + t12;
t5 = t38 * t10 + t41 * t8;
t11 = t16 * t40 - t39 * t20;
t7 = -t27 * pkin(2) - t11;
t4 = t41 * t10 - t38 * t8;
t21 = qJD(3) - t22;
t14 = t41 * t23 + t38 * t27;
t13 = t23 * t38 - t41 * t27;
t3 = t13 * pkin(3) - t14 * qJ(4) + t7;
t2 = qJ(4) * t21 + t5;
t1 = -t21 * pkin(3) + qJD(4) - t4;
t6 = m(2) * (t24 ^ 2 + t25 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t11 ^ 2 + t12 ^ 2 + t19 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(4) * (t4 ^ 2 + t5 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t35 * mrSges(2,2) - t24 * mrSges(2,3) + Ifges(2,1) * t29 / 0.2e1) * t29 + (t11 * mrSges(3,1) - t12 * mrSges(3,2) + Ifges(3,3) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t35 * mrSges(2,1) + t25 * mrSges(2,3) + Ifges(2,4) * t29 + Ifges(2,2) * t28 / 0.2e1) * t28 + (t19 * mrSges(3,2) - t11 * mrSges(3,3) + Ifges(3,5) * t27 + Ifges(3,1) * t23 / 0.2e1) * t23 + (-t19 * mrSges(3,1) + t12 * mrSges(3,3) + Ifges(3,4) * t23 + Ifges(3,6) * t27 + Ifges(3,2) * t22 / 0.2e1) * t22 + (t4 * mrSges(4,1) - t1 * mrSges(5,1) - t5 * mrSges(4,2) + t2 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t21) * t21 + (t7 * mrSges(4,2) + t1 * mrSges(5,2) - t4 * mrSges(4,3) - t3 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t14 + (Ifges(5,4) + Ifges(4,5)) * t21) * t14 + (V_base(2) * mrSges(1,1) + t24 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t25 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t29 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t28 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t7 * mrSges(4,1) + t3 * mrSges(5,1) - t2 * mrSges(5,2) - t5 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t13 + (-Ifges(4,6) + Ifges(5,6)) * t21 + (-Ifges(4,4) + Ifges(5,5)) * t14) * t13;
T = t6;
