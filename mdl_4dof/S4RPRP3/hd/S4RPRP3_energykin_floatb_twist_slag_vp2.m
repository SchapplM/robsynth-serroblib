% Calculate kinetic energy for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:29
% EndTime: 2019-12-31 16:42:29
% DurationCPUTime: 0.46s
% Computational Cost: add. (711->104), mult. (1031->141), div. (0->0), fcn. (720->6), ass. (0->33)
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t27 = -t38 * V_base(4) + t40 * V_base(5);
t28 = t38 * V_base(5) + t40 * V_base(4);
t35 = sin(pkin(6));
t36 = cos(pkin(6));
t22 = t27 * t36 - t28 * t35;
t23 = t27 * t35 + t28 * t36;
t26 = -pkin(1) * t27 + qJD(2) + V_base(3);
t11 = -pkin(2) * t22 - pkin(5) * t23 + t26;
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t32 = V_base(5) * pkin(4) + V_base(1);
t33 = -V_base(4) * pkin(4) + V_base(2);
t24 = -t32 * t38 + t40 * t33;
t34 = V_base(6) + qJD(1);
t17 = pkin(1) * t34 - qJ(2) * t28 + t24;
t25 = t40 * t32 + t38 * t33;
t20 = qJ(2) * t27 + t25;
t13 = t35 * t17 + t36 * t20;
t8 = pkin(5) * t34 + t13;
t4 = t37 * t11 + t39 * t8;
t3 = t39 * t11 - t37 * t8;
t12 = t17 * t36 - t35 * t20;
t7 = -pkin(2) * t34 - t12;
t41 = V_base(3) ^ 2;
t21 = qJD(3) - t22;
t16 = t23 * t39 + t34 * t37;
t15 = -t23 * t37 + t34 * t39;
t5 = -pkin(3) * t15 + qJD(4) + t7;
t2 = qJ(4) * t15 + t4;
t1 = pkin(3) * t21 - qJ(4) * t16 + t3;
t6 = m(2) * (t24 ^ 2 + t25 ^ 2 + t41) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t41) / 0.2e1 + m(3) * (t12 ^ 2 + t13 ^ 2 + t26 ^ 2) / 0.2e1 + m(4) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t24 * mrSges(2,3) + Ifges(2,1) * t28 / 0.2e1) * t28 + (t26 * mrSges(3,2) - t12 * mrSges(3,3) + Ifges(3,1) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t25 * mrSges(2,3) + Ifges(2,4) * t28 + Ifges(2,2) * t27 / 0.2e1) * t27 + (-t26 * mrSges(3,1) + t13 * mrSges(3,3) + Ifges(3,4) * t23 + Ifges(3,2) * t22 / 0.2e1) * t22 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t21) * t21 + (t7 * mrSges(4,2) + t5 * mrSges(5,2) - t3 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t16 + (Ifges(4,5) + Ifges(5,5)) * t21) * t16 + (t24 * mrSges(2,1) + t12 * mrSges(3,1) - t25 * mrSges(2,2) - t13 * mrSges(3,2) + Ifges(2,5) * t28 + Ifges(3,5) * t23 + Ifges(2,6) * t27 + Ifges(3,6) * t22 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t34) * t34 + (-t7 * mrSges(4,1) - t5 * mrSges(5,1) + t4 * mrSges(4,3) + t2 * mrSges(5,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t15 + (Ifges(4,6) + Ifges(5,6)) * t21 + (Ifges(4,4) + Ifges(5,4)) * t16) * t15;
T = t6;
