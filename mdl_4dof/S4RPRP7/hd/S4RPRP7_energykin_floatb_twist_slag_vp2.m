% Calculate kinetic energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:00
% EndTime: 2019-12-31 16:47:00
% DurationCPUTime: 0.39s
% Computational Cost: add. (457->101), mult. (585->126), div. (0->0), fcn. (328->4), ass. (0->30)
t37 = pkin(1) + pkin(5);
t30 = sin(qJ(3));
t35 = cos(qJ(3));
t31 = sin(qJ(1));
t36 = cos(qJ(1));
t21 = t31 * V_base(5) + t36 * V_base(4);
t29 = V_base(6) + qJD(1);
t25 = V_base(5) * pkin(4) + V_base(1);
t26 = -V_base(4) * pkin(4) + V_base(2);
t16 = -t31 * t25 + t36 * t26;
t33 = qJD(2) - t16;
t7 = t21 * pkin(2) - t37 * t29 + t33;
t20 = t31 * V_base(4) - t36 * V_base(5);
t34 = -qJ(2) * t21 + V_base(3);
t9 = t37 * t20 + t34;
t4 = t30 * t7 + t35 * t9;
t17 = t36 * t25 + t31 * t26;
t13 = -t29 * qJ(2) - t17;
t10 = -pkin(2) * t20 - t13;
t3 = -t30 * t9 + t35 * t7;
t32 = V_base(3) ^ 2;
t19 = qJD(3) + t21;
t15 = t30 * t20 + t35 * t29;
t14 = -t35 * t20 + t29 * t30;
t12 = -t29 * pkin(1) + t33;
t11 = pkin(1) * t20 + t34;
t5 = pkin(3) * t14 - qJ(4) * t15 + t10;
t2 = qJ(4) * t19 + t4;
t1 = -t19 * pkin(3) + qJD(4) - t3;
t6 = m(2) * (t16 ^ 2 + t17 ^ 2 + t32) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t32) / 0.2e1 + m(4) * (t10 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(3) * (t11 ^ 2 + t12 ^ 2 + t13 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t16 * mrSges(2,1) - t17 * mrSges(2,2) + t12 * mrSges(3,2) - t13 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t29) * t29 + (t3 * mrSges(4,1) - t1 * mrSges(5,1) - t4 * mrSges(4,2) + t2 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t19) * t19 + (t12 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t16 * mrSges(2,3) - t11 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t21 + (-Ifges(3,4) + Ifges(2,5)) * t29) * t21 + (t10 * mrSges(4,2) + t1 * mrSges(5,2) - t3 * mrSges(4,3) - t5 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t15 + (Ifges(5,4) + Ifges(4,5)) * t19) * t15 + (V_base(3) * mrSges(2,1) + t13 * mrSges(3,1) - t11 * mrSges(3,2) - t17 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t20 + (Ifges(3,5) - Ifges(2,6)) * t29 + (-Ifges(2,4) - Ifges(3,6)) * t21) * t20 + (t10 * mrSges(4,1) + t5 * mrSges(5,1) - t2 * mrSges(5,2) - t4 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t14 + (-Ifges(4,6) + Ifges(5,6)) * t19 + (-Ifges(4,4) + Ifges(5,5)) * t15) * t14;
T = t6;
