% Calculate kinetic energy for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:47
% EndTime: 2019-12-31 16:38:47
% DurationCPUTime: 0.44s
% Computational Cost: add. (619->105), mult. (899->141), div. (0->0), fcn. (604->6), ass. (0->36)
t42 = pkin(2) + pkin(5);
t30 = V_base(5) * pkin(4) + V_base(1);
t31 = -V_base(4) * pkin(4) + V_base(2);
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t21 = -t30 * t35 + t37 * t31;
t26 = t35 * V_base(5) + t37 * V_base(4);
t32 = V_base(6) + qJD(1);
t14 = pkin(1) * t32 - qJ(2) * t26 + t21;
t22 = t37 * t30 + t35 * t31;
t25 = -t35 * V_base(4) + t37 * V_base(5);
t17 = qJ(2) * t25 + t22;
t33 = sin(pkin(6));
t41 = cos(pkin(6));
t10 = t33 * t14 + t41 * t17;
t7 = -qJ(3) * t32 - t10;
t9 = t41 * t14 - t33 * t17;
t23 = -pkin(1) * t25 + qJD(2) + V_base(3);
t40 = qJD(3) - t9;
t20 = t33 * t25 + t41 * t26;
t39 = -qJ(3) * t20 + t23;
t38 = V_base(3) ^ 2;
t36 = cos(qJ(4));
t34 = sin(qJ(4));
t19 = -t41 * t25 + t26 * t33;
t18 = qJD(4) + t20;
t13 = t19 * t34 + t32 * t36;
t12 = t19 * t36 - t32 * t34;
t8 = pkin(2) * t19 + t39;
t6 = -t32 * pkin(2) + t40;
t5 = t42 * t19 + t39;
t4 = -pkin(3) * t19 - t7;
t3 = t20 * pkin(3) - t42 * t32 + t40;
t2 = t3 * t34 + t36 * t5;
t1 = t3 * t36 - t34 * t5;
t11 = m(2) * (t21 ^ 2 + t22 ^ 2 + t38) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t38) / 0.2e1 + m(3) * (t10 ^ 2 + t23 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(4) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t21 * mrSges(2,3) + Ifges(2,1) * t26 / 0.2e1) * t26 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,3) * t18 / 0.2e1) * t18 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t22 * mrSges(2,3) + Ifges(2,4) * t26 + Ifges(2,2) * t25 / 0.2e1) * t25 + (t4 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,5) * t18 + Ifges(5,1) * t13 / 0.2e1) * t13 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t4 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t13 + Ifges(5,6) * t18 + Ifges(5,2) * t12 / 0.2e1) * t12 + (t6 * mrSges(4,1) + t23 * mrSges(3,2) - t9 * mrSges(3,3) - t8 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t20) * t20 + (t23 * mrSges(3,1) + t7 * mrSges(4,1) - t8 * mrSges(4,2) - t10 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t19 + (-Ifges(3,4) - Ifges(4,6)) * t20) * t19 + (t21 * mrSges(2,1) + t9 * mrSges(3,1) - t22 * mrSges(2,2) - t10 * mrSges(3,2) + t6 * mrSges(4,2) - t7 * mrSges(4,3) + Ifges(2,5) * t26 + Ifges(2,6) * t25 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t32 + (-Ifges(4,4) + Ifges(3,5)) * t20 + (Ifges(4,5) - Ifges(3,6)) * t19) * t32;
T = t11;
