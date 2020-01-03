% Calculate kinetic energy for
% S4RPPR6
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:34
% EndTime: 2019-12-31 16:40:35
% DurationCPUTime: 0.55s
% Computational Cost: add. (609->105), mult. (821->141), div. (0->0), fcn. (536->6), ass. (0->36)
t46 = -pkin(2) - pkin(3);
t45 = cos(qJ(1));
t38 = sin(qJ(1));
t25 = t38 * V_base(4) - t45 * V_base(5);
t26 = t38 * V_base(5) + t45 * V_base(4);
t14 = pkin(1) * t25 - qJ(2) * t26 + V_base(3);
t32 = V_base(5) * pkin(4) + V_base(1);
t33 = -V_base(4) * pkin(4) + V_base(2);
t22 = t45 * t32 + t38 * t33;
t35 = V_base(6) + qJD(1);
t18 = qJ(2) * t35 + t22;
t36 = sin(pkin(6));
t44 = cos(pkin(6));
t10 = t36 * t14 + t44 * t18;
t21 = -t38 * t32 + t45 * t33;
t8 = t25 * qJ(3) + t10;
t9 = t44 * t14 - t36 * t18;
t43 = pkin(1) * t35 - qJD(2) + t21;
t42 = qJD(3) - t9;
t20 = t44 * t26 + t36 * t35;
t41 = qJ(3) * t20 + t43;
t40 = V_base(3) ^ 2;
t39 = cos(qJ(4));
t37 = sin(qJ(4));
t24 = qJD(4) - t25;
t19 = t26 * t36 - t44 * t35;
t12 = t19 * t37 + t20 * t39;
t11 = t19 * t39 - t20 * t37;
t7 = -t25 * pkin(2) + t42;
t6 = pkin(2) * t19 - t41;
t5 = pkin(5) * t19 + t8;
t4 = t46 * t19 + t41;
t3 = -t20 * pkin(5) + t46 * t25 + t42;
t2 = t3 * t37 + t39 * t5;
t1 = t3 * t39 - t37 * t5;
t13 = m(2) * (t21 ^ 2 + t22 ^ 2 + t40) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t40) / 0.2e1 + m(3) * (t10 ^ 2 + t43 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(4) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t21 * mrSges(2,1) - t22 * mrSges(2,2) + Ifges(2,3) * t35 / 0.2e1) * t35 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,3) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t21 * mrSges(2,3) + Ifges(2,5) * t35 + Ifges(2,1) * t26 / 0.2e1) * t26 + (t4 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,5) * t24 + Ifges(5,1) * t12 / 0.2e1) * t12 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t4 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t12 + Ifges(5,6) * t24 + Ifges(5,2) * t11 / 0.2e1) * t11 + (-t43 * mrSges(3,2) + t7 * mrSges(4,2) - t9 * mrSges(3,3) - t6 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t20) * t20 + (-t43 * mrSges(3,1) + t6 * mrSges(4,1) - t8 * mrSges(4,2) - t10 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t19 + (-Ifges(3,4) + Ifges(4,5)) * t20) * t19 + (V_base(3) * mrSges(2,1) + t9 * mrSges(3,1) - t7 * mrSges(4,1) - t10 * mrSges(3,2) - t22 * mrSges(2,3) + t8 * mrSges(4,3) - Ifges(2,4) * t26 - Ifges(2,6) * t35 + (Ifges(2,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t25 + (Ifges(4,4) + Ifges(3,5)) * t20 + (-Ifges(3,6) + Ifges(4,6)) * t19) * t25;
T = t13;
