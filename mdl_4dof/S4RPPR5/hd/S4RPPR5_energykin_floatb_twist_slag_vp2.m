% Calculate kinetic energy for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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

function T = S4RPPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:41
% EndTime: 2019-12-31 16:39:41
% DurationCPUTime: 0.44s
% Computational Cost: add. (563->105), mult. (759->141), div. (0->0), fcn. (468->6), ass. (0->34)
t30 = V_base(5) * pkin(4) + V_base(1);
t31 = -V_base(4) * pkin(4) + V_base(2);
t38 = sin(qJ(1));
t42 = cos(qJ(1));
t22 = t42 * t30 + t38 * t31;
t34 = V_base(6) + qJD(1);
t20 = t34 * qJ(2) + t22;
t25 = t38 * V_base(4) - t42 * V_base(5);
t14 = qJ(3) * t25 + t20;
t35 = sin(pkin(6));
t36 = cos(pkin(6));
t26 = t38 * V_base(5) + t42 * V_base(4);
t21 = -t38 * t30 + t42 * t31;
t41 = qJD(2) - t21;
t9 = -t26 * qJ(3) + (-pkin(1) - pkin(2)) * t34 + t41;
t7 = t36 * t14 + t35 * t9;
t18 = t25 * pkin(1) - t26 * qJ(2) + V_base(3);
t6 = -t14 * t35 + t36 * t9;
t16 = t25 * t36 - t26 * t35;
t10 = -pkin(2) * t25 + qJD(3) - t18;
t40 = V_base(3) ^ 2;
t39 = cos(qJ(4));
t37 = sin(qJ(4));
t19 = -t34 * pkin(1) + t41;
t17 = t25 * t35 + t26 * t36;
t15 = qJD(4) - t16;
t13 = t17 * t39 - t34 * t37;
t12 = -t17 * t37 - t34 * t39;
t5 = -pkin(5) * t34 + t7;
t4 = pkin(3) * t34 - t6;
t3 = -pkin(3) * t16 - pkin(5) * t17 + t10;
t2 = t3 * t37 + t39 * t5;
t1 = t3 * t39 - t37 * t5;
t8 = m(2) * (t21 ^ 2 + t22 ^ 2 + t40) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t40) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t20 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t10 * mrSges(4,2) - t6 * mrSges(4,3) + Ifges(4,1) * t17 / 0.2e1) * t17 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,3) * t15 / 0.2e1) * t15 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t10 * mrSges(4,1) + t7 * mrSges(4,3) + Ifges(4,4) * t17 + Ifges(4,2) * t16 / 0.2e1) * t16 + (t4 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,5) * t15 + Ifges(5,1) * t13 / 0.2e1) * t13 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t4 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t13 + Ifges(5,6) * t15 + Ifges(5,2) * t12 / 0.2e1) * t12 + (V_base(3) * mrSges(2,2) + t19 * mrSges(3,2) - t21 * mrSges(2,3) - t18 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t26) * t26 + (V_base(3) * mrSges(2,1) + t18 * mrSges(3,1) - t20 * mrSges(3,2) - t22 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t25 + (-Ifges(2,4) + Ifges(3,5)) * t26) * t25 + (t21 * mrSges(2,1) - t19 * mrSges(3,1) - t6 * mrSges(4,1) - t22 * mrSges(2,2) + t7 * mrSges(4,2) + t20 * mrSges(3,3) - Ifges(4,5) * t17 - Ifges(4,6) * t16 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * t34 + (Ifges(3,4) + Ifges(2,5)) * t26 + (-Ifges(2,6) + Ifges(3,6)) * t25) * t34;
T = t8;
