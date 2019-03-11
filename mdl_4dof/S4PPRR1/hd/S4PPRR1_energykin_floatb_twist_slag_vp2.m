% Calculate kinetic energy for
% S4PPRR1
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:37
% EndTime: 2019-03-08 18:15:37
% DurationCPUTime: 0.39s
% Computational Cost: add. (553->105), mult. (827->142), div. (0->0), fcn. (524->6), ass. (0->34)
t37 = sin(pkin(6));
t43 = cos(pkin(6));
t26 = t37 * V_base(5) + t43 * V_base(4);
t30 = V_base(5) * qJ(1) + V_base(1);
t31 = -V_base(4) * qJ(1) + V_base(2);
t21 = -t37 * t30 + t43 * t31;
t42 = qJD(2) - t21;
t13 = -t26 * pkin(4) + (-pkin(1) - pkin(2)) * V_base(6) + t42;
t22 = t43 * t30 + t37 * t31;
t20 = V_base(6) * qJ(2) + t22;
t25 = t37 * V_base(4) - t43 * V_base(5);
t15 = pkin(4) * t25 + t20;
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t6 = t39 * t13 + t41 * t15;
t36 = V_base(3) + qJD(1);
t5 = t41 * t13 - t15 * t39;
t34 = -V_base(6) + qJD(3);
t16 = t25 * pkin(1) - t26 * qJ(2) + t36;
t12 = -pkin(2) * t25 - t16;
t40 = cos(qJ(4));
t38 = sin(qJ(4));
t33 = qJD(4) + t34;
t19 = -V_base(6) * pkin(1) + t42;
t18 = t39 * t25 + t26 * t41;
t17 = t25 * t41 - t39 * t26;
t9 = t17 * t38 + t18 * t40;
t8 = t17 * t40 - t18 * t38;
t7 = -pkin(3) * t17 + t12;
t4 = pkin(5) * t17 + t6;
t3 = pkin(3) * t34 - pkin(5) * t18 + t5;
t2 = t3 * t38 + t4 * t40;
t1 = t3 * t40 - t38 * t4;
t10 = m(2) * (t21 ^ 2 + t22 ^ 2 + t36 ^ 2) / 0.2e1 + m(3) * (t16 ^ 2 + t19 ^ 2 + t20 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t7 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t9 / 0.2e1) * t9 + (t5 * mrSges(4,1) - t6 * mrSges(4,2) + Ifges(4,3) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,2) * t8 / 0.2e1) * t8 + (t12 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,5) * t34 + Ifges(4,1) * t18 / 0.2e1) * t18 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t8 + Ifges(5,3) * t33 / 0.2e1) * t33 + (-t12 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t18 + Ifges(4,6) * t34 + Ifges(4,2) * t17 / 0.2e1) * t17 + (t36 * mrSges(2,2) + t19 * mrSges(3,2) - t21 * mrSges(2,3) - t16 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t26) * t26 + (t36 * mrSges(2,1) + t16 * mrSges(3,1) - t20 * mrSges(3,2) - t22 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t25 + (-Ifges(2,4) + Ifges(3,5)) * t26) * t25 + (V_base(2) * mrSges(1,1) + t21 * mrSges(2,1) - t19 * mrSges(3,1) - V_base(1) * mrSges(1,2) - t22 * mrSges(2,2) + t20 * mrSges(3,3) + Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6) + (Ifges(3,4) + Ifges(2,5)) * t26 + (-Ifges(2,6) + Ifges(3,6)) * t25) * V_base(6);
T  = t10;
