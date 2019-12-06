% Calculate kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:33
% EndTime: 2019-12-05 18:08:34
% DurationCPUTime: 0.71s
% Computational Cost: add. (637->112), mult. (854->159), div. (0->0), fcn. (636->8), ass. (0->36)
t35 = sin(qJ(1));
t40 = cos(qJ(1));
t27 = t35 * V_base(2) + t40 * V_base(1);
t31 = V_base(6) + qJD(1);
t18 = qJ(2) * t31 + t27;
t25 = t35 * V_base(5) + t40 * V_base(4);
t19 = -qJ(2) * t25 + V_base(3);
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t13 = t18 * t38 + t19 * t34;
t26 = -t35 * V_base(1) + t40 * V_base(2);
t23 = qJD(2) - t26;
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t5 = t13 * t33 - t37 * t23;
t42 = t5 ^ 2;
t11 = t18 * t34 - t38 * t19;
t41 = t11 ^ 2;
t24 = t35 * V_base(4) - t40 * V_base(5);
t16 = t25 * t38 + t31 * t34;
t21 = qJD(3) + t24;
t9 = -t16 * t33 + t21 * t37;
t15 = -t25 * t34 + t31 * t38;
t39 = V_base(3) ^ 2;
t36 = cos(qJ(5));
t32 = sin(qJ(5));
t22 = t23 ^ 2;
t14 = qJD(4) - t15;
t10 = t16 * t37 + t21 * t33;
t8 = qJD(5) - t9;
t7 = t13 * t37 + t23 * t33;
t4 = t10 * t36 + t14 * t32;
t3 = -t10 * t32 + t14 * t36;
t2 = t11 * t32 + t36 * t7;
t1 = t11 * t36 - t32 * t7;
t6 = m(2) * (t26 ^ 2 + t27 ^ 2 + t39) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t39) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t22) / 0.2e1 + m(4) * (t13 ^ 2 + t22 + t41) / 0.2e1 + m(5) * (t7 ^ 2 + t41 + t42) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t42) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t11 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,2) * t9 / 0.2e1) * t9 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t8 / 0.2e1) * t8 + (-t11 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,3) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t5 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t8 + Ifges(6,1) * t4 / 0.2e1) * t4 + (t23 * mrSges(4,2) + t11 * mrSges(4,3) + Ifges(4,5) * t21 + Ifges(4,1) * t16 / 0.2e1) * t16 + (-t5 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,6) * t9 + Ifges(5,3) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t5 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t4 + Ifges(6,6) * t8 + Ifges(6,2) * t3 / 0.2e1) * t3 + (-t23 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t16 + Ifges(4,6) * t21 + Ifges(4,2) * t15 / 0.2e1) * t15 + (t11 * mrSges(5,2) + t5 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,5) * t14 + Ifges(5,1) * t10 / 0.2e1) * t10 + (t26 * mrSges(2,1) - t23 * mrSges(3,1) - t27 * mrSges(2,2) + t18 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * t31) * t31 + (V_base(3) * mrSges(2,2) + t23 * mrSges(3,2) - t26 * mrSges(2,3) - t19 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t25 + (Ifges(3,4) + Ifges(2,5)) * t31) * t25 + (V_base(3) * mrSges(2,1) + t19 * mrSges(3,1) - t18 * mrSges(3,2) - t27 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t24 + (-Ifges(2,6) + Ifges(3,6)) * t31 + (-Ifges(2,4) + Ifges(3,5)) * t25) * t24;
T = t6;
