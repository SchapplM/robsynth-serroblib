% Calculate kinetic energy for
% S4PRRP4
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:43
% EndTime: 2019-12-31 16:27:43
% DurationCPUTime: 0.46s
% Computational Cost: add. (689->104), mult. (1027->142), div. (0->0), fcn. (716->6), ass. (0->33)
t36 = sin(pkin(6));
t37 = cos(pkin(6));
t26 = -t36 * V_base(4) + t37 * V_base(5);
t27 = t36 * V_base(5) + t37 * V_base(4);
t39 = sin(qJ(2));
t40 = cos(qJ(2));
t21 = t26 * t40 - t39 * t27;
t22 = t39 * t26 + t27 * t40;
t35 = V_base(3) + qJD(1);
t25 = -pkin(1) * t26 + t35;
t10 = -pkin(2) * t21 - pkin(5) * t22 + t25;
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t32 = V_base(5) * qJ(1) + V_base(1);
t33 = -V_base(4) * qJ(1) + V_base(2);
t23 = -t32 * t36 + t37 * t33;
t16 = V_base(6) * pkin(1) - pkin(4) * t27 + t23;
t24 = t37 * t32 + t36 * t33;
t19 = pkin(4) * t26 + t24;
t12 = t39 * t16 + t40 * t19;
t34 = V_base(6) + qJD(2);
t9 = pkin(5) * t34 + t12;
t4 = t38 * t10 + t41 * t9;
t11 = t16 * t40 - t39 * t19;
t8 = -t34 * pkin(2) - t11;
t3 = t41 * t10 - t38 * t9;
t20 = qJD(3) - t21;
t14 = t41 * t22 + t38 * t34;
t13 = t22 * t38 - t41 * t34;
t5 = t13 * pkin(3) - t14 * qJ(4) + t8;
t2 = qJ(4) * t20 + t4;
t1 = -t20 * pkin(3) + qJD(4) - t3;
t6 = m(2) * (t23 ^ 2 + t24 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(4) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t11 * mrSges(3,1) - t12 * mrSges(3,2) + Ifges(3,3) * t34 / 0.2e1) * t34 + (t35 * mrSges(2,2) - t23 * mrSges(2,3) + Ifges(2,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t35 * mrSges(2,1) + t24 * mrSges(2,3) + Ifges(2,4) * t27 + Ifges(2,2) * t26 / 0.2e1) * t26 + (t25 * mrSges(3,2) - t11 * mrSges(3,3) + Ifges(3,5) * t34 + Ifges(3,1) * t22 / 0.2e1) * t22 + (-t25 * mrSges(3,1) + t12 * mrSges(3,3) + Ifges(3,4) * t22 + Ifges(3,6) * t34 + Ifges(3,2) * t21 / 0.2e1) * t21 + (t3 * mrSges(4,1) - t1 * mrSges(5,1) - t4 * mrSges(4,2) + t2 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t20) * t20 + (t8 * mrSges(4,2) + t1 * mrSges(5,2) - t3 * mrSges(4,3) - t5 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t14 + (Ifges(5,4) + Ifges(4,5)) * t20) * t14 + (V_base(2) * mrSges(1,1) + t23 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t24 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t27 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t26 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t8 * mrSges(4,1) + t5 * mrSges(5,1) - t2 * mrSges(5,2) - t4 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t13 + (-Ifges(4,6) + Ifges(5,6)) * t20 + (-Ifges(4,4) + Ifges(5,5)) * t14) * t13;
T = t6;
