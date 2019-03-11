% Calculate kinetic energy for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:58
% EndTime: 2019-03-08 18:23:58
% DurationCPUTime: 0.36s
% Computational Cost: add. (437->100), mult. (577->127), div. (0->0), fcn. (300->4), ass. (0->27)
t29 = V_base(2) + qJD(1);
t21 = V_base(6) * pkin(1) - V_base(4) * pkin(4) + t29;
t24 = -V_base(6) * qJ(1) + V_base(1);
t22 = V_base(5) * pkin(4) + t24;
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t14 = t31 * t21 + t33 * t22;
t19 = -t31 * V_base(4) + t33 * V_base(5);
t10 = pkin(5) * t19 + t14;
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t13 = t33 * t21 - t22 * t31;
t20 = t31 * V_base(5) + t33 * V_base(4);
t27 = V_base(6) + qJD(2);
t8 = pkin(2) * t27 - pkin(5) * t20 + t13;
t4 = t32 * t10 + t30 * t8;
t25 = -V_base(4) * qJ(1) - V_base(3);
t3 = -t10 * t30 + t32 * t8;
t23 = -V_base(5) * pkin(1) - t25;
t15 = -pkin(2) * t19 + t23;
t26 = qJD(3) + t27;
t12 = t19 * t30 + t20 * t32;
t11 = t19 * t32 - t20 * t30;
t5 = -pkin(3) * t11 + qJD(4) + t15;
t2 = qJ(4) * t11 + t4;
t1 = pkin(3) * t26 - qJ(4) * t12 + t3;
t6 = m(2) * (t24 ^ 2 + t25 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t23 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t13 * mrSges(3,1) - t14 * mrSges(3,2) + Ifges(3,3) * t27 / 0.2e1) * t27 + (t23 * mrSges(3,2) - t13 * mrSges(3,3) + Ifges(3,5) * t27 + Ifges(3,1) * t20 / 0.2e1) * t20 + (-t23 * mrSges(3,1) + t14 * mrSges(3,3) + Ifges(3,4) * t20 + Ifges(3,6) * t27 + Ifges(3,2) * t19 / 0.2e1) * t19 + (V_base(2) * mrSges(1,1) + t29 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t24 * mrSges(2,3) + (Ifges(1,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(6)) * V_base(6) + (t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t26) * t26 + (-V_base(3) * mrSges(1,1) + t25 * mrSges(2,1) - t24 * mrSges(2,2) + V_base(1) * mrSges(1,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5) + (Ifges(1,6) - Ifges(2,6)) * V_base(6)) * V_base(5) + (t15 * mrSges(4,2) + t5 * mrSges(5,2) - t3 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t12 + (Ifges(4,5) + Ifges(5,5)) * t26) * t12 + (V_base(3) * mrSges(1,2) + t29 * mrSges(2,2) - V_base(2) * mrSges(1,3) - t25 * mrSges(2,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(4) + (-Ifges(2,4) + Ifges(1,5)) * V_base(6) + (Ifges(1,4) + Ifges(2,5)) * V_base(5)) * V_base(4) + (-t15 * mrSges(4,1) - t5 * mrSges(5,1) + t4 * mrSges(4,3) + t2 * mrSges(5,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t11 + (Ifges(4,6) + Ifges(5,6)) * t26 + (Ifges(4,4) + Ifges(5,4)) * t12) * t11;
T  = t6;
