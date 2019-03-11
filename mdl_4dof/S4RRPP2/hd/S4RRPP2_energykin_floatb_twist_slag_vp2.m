% Calculate kinetic energy for
% S4RRPP2
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:58
% EndTime: 2019-03-08 18:33:58
% DurationCPUTime: 0.33s
% Computational Cost: add. (523->101), mult. (731->126), div. (0->0), fcn. (464->4), ass. (0->30)
t37 = -pkin(2) - pkin(3);
t26 = V_base(5) * pkin(4) + V_base(1);
t27 = -V_base(4) * pkin(4) + V_base(2);
t31 = sin(qJ(1));
t32 = cos(qJ(1));
t16 = -t26 * t31 + t32 * t27;
t21 = t31 * V_base(5) + t32 * V_base(4);
t29 = V_base(6) + qJD(1);
t10 = pkin(1) * t29 - pkin(5) * t21 + t16;
t17 = t32 * t26 + t31 * t27;
t20 = -t31 * V_base(4) + t32 * V_base(5);
t13 = pkin(5) * t20 + t17;
t30 = sin(qJ(2));
t36 = cos(qJ(2));
t8 = t30 * t10 + t36 * t13;
t28 = qJD(2) + t29;
t5 = t28 * qJ(3) + t8;
t18 = -pkin(1) * t20 + V_base(3);
t7 = t36 * t10 - t30 * t13;
t35 = qJD(3) - t7;
t15 = t30 * t20 + t36 * t21;
t34 = qJ(3) * t15 - t18;
t33 = V_base(3) ^ 2;
t14 = -t36 * t20 + t21 * t30;
t6 = pkin(2) * t14 - t34;
t4 = -t28 * pkin(2) + t35;
t3 = t37 * t14 + qJD(4) + t34;
t2 = qJ(4) * t14 + t5;
t1 = -t15 * qJ(4) + t37 * t28 + t35;
t9 = m(2) * (t16 ^ 2 + t17 ^ 2 + t33) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t33) / 0.2e1 + m(3) * (t18 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(4) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t16 * mrSges(2,1) - t17 * mrSges(2,2) + Ifges(2,3) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t16 * mrSges(2,3) + Ifges(2,5) * t29 + Ifges(2,1) * t21 / 0.2e1) * t21 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t17 * mrSges(2,3) + Ifges(2,4) * t21 + Ifges(2,6) * t29 + Ifges(2,2) * t20 / 0.2e1) * t20 + (t7 * mrSges(3,1) - t4 * mrSges(4,1) - t1 * mrSges(5,1) - t8 * mrSges(3,2) + t2 * mrSges(5,2) + t5 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t28) * t28 + (t18 * mrSges(3,2) + t4 * mrSges(4,2) + t3 * mrSges(5,2) - t7 * mrSges(3,3) - t6 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t15 + (Ifges(4,4) + Ifges(3,5) - Ifges(5,5)) * t28) * t15 + (t18 * mrSges(3,1) + t6 * mrSges(4,1) - t3 * mrSges(5,1) - t5 * mrSges(4,2) - t8 * mrSges(3,3) + t2 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t14 + (-Ifges(3,6) + Ifges(4,6) - Ifges(5,6)) * t28 + (-Ifges(3,4) + Ifges(5,4) + Ifges(4,5)) * t15) * t14;
T  = t9;
