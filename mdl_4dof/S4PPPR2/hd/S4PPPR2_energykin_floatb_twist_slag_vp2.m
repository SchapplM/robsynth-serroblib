% Calculate kinetic energy for
% S4PPPR2
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
%   pkin=[a2,a3,a4,d4,theta2]';
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
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:06
% EndTime: 2019-03-08 18:10:06
% DurationCPUTime: 0.38s
% Computational Cost: add. (337->101), mult. (455->125), div. (0->0), fcn. (196->4), ass. (0->29)
t35 = -pkin(2) - pkin(3);
t27 = V_base(2) + qJD(1);
t18 = V_base(6) * pkin(1) - V_base(4) * qJ(2) + t27;
t21 = -V_base(6) * qJ(1) + V_base(1);
t20 = V_base(5) * qJ(2) + t21;
t28 = sin(pkin(5));
t34 = cos(pkin(5));
t12 = t28 * t18 + t34 * t20;
t22 = -V_base(4) * qJ(1) - V_base(3);
t8 = V_base(6) * qJ(3) + t12;
t11 = t34 * t18 - t28 * t20;
t33 = qJD(3) - t11;
t32 = V_base(5) * pkin(1) - qJD(2) + t22;
t17 = t28 * V_base(5) + t34 * V_base(4);
t31 = qJ(3) * t17 + t32;
t30 = cos(qJ(4));
t29 = sin(qJ(4));
t24 = -V_base(6) + qJD(4);
t16 = t28 * V_base(4) - t34 * V_base(5);
t10 = t29 * t16 + t17 * t30;
t9 = t16 * t30 - t29 * t17;
t7 = -V_base(6) * pkin(2) + t33;
t6 = pkin(2) * t16 - t31;
t5 = pkin(4) * t16 + t8;
t4 = -t17 * pkin(4) + t35 * V_base(6) + t33;
t3 = t35 * t16 + t31;
t2 = t29 * t4 + t30 * t5;
t1 = -t29 * t5 + t30 * t4;
t13 = m(2) * (t21 ^ 2 + t22 ^ 2 + t27 ^ 2) / 0.2e1 + m(3) * (t11 ^ 2 + t12 ^ 2 + t32 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(4) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-t3 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,2) * t9 / 0.2e1) * t9 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,6) * t9 + Ifges(5,3) * t24 / 0.2e1) * t24 + (t3 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,5) * t24 + Ifges(5,1) * t10 / 0.2e1) * t10 + (-V_base(3) * mrSges(1,1) + t22 * mrSges(2,1) - t21 * mrSges(2,2) + V_base(1) * mrSges(1,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5)) * V_base(5) + (-t32 * mrSges(3,2) + t7 * mrSges(4,2) - t11 * mrSges(3,3) - t6 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t17) * t17 + (V_base(3) * mrSges(1,2) + t27 * mrSges(2,2) - V_base(2) * mrSges(1,3) - t22 * mrSges(2,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(4) + (Ifges(1,4) + Ifges(2,5)) * V_base(5)) * V_base(4) + (-t32 * mrSges(3,1) + t6 * mrSges(4,1) - t8 * mrSges(4,2) - t12 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t16 + (-Ifges(3,4) + Ifges(4,5)) * t17) * t16 + (V_base(2) * mrSges(1,1) + t27 * mrSges(2,1) + t11 * mrSges(3,1) - t7 * mrSges(4,1) - V_base(1) * mrSges(1,2) - t12 * mrSges(3,2) - t21 * mrSges(2,3) + t8 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(1,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(6) + (Ifges(1,6) - Ifges(2,6)) * V_base(5) + (-Ifges(2,4) + Ifges(1,5)) * V_base(4) + (Ifges(4,4) + Ifges(3,5)) * t17 + (-Ifges(3,6) + Ifges(4,6)) * t16) * V_base(6);
T  = t13;
