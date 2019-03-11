% Calculate kinetic energy for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_energykin_floatb_twist_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 17:59:54
% EndTime: 2019-03-08 17:59:54
% DurationCPUTime: 0.19s
% Computational Cost: add. (149->56), mult. (215->85), div. (0->0), fcn. (120->4), ass. (0->19)
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t9 = t16 * V_base(6) - t18 * V_base(4);
t11 = -t16 * V_base(1) - t18 * V_base(3);
t19 = V_base(2) ^ 2;
t17 = cos(qJ(2));
t15 = sin(qJ(2));
t14 = V_base(5) + qJD(1);
t12 = t16 * V_base(3) - t18 * V_base(1);
t10 = t12 ^ 2;
t8 = -t16 * V_base(4) - t18 * V_base(6);
t7 = qJD(2) + t9;
t6 = t8 * pkin(1) + V_base(2);
t5 = -t14 * pkin(1) + t11;
t4 = -t15 * t14 + t17 * t8;
t3 = -t17 * t14 - t15 * t8;
t2 = -t15 * t6 + t17 * t5;
t1 = -t15 * t5 - t17 * t6;
t13 = m(1) * (V_base(1) ^ 2 + V_base(3) ^ 2 + t19) / 0.2e1 + m(2) * (t11 ^ 2 + t10 + t19) / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t10) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(2) * mrSges(2,1) + t11 * mrSges(2,3) + Ifges(2,2) * t9 / 0.2e1) * t9 + (t1 * mrSges(3,1) - t2 * mrSges(3,2) + Ifges(3,3) * t7 / 0.2e1) * t7 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(2) * mrSges(2,2) - t12 * mrSges(2,3) + Ifges(2,4) * t9 + Ifges(2,1) * t8 / 0.2e1) * t8 + (t12 * mrSges(3,2) - t1 * mrSges(3,3) + Ifges(3,5) * t7 + Ifges(3,1) * t4 / 0.2e1) * t4 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t12 * mrSges(3,1) + t2 * mrSges(3,3) + Ifges(3,4) * t4 + Ifges(3,6) * t7 + Ifges(3,2) * t3 / 0.2e1) * t3 + (t12 * mrSges(2,1) - t11 * mrSges(2,2) + Ifges(2,5) * t8 + Ifges(2,6) * t9 + Ifges(2,3) * t14 / 0.2e1) * t14;
T  = t13;
