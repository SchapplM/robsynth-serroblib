% Calculate kinetic energy for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3PPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:01:49
% EndTime: 2019-03-08 18:01:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (157->80), mult. (198->92), div. (0->0), fcn. (32->2), ass. (0->20)
t21 = -pkin(2) - qJ(1);
t16 = V_base(2) + qJD(1);
t20 = V_base(4) * qJ(1);
t6 = V_base(4) * pkin(1) + V_base(6) * qJ(2) + t16;
t11 = V_base(6) * qJ(1) - V_base(1);
t19 = V_base(5) * qJ(2) - V_base(3);
t18 = cos(qJ(3));
t17 = sin(qJ(3));
t13 = -V_base(6) + qJD(3);
t12 = -V_base(3) - t20;
t10 = -t19 + t20;
t9 = t17 * V_base(4) + t18 * V_base(5);
t8 = -t17 * V_base(5) + t18 * V_base(4);
t7 = -V_base(5) * pkin(1) + qJD(2) - t11;
t5 = t21 * V_base(4) + t19;
t4 = V_base(4) * pkin(3) + t6;
t3 = qJD(2) + V_base(1) + t21 * V_base(6) + (-pkin(1) - pkin(3)) * V_base(5);
t2 = t17 * t3 + t18 * t4;
t1 = -t17 * t4 + t18 * t3;
t14 = m(3) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(2) * (t11 ^ 2 + t12 ^ 2 + t16 ^ 2) / 0.2e1 + m(4) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t5 * mrSges(4,2) - t1 * mrSges(4,3) + Ifges(4,1) * t9 / 0.2e1) * t9 + (-t5 * mrSges(4,1) + t2 * mrSges(4,3) + Ifges(4,4) * t9 + Ifges(4,2) * t8 / 0.2e1) * t8 + (t1 * mrSges(4,1) - t2 * mrSges(4,2) + Ifges(4,5) * t9 + Ifges(4,6) * t8 + Ifges(4,3) * t13 / 0.2e1) * t13 + (V_base(2) * mrSges(1,1) - t7 * mrSges(3,1) - V_base(1) * mrSges(1,2) - t16 * mrSges(2,2) + t11 * mrSges(2,3) + t6 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(1,3) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) + t11 * mrSges(2,1) - t12 * mrSges(2,2) + t7 * mrSges(3,2) + V_base(1) * mrSges(1,3) - t10 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5) + (Ifges(3,4) - Ifges(2,5) + Ifges(1,6)) * V_base(6)) * V_base(5) + (t16 * mrSges(2,1) + t10 * mrSges(3,1) + V_base(3) * mrSges(1,2) - t6 * mrSges(3,2) - V_base(2) * mrSges(1,3) - t12 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(1,1) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(4) + (Ifges(2,4) + Ifges(1,5) + Ifges(3,6)) * V_base(6) + (Ifges(1,4) + Ifges(3,5) - Ifges(2,6)) * V_base(5)) * V_base(4);
T  = t14;
