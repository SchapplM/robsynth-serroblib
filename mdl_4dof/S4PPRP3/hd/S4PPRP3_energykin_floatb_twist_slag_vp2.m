% Calculate kinetic energy for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:19
% EndTime: 2019-03-08 18:14:19
% DurationCPUTime: 0.27s
% Computational Cost: add. (257->100), mult. (313->110), div. (0->0), fcn. (76->2), ass. (0->22)
t23 = V_base(2) + qJD(1);
t12 = V_base(4) * pkin(1) + V_base(6) * qJ(2) + t23;
t10 = V_base(4) * pkin(4) + t12;
t24 = sin(qJ(3));
t25 = cos(qJ(3));
t26 = -pkin(2) - qJ(1);
t8 = qJD(2) + V_base(1) + t26 * V_base(6) + (-pkin(1) - pkin(4)) * V_base(5);
t4 = t25 * t10 + t24 * t8;
t3 = -t24 * t10 + t25 * t8;
t18 = -V_base(4) * qJ(1) - V_base(3);
t17 = V_base(6) * qJ(1) - V_base(1);
t21 = V_base(5) * qJ(2);
t11 = t26 * V_base(4) + t21 - V_base(3);
t19 = -V_base(6) + qJD(3);
t16 = -t18 - t21;
t15 = t24 * V_base(4) + t25 * V_base(5);
t14 = -t24 * V_base(5) + t25 * V_base(4);
t13 = -V_base(5) * pkin(1) + qJD(2) - t17;
t5 = -t14 * pkin(3) + qJD(4) + t11;
t2 = t14 * qJ(4) + t4;
t1 = t19 * pkin(3) - t15 * qJ(4) + t3;
t6 = m(2) * (t17 ^ 2 + t18 ^ 2 + t23 ^ 2) / 0.2e1 + m(3) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t19) * t19 + (t11 * mrSges(4,2) + t5 * mrSges(5,2) - t3 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t15 + (Ifges(4,5) + Ifges(5,5)) * t19) * t15 + (-t11 * mrSges(4,1) - t5 * mrSges(5,1) + t4 * mrSges(4,3) + t2 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t14 + (Ifges(4,6) + Ifges(5,6)) * t19 + (Ifges(4,4) + Ifges(5,4)) * t15) * t14 + (V_base(2) * mrSges(1,1) - t13 * mrSges(3,1) - V_base(1) * mrSges(1,2) - t23 * mrSges(2,2) + t17 * mrSges(2,3) + t12 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(1,3) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) + t17 * mrSges(2,1) - t18 * mrSges(2,2) + t13 * mrSges(3,2) + V_base(1) * mrSges(1,3) - t16 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5) + (Ifges(3,4) - Ifges(2,5) + Ifges(1,6)) * V_base(6)) * V_base(5) + (t23 * mrSges(2,1) + t16 * mrSges(3,1) + V_base(3) * mrSges(1,2) - t12 * mrSges(3,2) - V_base(2) * mrSges(1,3) - t18 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(1,1) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(4) + (Ifges(2,4) + Ifges(1,5) + Ifges(3,6)) * V_base(6) + (Ifges(1,4) + Ifges(3,5) - Ifges(2,6)) * V_base(5)) * V_base(4);
T  = t6;
