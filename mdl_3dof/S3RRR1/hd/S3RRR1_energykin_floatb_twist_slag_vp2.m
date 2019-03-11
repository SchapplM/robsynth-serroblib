% Calculate kinetic energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3RRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:46
% EndTime: 2019-03-08 18:07:47
% DurationCPUTime: 0.36s
% Computational Cost: add. (471->84), mult. (694->125), div. (0->0), fcn. (468->6), ass. (0->31)
t26 = V_base(5) * pkin(3) + V_base(1);
t27 = -V_base(4) * pkin(3) + V_base(2);
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t17 = -t26 * t32 + t35 * t27;
t21 = t32 * V_base(5) + t35 * V_base(4);
t29 = V_base(6) + qJD(1);
t12 = pkin(1) * t29 - pkin(4) * t21 + t17;
t18 = t35 * t26 + t32 * t27;
t20 = -t32 * V_base(4) + t35 * V_base(5);
t14 = pkin(4) * t20 + t18;
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t6 = t31 * t12 + t34 * t14;
t5 = t34 * t12 - t14 * t31;
t19 = -pkin(1) * t20 + V_base(3);
t28 = qJD(2) + t29;
t36 = V_base(3) ^ 2;
t33 = cos(qJ(3));
t30 = sin(qJ(3));
t25 = qJD(3) + t28;
t16 = t20 * t31 + t21 * t34;
t15 = t20 * t34 - t21 * t31;
t9 = -pkin(2) * t15 + t19;
t8 = t15 * t30 + t16 * t33;
t7 = t15 * t33 - t16 * t30;
t4 = pkin(5) * t15 + t6;
t3 = pkin(2) * t28 - pkin(5) * t16 + t5;
t2 = t3 * t30 + t33 * t4;
t1 = t3 * t33 - t30 * t4;
t10 = m(2) * (t17 ^ 2 + t18 ^ 2 + t36) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t36) / 0.2e1 + m(3) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t9 * mrSges(4,2) - t1 * mrSges(4,3) + Ifges(4,1) * t8 / 0.2e1) * t8 + (t17 * mrSges(2,1) - t18 * mrSges(2,2) + Ifges(2,3) * t29 / 0.2e1) * t29 + (t5 * mrSges(3,1) - t6 * mrSges(3,2) + Ifges(3,3) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t9 * mrSges(4,1) + t2 * mrSges(4,3) + Ifges(4,4) * t8 + Ifges(4,2) * t7 / 0.2e1) * t7 + (V_base(3) * mrSges(2,2) - t17 * mrSges(2,3) + Ifges(2,5) * t29 + Ifges(2,1) * t21 / 0.2e1) * t21 + (t19 * mrSges(3,2) - t5 * mrSges(3,3) + Ifges(3,5) * t28 + Ifges(3,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(4,1) - t2 * mrSges(4,2) + Ifges(4,5) * t8 + Ifges(4,6) * t7 + Ifges(4,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(2,1) + t18 * mrSges(2,3) + Ifges(2,4) * t21 + Ifges(2,6) * t29 + Ifges(2,2) * t20 / 0.2e1) * t20 + (-t19 * mrSges(3,1) + t6 * mrSges(3,3) + Ifges(3,4) * t16 + Ifges(3,6) * t28 + Ifges(3,2) * t15 / 0.2e1) * t15;
T  = t10;
