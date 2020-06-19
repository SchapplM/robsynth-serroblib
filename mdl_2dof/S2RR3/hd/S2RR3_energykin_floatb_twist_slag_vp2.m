% Calculate kinetic energy for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
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
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_energykin_floatb_twist_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:21
% EndTime: 2020-06-19 09:14:22
% DurationCPUTime: 0.42s
% Computational Cost: add. (187->60), mult. (273->90), div. (0->0), fcn. (148->4), ass. (0->21)
t15 = V_base(5) * pkin(2) + V_base(1);
t16 = -V_base(4) * pkin(2) + V_base(2);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t8 = t22 * t15 + t20 * t16;
t7 = -t20 * t15 + t22 * t16;
t18 = V_base(6) + qJD(1);
t23 = V_base(3) ^ 2;
t21 = cos(qJ(2));
t19 = sin(qJ(2));
t17 = qJD(2) + t18;
t11 = t20 * V_base(5) + t22 * V_base(4);
t10 = -t20 * V_base(4) + t22 * V_base(5);
t9 = -t10 * pkin(1) + V_base(3);
t6 = t19 * t10 + t21 * t11;
t5 = t21 * t10 - t19 * t11;
t4 = t10 * pkin(3) + t8;
t3 = t18 * pkin(1) - t11 * pkin(3) + t7;
t2 = t19 * t3 + t21 * t4;
t1 = -t19 * t4 + t21 * t3;
t12 = m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t23) / 0.2e1 + m(2) * (t7 ^ 2 + t8 ^ 2 + t23) / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t9 * mrSges(3,2) - t1 * mrSges(3,3) + Ifges(3,1) * t6 / 0.2e1) * t6 + (t7 * mrSges(2,1) - t8 * mrSges(2,2) + Ifges(2,3) * t18 / 0.2e1) * t18 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t9 * mrSges(3,1) + t2 * mrSges(3,3) + Ifges(3,4) * t6 + Ifges(3,2) * t5 / 0.2e1) * t5 + (V_base(3) * mrSges(2,2) - t7 * mrSges(2,3) + Ifges(2,5) * t18 + Ifges(2,1) * t11 / 0.2e1) * t11 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(3,1) - t2 * mrSges(3,2) + Ifges(3,5) * t6 + Ifges(3,6) * t5 + Ifges(3,3) * t17 / 0.2e1) * t17 + (-V_base(3) * mrSges(2,1) + t8 * mrSges(2,3) + Ifges(2,4) * t11 + Ifges(2,6) * t18 + Ifges(2,2) * t10 / 0.2e1) * t10;
T = t12;
