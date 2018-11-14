% Calculate kinetic energy for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:56:22
% EndTime: 2018-11-14 13:56:22
% DurationCPUTime: 0.32s
% Computational Cost: add. (329->104), mult. (419->125), div. (0->0), fcn. (148->4), ass. (0->28)
t32 = -pkin(2) - qJ(1);
t12 = qJD(2) + V_base(1) + t32 * V_base(6) + (-pkin(1) - qJ(3)) * V_base(5);
t27 = V_base(2) + qJD(1);
t18 = V_base(4) * pkin(1) + V_base(6) * qJ(2) + t27;
t15 = V_base(4) * qJ(3) + t18;
t28 = sin(pkin(5));
t29 = cos(pkin(5));
t6 = t28 * t12 + t29 * t15;
t5 = t29 * t12 - t28 * t15;
t22 = -V_base(4) * qJ(1) - V_base(3);
t21 = V_base(6) * qJ(1) - V_base(1);
t25 = V_base(5) * qJ(2);
t14 = t32 * V_base(4) + qJD(3) + t25 - V_base(3);
t31 = cos(qJ(4));
t30 = sin(qJ(4));
t23 = -V_base(6) + qJD(4);
t20 = -t22 - t25;
t19 = -V_base(5) * pkin(1) + qJD(2) - t21;
t17 = t28 * V_base(4) + t29 * V_base(5);
t16 = -t28 * V_base(5) + t29 * V_base(4);
t9 = -t16 * pkin(3) + t14;
t8 = t30 * t16 + t31 * t17;
t7 = t31 * t16 - t30 * t17;
t4 = t16 * pkin(4) + t6;
t3 = -V_base(6) * pkin(3) - t17 * pkin(4) + t5;
t2 = t30 * t3 + t31 * t4;
t1 = t31 * t3 - t30 * t4;
t10 = m(2) * (t21 ^ 2 + t22 ^ 2 + t27 ^ 2) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t20 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t9 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t8 / 0.2e1) * t8 + (t14 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,1) * t17 / 0.2e1) * t17 + (-t9 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t8 + Ifges(5,2) * t7 / 0.2e1) * t7 + (-t14 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t17 + Ifges(4,2) * t16 / 0.2e1) * t16 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t8 + Ifges(5,6) * t7 + Ifges(5,3) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + t21 * mrSges(2,1) - t22 * mrSges(2,2) + t19 * mrSges(3,2) + V_base(1) * mrSges(1,3) - t20 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5)) * V_base(5) + (t27 * mrSges(2,1) + t20 * mrSges(3,1) + V_base(3) * mrSges(1,2) - t18 * mrSges(3,2) - V_base(2) * mrSges(1,3) - t22 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(1,1) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(4) + (Ifges(1,4) + Ifges(3,5) - Ifges(2,6)) * V_base(5)) * V_base(4) + (V_base(2) * mrSges(1,1) - t19 * mrSges(3,1) - t5 * mrSges(4,1) - V_base(1) * mrSges(1,2) - t27 * mrSges(2,2) + t6 * mrSges(4,2) + t21 * mrSges(2,3) + t18 * mrSges(3,3) - Ifges(4,5) * t17 - Ifges(4,6) * t16 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(1,3) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(6) + (Ifges(3,4) - Ifges(2,5) + Ifges(1,6)) * V_base(5) + (Ifges(2,4) + Ifges(1,5) + Ifges(3,6)) * V_base(4)) * V_base(6);
T  = t10;
