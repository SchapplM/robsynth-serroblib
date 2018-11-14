% Calculate kinetic energy for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2018-11-14 13:42
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:58
% EndTime: 2018-11-14 13:40:58
% DurationCPUTime: 0.34s
% Computational Cost: add. (485->101), mult. (731->125), div. (0->0), fcn. (464->4), ass. (0->29)
t25 = V_base(5) * qJ(1) + V_base(1);
t26 = -V_base(4) * qJ(1) + V_base(2);
t29 = sin(pkin(5));
t30 = cos(pkin(5));
t16 = -t25 * t29 + t30 * t26;
t21 = t29 * V_base(5) + t30 * V_base(4);
t10 = V_base(6) * pkin(1) - pkin(4) * t21 + t16;
t17 = t30 * t25 + t29 * t26;
t20 = -t29 * V_base(4) + t30 * V_base(5);
t13 = pkin(4) * t20 + t17;
t31 = sin(qJ(2));
t35 = cos(qJ(2));
t8 = t31 * t10 + t35 * t13;
t34 = pkin(2) + qJ(4);
t28 = V_base(3) + qJD(1);
t27 = V_base(6) + qJD(2);
t5 = -qJ(3) * t27 - t8;
t7 = t35 * t10 - t31 * t13;
t18 = -pkin(1) * t20 + t28;
t33 = qJD(3) - t7;
t15 = t31 * t20 + t35 * t21;
t32 = -qJ(3) * t15 + t18;
t14 = -t35 * t20 + t21 * t31;
t6 = pkin(2) * t14 + t32;
t4 = -t27 * pkin(2) + t33;
t3 = t34 * t14 + t32;
t2 = -pkin(3) * t14 + qJD(4) - t5;
t1 = t15 * pkin(3) - t34 * t27 + t33;
t9 = m(2) * (t16 ^ 2 + t17 ^ 2 + t28 ^ 2) / 0.2e1 + m(3) * (t18 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t28 * mrSges(2,2) - t16 * mrSges(2,3) + Ifges(2,1) * t21 / 0.2e1) * t21 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t28 * mrSges(2,1) + t17 * mrSges(2,3) + Ifges(2,4) * t21 + Ifges(2,2) * t20 / 0.2e1) * t20 + (V_base(2) * mrSges(1,1) + t16 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t17 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t21 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t20 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t7 * mrSges(3,1) - t8 * mrSges(3,2) + t4 * mrSges(4,2) + t2 * mrSges(5,2) - t5 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t27) * t27 + (t4 * mrSges(4,1) + t1 * mrSges(5,1) + t18 * mrSges(3,2) - t3 * mrSges(5,2) - t7 * mrSges(3,3) - t6 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t15 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t27) * t15 + (t18 * mrSges(3,1) + t5 * mrSges(4,1) - t2 * mrSges(5,1) - t6 * mrSges(4,2) - t8 * mrSges(3,3) + t3 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t14 + (Ifges(5,4) + Ifges(4,5) - Ifges(3,6)) * t27 + (-Ifges(3,4) - Ifges(4,6) + Ifges(5,6)) * t15) * t14;
T  = t9;
