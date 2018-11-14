% Calculate kinetic energy for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:02
% EndTime: 2018-11-14 13:38:03
% DurationCPUTime: 0.35s
% Computational Cost: add. (355->104), mult. (517->125), div. (0->0), fcn. (268->4), ass. (0->33)
t27 = sin(pkin(5));
t34 = cos(pkin(5));
t16 = t27 * V_base(4) - t34 * V_base(5);
t35 = pkin(1) + qJ(3);
t37 = t35 * t16;
t36 = pkin(2) + pkin(4);
t21 = V_base(5) * qJ(1) + V_base(1);
t22 = -V_base(4) * qJ(1) + V_base(2);
t15 = t34 * t21 + t27 * t22;
t13 = -V_base(6) * qJ(2) - t15;
t26 = V_base(3) + qJD(1);
t33 = qJD(3) - t13;
t14 = -t27 * t21 + t34 * t22;
t17 = t27 * V_base(5) + t34 * V_base(4);
t32 = t17 * qJ(2) - t26;
t31 = qJD(2) - t14;
t30 = -t35 * V_base(6) + t31;
t29 = cos(qJ(4));
t28 = sin(qJ(4));
t24 = V_base(6) + qJD(4);
t12 = -V_base(6) * pkin(1) + t31;
t11 = t29 * t16 + t28 * t17;
t10 = -t28 * t16 + t29 * t17;
t9 = t16 * pkin(1) - t32;
t8 = -t16 * pkin(2) + t33;
t7 = t17 * pkin(2) + t30;
t6 = t32 - t37;
t5 = V_base(6) * pkin(3) - t36 * t16 + t33;
t4 = t36 * t17 + t30;
t3 = (-pkin(3) - qJ(2)) * t17 + t37 + t26;
t2 = t28 * t5 + t29 * t4;
t1 = -t28 * t4 + t29 * t5;
t18 = m(2) * (t14 ^ 2 + t15 ^ 2 + t26 ^ 2) / 0.2e1 + m(3) * (t12 ^ 2 + t13 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(4) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,3) * t24 / 0.2e1) * t24 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t3 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,5) * t24 + Ifges(5,1) * t11 / 0.2e1) * t11 + (-t3 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t11 + Ifges(5,6) * t24 + Ifges(5,2) * t10 / 0.2e1) * t10 + (t12 * mrSges(3,1) + t6 * mrSges(4,1) + t26 * mrSges(2,2) - t7 * mrSges(4,2) - t14 * mrSges(2,3) - t9 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t17) * t17 + (t26 * mrSges(2,1) + t13 * mrSges(3,1) - t9 * mrSges(3,2) + t8 * mrSges(4,2) - t15 * mrSges(2,3) - t6 * mrSges(4,3) + (Ifges(2,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t16 + (-Ifges(2,4) + Ifges(4,5) - Ifges(3,6)) * t17) * t16 + (V_base(2) * mrSges(1,1) + t14 * mrSges(2,1) + t8 * mrSges(4,1) - V_base(1) * mrSges(1,2) - t15 * mrSges(2,2) + t12 * mrSges(3,2) - t13 * mrSges(3,3) - t7 * mrSges(4,3) + Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + (Ifges(2,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6) + (-Ifges(3,4) + Ifges(2,5) - Ifges(4,6)) * t17 + (-Ifges(4,4) + Ifges(3,5) - Ifges(2,6)) * t16) * V_base(6);
T  = t18;
