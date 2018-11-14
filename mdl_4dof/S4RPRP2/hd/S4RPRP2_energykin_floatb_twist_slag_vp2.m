% Calculate kinetic energy for
% S4RPRP2
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
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:21
% EndTime: 2018-11-14 13:49:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (463->101), mult. (617->126), div. (0->0), fcn. (356->4), ass. (0->28)
t26 = V_base(5) * pkin(4) + V_base(1);
t27 = -V_base(4) * pkin(4) + V_base(2);
t33 = sin(qJ(1));
t37 = cos(qJ(1));
t18 = t37 * t26 + t33 * t27;
t31 = V_base(6) + qJD(1);
t16 = t31 * qJ(2) + t18;
t21 = t33 * V_base(4) - t37 * V_base(5);
t11 = pkin(5) * t21 + t16;
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t22 = t33 * V_base(5) + t37 * V_base(4);
t17 = -t33 * t26 + t37 * t27;
t36 = qJD(2) - t17;
t8 = -t22 * pkin(5) + (-pkin(1) - pkin(2)) * t31 + t36;
t4 = t34 * t11 + t32 * t8;
t3 = -t11 * t32 + t34 * t8;
t14 = t21 * pkin(1) - t22 * qJ(2) + V_base(3);
t9 = -pkin(2) * t21 - t14;
t35 = V_base(3) ^ 2;
t29 = qJD(3) - t31;
t15 = -t31 * pkin(1) + t36;
t13 = t21 * t32 + t22 * t34;
t12 = t21 * t34 - t22 * t32;
t5 = -pkin(3) * t12 + qJD(4) + t9;
t2 = qJ(4) * t12 + t4;
t1 = pkin(3) * t29 - qJ(4) * t13 + t3;
t6 = m(2) * (t17 ^ 2 + t18 ^ 2 + t35) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t35) / 0.2e1 + m(3) * (t14 ^ 2 + t15 ^ 2 + t16 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(4) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t17 * mrSges(2,1) - t15 * mrSges(3,1) - t18 * mrSges(2,2) + t16 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * t31) * t31 + (t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29) * t29 + (V_base(3) * mrSges(2,2) + t15 * mrSges(3,2) - t17 * mrSges(2,3) - t14 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t22 + (Ifges(3,4) + Ifges(2,5)) * t31) * t22 + (t9 * mrSges(4,2) + t5 * mrSges(5,2) - t3 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t13 + (Ifges(4,5) + Ifges(5,5)) * t29) * t13 + (V_base(3) * mrSges(2,1) + t14 * mrSges(3,1) - t16 * mrSges(3,2) - t18 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t21 + (-Ifges(2,6) + Ifges(3,6)) * t31 + (-Ifges(2,4) + Ifges(3,5)) * t22) * t21 + (-t9 * mrSges(4,1) - t5 * mrSges(5,1) + t4 * mrSges(4,3) + t2 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t12 + (Ifges(4,6) + Ifges(5,6)) * t29 + (Ifges(4,4) + Ifges(5,4)) * t13) * t12;
T  = t6;
