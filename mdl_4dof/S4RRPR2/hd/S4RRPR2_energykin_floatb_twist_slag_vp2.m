% Calculate kinetic energy for
% S4RRPR2
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:28
% EndTime: 2019-07-18 18:16:29
% DurationCPUTime: 0.47s
% Computational Cost: add. (635->103), mult. (881->141), div. (0->0), fcn. (584->6), ass. (0->36)
t42 = -pkin(2) - pkin(3);
t41 = cos(qJ(2));
t29 = V_base(5) * pkin(4) + V_base(1);
t30 = -V_base(4) * pkin(4) + V_base(2);
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t19 = -t29 * t35 + t37 * t30;
t24 = t35 * V_base(5) + t37 * V_base(4);
t32 = V_base(6) + qJD(1);
t13 = pkin(1) * t32 - pkin(5) * t24 + t19;
t20 = t37 * t29 + t35 * t30;
t23 = -t35 * V_base(4) + t37 * V_base(5);
t16 = pkin(5) * t23 + t20;
t34 = sin(qJ(2));
t9 = t34 * t13 + t41 * t16;
t21 = -pkin(1) * t23 + V_base(3);
t8 = t41 * t13 - t34 * t16;
t31 = qJD(2) + t32;
t40 = qJD(3) - t8;
t18 = t34 * t23 + t41 * t24;
t39 = qJ(3) * t18 - t21;
t38 = V_base(3) ^ 2;
t36 = cos(qJ(4));
t33 = sin(qJ(4));
t28 = qJD(4) - t31;
t17 = -t41 * t23 + t24 * t34;
t11 = t17 * t33 + t18 * t36;
t10 = t17 * t36 - t18 * t33;
t7 = pkin(2) * t17 - t39;
t6 = qJ(3) * t31 + t9;
t5 = -t31 * pkin(2) + t40;
t4 = t42 * t31 + t40;
t3 = t42 * t17 + t39;
t2 = t33 * t4 + t36 * t6;
t1 = -t33 * t6 + t36 * t4;
t12 = m(2) * (t19 ^ 2 + t20 ^ 2 + t38) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t38) / 0.2e1 + m(3) * (t21 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(4) * (t5 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t19 * mrSges(2,1) - t20 * mrSges(2,2) + Ifges(2,3) * t32 / 0.2e1) * t32 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t19 * mrSges(2,3) + Ifges(2,5) * t32 + Ifges(2,1) * t24 / 0.2e1) * t24 + (t3 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t11 / 0.2e1) * t11 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t20 * mrSges(2,3) + Ifges(2,4) * t24 + Ifges(2,6) * t32 + Ifges(2,2) * t23 / 0.2e1) * t23 + (-t3 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t11 + Ifges(5,6) * t28 + Ifges(5,2) * t10 / 0.2e1) * t10 + (t8 * mrSges(3,1) - t5 * mrSges(4,1) - t9 * mrSges(3,2) + t6 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t31) * t31 + (t21 * mrSges(3,2) + t5 * mrSges(4,2) - t8 * mrSges(3,3) - t7 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t18 + (Ifges(4,4) + Ifges(3,5)) * t31) * t18 + (t21 * mrSges(3,1) + t7 * mrSges(4,1) - t6 * mrSges(4,2) - t9 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t17 + (-Ifges(3,6) + Ifges(4,6)) * t31 + (-Ifges(3,4) + Ifges(4,5)) * t18) * t17;
T  = t12;
