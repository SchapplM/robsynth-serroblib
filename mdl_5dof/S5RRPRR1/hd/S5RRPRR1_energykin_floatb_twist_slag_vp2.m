% Calculate kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:37
% EndTime: 2019-07-18 17:20:38
% DurationCPUTime: 0.72s
% Computational Cost: add. (941->122), mult. (1258->168), div. (0->0), fcn. (996->8), ass. (0->42)
t44 = sin(qJ(1));
t48 = cos(qJ(1));
t31 = t44 * V_base(5) + t48 * V_base(4);
t40 = V_base(6) + qJD(1);
t43 = sin(qJ(2));
t47 = cos(qJ(2));
t22 = -t31 * t43 + t40 * t47;
t34 = t44 * V_base(2) + t48 * V_base(1);
t25 = t47 * t34 + t43 * V_base(3);
t19 = t22 * qJ(3) + t25;
t13 = pkin(3) * t22 + t19;
t42 = sin(qJ(4));
t46 = cos(qJ(4));
t23 = t31 * t47 + t40 * t43;
t30 = -t44 * V_base(4) + t48 * V_base(5);
t29 = qJD(2) - t30;
t24 = -t34 * t43 + t47 * V_base(3);
t51 = t29 * pkin(1) + t24;
t50 = pkin(2) * t29 + (-pkin(3) - qJ(3)) * t23 + t51;
t4 = t13 * t42 - t46 * t50;
t54 = t4 ^ 2;
t32 = t44 * V_base(1) - t48 * V_base(2);
t53 = t32 ^ 2;
t6 = t46 * t13 + t42 * t50;
t52 = qJD(3) + t32;
t17 = t22 * t46 - t23 * t42;
t14 = (-pkin(1) - pkin(2)) * t22 + t52;
t49 = V_base(3) ^ 2;
t45 = cos(qJ(5));
t41 = sin(qJ(5));
t27 = qJD(4) + t29;
t20 = -pkin(1) * t22 + t52;
t18 = t22 * t42 + t23 * t46;
t16 = qJD(5) - t17;
t15 = -qJ(3) * t23 + t51;
t11 = t18 * t45 + t27 * t41;
t10 = -t18 * t41 + t27 * t45;
t7 = -pkin(4) * t18 + t14;
t3 = pkin(4) * t27 + t6;
t2 = t3 * t45 + t41 * t7;
t1 = -t3 * t41 + t45 * t7;
t5 = m(2) * (t34 ^ 2 + t49 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t49) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t53) / 0.2e1 + m(5) * (t14 ^ 2 + t6 ^ 2 + t54) / 0.2e1 + m(4) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t54) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t32 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t40 / 0.2e1) * t40 + (-t4 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t27 / 0.2e1) * t27 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t16 / 0.2e1) * t16 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) + t32 * mrSges(2,3) + Ifges(2,5) * t40 + Ifges(2,1) * t31 / 0.2e1) * t31 + (t14 * mrSges(5,2) + t4 * mrSges(5,3) + Ifges(5,5) * t27 + Ifges(5,1) * t18 / 0.2e1) * t18 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t16 + Ifges(6,1) * t11 / 0.2e1) * t11 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t31 + Ifges(2,6) * t40 + Ifges(2,2) * t30 / 0.2e1) * t30 + (-t14 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t27 + Ifges(5,2) * t17 / 0.2e1) * t17 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t11 + Ifges(6,6) * t16 + Ifges(6,2) * t10 / 0.2e1) * t10 + (t24 * mrSges(3,1) + t15 * mrSges(4,1) - t25 * mrSges(3,2) - t19 * mrSges(4,2) + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t29) * t29 + (t32 * mrSges(3,2) + t20 * mrSges(4,2) - t24 * mrSges(3,3) - t15 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t23 + (Ifges(3,5) + Ifges(4,5)) * t29) * t23 + (-t32 * mrSges(3,1) - t20 * mrSges(4,1) + t25 * mrSges(3,3) + t19 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t22 + (Ifges(3,6) + Ifges(4,6)) * t29 + (Ifges(3,4) + Ifges(4,4)) * t23) * t22;
T  = t5;
