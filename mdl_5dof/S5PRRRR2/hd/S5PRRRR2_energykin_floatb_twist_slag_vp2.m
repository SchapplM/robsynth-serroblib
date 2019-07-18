% Calculate kinetic energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:08
% EndTime: 2019-07-18 13:30:09
% DurationCPUTime: 0.73s
% Computational Cost: add. (1069->124), mult. (1512->174), div. (0->0), fcn. (1044->8), ass. (0->43)
t39 = -V_base(4) * qJ(1) + V_base(2);
t34 = V_base(6) * pkin(1) + t39;
t38 = V_base(5) * qJ(1) + V_base(1);
t46 = sin(qJ(2));
t50 = cos(qJ(2));
t27 = t50 * t34 - t38 * t46;
t33 = t46 * V_base(5) + t50 * V_base(4);
t41 = V_base(6) + qJD(2);
t22 = pkin(2) * t41 - pkin(4) * t33 + t27;
t28 = t46 * t34 + t50 * t38;
t32 = -t46 * V_base(4) + t50 * V_base(5);
t24 = pkin(4) * t32 + t28;
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t13 = t45 * t22 + t49 * t24;
t25 = t32 * t49 - t33 * t45;
t11 = pkin(5) * t25 + t13;
t44 = sin(qJ(4));
t48 = cos(qJ(4));
t12 = t49 * t22 - t24 * t45;
t26 = t32 * t45 + t33 * t49;
t40 = qJD(3) + t41;
t51 = pkin(3) * t40 - pkin(5) * t26 + t12;
t4 = t11 * t44 - t48 * t51;
t52 = t4 ^ 2;
t6 = t48 * t11 + t44 * t51;
t42 = V_base(3) + qJD(1);
t17 = t25 * t48 - t26 * t44;
t36 = -V_base(5) * pkin(1) + t42;
t29 = -pkin(2) * t32 + t36;
t19 = -pkin(3) * t25 + t29;
t47 = cos(qJ(5));
t43 = sin(qJ(5));
t37 = qJD(4) + t40;
t18 = t25 * t44 + t26 * t48;
t16 = qJD(5) - t17;
t15 = t18 * t47 + t37 * t43;
t14 = -t18 * t43 + t37 * t47;
t10 = -pkin(6) * t18 + t19;
t3 = pkin(6) * t37 + t6;
t2 = t10 * t43 + t3 * t47;
t1 = t10 * t47 - t3 * t43;
t5 = m(2) * (t38 ^ 2 + t39 ^ 2 + t42 ^ 2) / 0.2e1 + m(3) * (t27 ^ 2 + t28 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t29 ^ 2) / 0.2e1 + m(5) * (t19 ^ 2 + t6 ^ 2 + t52) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t52) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t27 * mrSges(3,1) - t28 * mrSges(3,2) + Ifges(3,3) * t41 / 0.2e1) * t41 + (t12 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,3) * t40 / 0.2e1) * t40 + (-t4 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t37 / 0.2e1) * t37 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t16 / 0.2e1) * t16 + (t36 * mrSges(3,2) - t27 * mrSges(3,3) + Ifges(3,5) * t41 + Ifges(3,1) * t33 / 0.2e1) * t33 + (t29 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,5) * t40 + Ifges(4,1) * t26 / 0.2e1) * t26 + (t19 * mrSges(5,2) + t4 * mrSges(5,3) + Ifges(5,5) * t37 + Ifges(5,1) * t18 / 0.2e1) * t18 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t16 + Ifges(6,1) * t15 / 0.2e1) * t15 + (-t36 * mrSges(3,1) + t28 * mrSges(3,3) + Ifges(3,4) * t33 + Ifges(3,6) * t41 + Ifges(3,2) * t32 / 0.2e1) * t32 + (-t29 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t26 + Ifges(4,6) * t40 + Ifges(4,2) * t25 / 0.2e1) * t25 + (-t19 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t37 + Ifges(5,2) * t17 / 0.2e1) * t17 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t16 + Ifges(6,2) * t14 / 0.2e1) * t14 + (V_base(2) * mrSges(1,1) + t39 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t38 * mrSges(2,2) + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) - t42 * mrSges(2,1) + V_base(1) * mrSges(1,3) + t38 * mrSges(2,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(5) + (Ifges(1,6) + Ifges(2,6)) * V_base(6)) * V_base(5) + (V_base(3) * mrSges(1,2) + t42 * mrSges(2,2) - V_base(2) * mrSges(1,3) - t39 * mrSges(2,3) + (Ifges(2,1) / 0.2e1 + Ifges(1,1) / 0.2e1) * V_base(4) + (Ifges(1,5) + Ifges(2,5)) * V_base(6) + (Ifges(1,4) + Ifges(2,4)) * V_base(5)) * V_base(4);
T  = t5;
