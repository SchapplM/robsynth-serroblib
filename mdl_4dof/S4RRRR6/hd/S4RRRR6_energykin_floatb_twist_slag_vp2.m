% Calculate kinetic energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:18
% EndTime: 2019-12-31 17:29:19
% DurationCPUTime: 0.75s
% Computational Cost: add. (1517->113), mult. (2285->170), div. (0->0), fcn. (1810->10), ass. (0->47)
t40 = V_base(5) * pkin(5) + V_base(1);
t41 = -V_base(4) * pkin(5) + V_base(2);
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t33 = -t40 * t48 + t52 * t41;
t42 = V_base(6) + qJD(1);
t44 = cos(pkin(4));
t36 = t48 * V_base(5) + t52 * V_base(4);
t59 = pkin(6) * t36;
t28 = pkin(1) * t42 - t44 * t59 + t33;
t35 = -t48 * V_base(4) + t52 * V_base(5);
t43 = sin(pkin(4));
t31 = -pkin(1) * t35 - t43 * t59 + V_base(3);
t60 = t28 * t44 + t31 * t43;
t34 = t52 * t40 + t48 * t41;
t54 = t35 * t44 + t42 * t43;
t25 = t54 * pkin(6) + t34;
t47 = sin(qJ(2));
t51 = cos(qJ(2));
t13 = -t47 * t25 + t51 * t60;
t26 = -t36 * t47 + t54 * t51;
t14 = t51 * t25 + t47 * t60;
t32 = -t35 * t43 + t42 * t44 + qJD(2);
t12 = pkin(7) * t32 + t14;
t46 = sin(qJ(3));
t50 = cos(qJ(3));
t17 = -t28 * t43 + t44 * t31;
t27 = t36 * t51 + t54 * t47;
t9 = -pkin(2) * t26 - pkin(7) * t27 + t17;
t6 = t50 * t12 + t46 * t9;
t5 = -t12 * t46 + t50 * t9;
t19 = -t27 * t46 + t32 * t50;
t11 = -pkin(2) * t32 - t13;
t53 = V_base(3) ^ 2;
t49 = cos(qJ(4));
t45 = sin(qJ(4));
t24 = qJD(3) - t26;
t20 = t27 * t50 + t32 * t46;
t18 = qJD(4) - t19;
t16 = t20 * t49 + t24 * t45;
t15 = -t20 * t45 + t24 * t49;
t7 = -pkin(3) * t19 - pkin(8) * t20 + t11;
t4 = pkin(8) * t24 + t6;
t3 = -pkin(3) * t24 - t5;
t2 = t4 * t49 + t45 * t7;
t1 = -t4 * t45 + t49 * t7;
t8 = m(2) * (t33 ^ 2 + t34 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t17 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t42 / 0.2e1) * t42 + (t13 * mrSges(3,1) - t14 * mrSges(3,2) + Ifges(3,3) * t32 / 0.2e1) * t32 + (t5 * mrSges(4,1) - t6 * mrSges(4,2) + Ifges(4,3) * t24 / 0.2e1) * t24 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,3) * t18 / 0.2e1) * t18 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,5) * t42 + Ifges(2,1) * t36 / 0.2e1) * t36 + (t17 * mrSges(3,2) - t13 * mrSges(3,3) + Ifges(3,5) * t32 + Ifges(3,1) * t27 / 0.2e1) * t27 + (t11 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,5) * t24 + Ifges(4,1) * t20 / 0.2e1) * t20 + (t3 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,5) * t18 + Ifges(5,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t36 + Ifges(2,6) * t42 + Ifges(2,2) * t35 / 0.2e1) * t35 + (-t17 * mrSges(3,1) + t14 * mrSges(3,3) + Ifges(3,4) * t27 + Ifges(3,6) * t32 + Ifges(3,2) * t26 / 0.2e1) * t26 + (-t11 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t20 + Ifges(4,6) * t24 + Ifges(4,2) * t19 / 0.2e1) * t19 + (-t3 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t16 + Ifges(5,6) * t18 + Ifges(5,2) * t15 / 0.2e1) * t15;
T = t8;
