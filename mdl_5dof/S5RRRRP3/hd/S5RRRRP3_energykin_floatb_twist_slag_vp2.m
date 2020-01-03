% Calculate kinetic energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:03
% EndTime: 2019-12-31 21:49:04
% DurationCPUTime: 0.78s
% Computational Cost: add. (1641->128), mult. (2392->178), div. (0->0), fcn. (1852->8), ass. (0->44)
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t36 = -t50 * V_base(4) + t53 * V_base(5);
t37 = t50 * V_base(5) + t53 * V_base(4);
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t31 = t36 * t52 - t37 * t49;
t32 = t36 * t49 + t37 * t52;
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t23 = t31 * t51 - t32 * t48;
t24 = t31 * t48 + t32 * t51;
t35 = -pkin(1) * t36 + V_base(3);
t25 = -pkin(2) * t31 + t35;
t12 = -pkin(3) * t23 - pkin(8) * t24 + t25;
t47 = sin(qJ(4));
t55 = cos(qJ(4));
t43 = V_base(5) * pkin(5) + V_base(1);
t44 = -V_base(4) * pkin(5) + V_base(2);
t33 = -t43 * t50 + t53 * t44;
t46 = V_base(6) + qJD(1);
t28 = pkin(1) * t46 - pkin(6) * t37 + t33;
t34 = t53 * t43 + t50 * t44;
t30 = pkin(6) * t36 + t34;
t20 = t52 * t28 - t30 * t49;
t45 = qJD(2) + t46;
t14 = pkin(2) * t45 - pkin(7) * t32 + t20;
t21 = t49 * t28 + t52 * t30;
t17 = pkin(7) * t31 + t21;
t10 = t48 * t14 + t51 * t17;
t42 = qJD(3) + t45;
t8 = pkin(8) * t42 + t10;
t4 = t47 * t12 + t55 * t8;
t9 = t14 * t51 - t48 * t17;
t7 = -pkin(3) * t42 - t9;
t3 = t55 * t12 - t47 * t8;
t54 = V_base(3) ^ 2;
t22 = qJD(4) - t23;
t19 = t55 * t24 + t47 * t42;
t18 = t24 * t47 - t55 * t42;
t5 = pkin(4) * t18 - qJ(5) * t19 + t7;
t2 = qJ(5) * t22 + t4;
t1 = -t22 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t33 ^ 2 + t34 ^ 2 + t54) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t54) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t35 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t25 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t46 / 0.2e1) * t46 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + Ifges(3,3) * t45 / 0.2e1) * t45 + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + Ifges(4,3) * t42 / 0.2e1) * t42 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,5) * t46 + Ifges(2,1) * t37 / 0.2e1) * t37 + (t35 * mrSges(3,2) - t20 * mrSges(3,3) + Ifges(3,5) * t45 + Ifges(3,1) * t32 / 0.2e1) * t32 + (t25 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,5) * t42 + Ifges(4,1) * t24 / 0.2e1) * t24 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t37 + Ifges(2,6) * t46 + Ifges(2,2) * t36 / 0.2e1) * t36 + (-t35 * mrSges(3,1) + t21 * mrSges(3,3) + Ifges(3,4) * t32 + Ifges(3,6) * t45 + Ifges(3,2) * t31 / 0.2e1) * t31 + (-t25 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t24 + Ifges(4,6) * t42 + Ifges(4,2) * t23 / 0.2e1) * t23 + (t3 * mrSges(5,1) - t1 * mrSges(6,1) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t22) * t22 + (t7 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t19 + (Ifges(6,4) + Ifges(5,5)) * t22) * t19 + (t7 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t18 + (-Ifges(5,6) + Ifges(6,6)) * t22 + (-Ifges(5,4) + Ifges(6,5)) * t19) * t18;
T = t6;
