% Calculate kinetic energy for
% S5RRRRP2
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:00
% EndTime: 2020-01-03 12:11:01
% DurationCPUTime: 0.76s
% Computational Cost: add. (1631->128), mult. (2224->178), div. (0->0), fcn. (1708->8), ass. (0->44)
t44 = V_base(6) * pkin(5) + V_base(2);
t45 = -V_base(5) * pkin(5) + V_base(3);
t51 = sin(qJ(1));
t55 = cos(qJ(1));
t36 = t55 * t44 + t51 * t45;
t39 = t51 * V_base(5) - t55 * V_base(6);
t47 = V_base(4) + qJD(1);
t27 = pkin(1) * t47 - pkin(6) * t39 + t36;
t35 = t51 * t44 - t45 * t55;
t40 = t51 * V_base(6) + t55 * V_base(5);
t31 = pkin(6) * t40 + t35;
t50 = sin(qJ(2));
t54 = cos(qJ(2));
t23 = t50 * t27 + t54 * t31;
t46 = qJD(2) + t47;
t18 = pkin(7) * t46 + t23;
t33 = -t50 * t39 + t40 * t54;
t34 = t39 * t54 + t40 * t50;
t37 = -pkin(1) * t40 + V_base(1);
t21 = -pkin(2) * t33 - pkin(7) * t34 + t37;
t49 = sin(qJ(3));
t53 = cos(qJ(3));
t12 = t53 * t18 + t49 * t21;
t24 = -t34 * t49 + t46 * t53;
t10 = pkin(8) * t24 + t12;
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t11 = -t18 * t49 + t53 * t21;
t25 = t34 * t53 + t46 * t49;
t32 = qJD(3) - t33;
t7 = pkin(3) * t32 - pkin(8) * t25 + t11;
t4 = t52 * t10 + t48 * t7;
t3 = -t10 * t48 + t52 * t7;
t22 = t27 * t54 - t50 * t31;
t17 = -pkin(2) * t46 - t22;
t13 = -pkin(3) * t24 + t17;
t56 = V_base(1) ^ 2;
t30 = qJD(4) + t32;
t15 = t24 * t48 + t25 * t52;
t14 = t24 * t52 - t25 * t48;
t8 = -pkin(4) * t14 + qJD(5) + t13;
t2 = qJ(5) * t14 + t4;
t1 = pkin(4) * t30 - qJ(5) * t15 + t3;
t5 = m(2) * (t35 ^ 2 + t36 ^ 2 + t56) / 0.2e1 + m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t56) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t37 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t17 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t36 * mrSges(2,1) - t35 * mrSges(2,2) + Ifges(2,3) * t47 / 0.2e1) * t47 + (t22 * mrSges(3,1) - t23 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(1) * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,6) * t47 + Ifges(2,2) * t40 / 0.2e1) * t40 + (t37 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t34 / 0.2e1) * t34 + (t17 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t32 + Ifges(4,1) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (V_base(1) * mrSges(2,2) - t36 * mrSges(2,3) + Ifges(2,4) * t40 + Ifges(2,5) * t47 + Ifges(2,1) * t39 / 0.2e1) * t39 + (-t37 * mrSges(3,1) + t23 * mrSges(3,3) + Ifges(3,4) * t34 + Ifges(3,6) * t46 + Ifges(3,2) * t33 / 0.2e1) * t33 + (-t17 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,6) * t32 + Ifges(4,2) * t24 / 0.2e1) * t24 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t30) * t30 + (t13 * mrSges(5,2) + t8 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t15 + (Ifges(5,5) + Ifges(6,5)) * t30) * t15 + (-t13 * mrSges(5,1) - t8 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t14 + (Ifges(5,6) + Ifges(6,6)) * t30 + (Ifges(5,4) + Ifges(6,4)) * t15) * t14;
T = t5;
