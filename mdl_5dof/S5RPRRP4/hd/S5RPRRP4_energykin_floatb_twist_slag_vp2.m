% Calculate kinetic energy for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:24
% EndTime: 2019-12-05 18:05:25
% DurationCPUTime: 0.73s
% Computational Cost: add. (1505->128), mult. (2052->176), div. (0->0), fcn. (1552->8), ass. (0->43)
t44 = V_base(6) * pkin(5) + V_base(2);
t45 = -V_base(5) * pkin(5) + V_base(3);
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t37 = -t54 * t44 - t51 * t45;
t46 = V_base(4) + qJD(1);
t30 = -pkin(1) * t46 + qJD(2) - t37;
t39 = -t51 * V_base(5) + t54 * V_base(6);
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t34 = -t47 * t39 + t46 * t48;
t35 = t39 * t48 + t46 * t47;
t19 = -pkin(2) * t34 - pkin(6) * t35 + t30;
t40 = t51 * V_base(6) + t54 * V_base(5);
t27 = pkin(1) * t40 - qJ(2) * t39 + V_base(1);
t36 = -t44 * t51 + t54 * t45;
t32 = qJ(2) * t46 + t36;
t23 = t47 * t27 + t48 * t32;
t21 = pkin(6) * t40 + t23;
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t12 = t50 * t19 + t53 * t21;
t24 = -t35 * t50 + t40 * t53;
t10 = pkin(7) * t24 + t12;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t11 = t53 * t19 - t21 * t50;
t25 = t35 * t53 + t40 * t50;
t33 = qJD(3) - t34;
t7 = pkin(3) * t33 - pkin(7) * t25 + t11;
t4 = t52 * t10 + t49 * t7;
t3 = -t10 * t49 + t52 * t7;
t22 = t27 * t48 - t47 * t32;
t20 = -pkin(2) * t40 - t22;
t13 = -pkin(3) * t24 + t20;
t55 = V_base(1) ^ 2;
t31 = qJD(4) + t33;
t15 = t24 * t49 + t25 * t52;
t14 = t24 * t52 - t25 * t49;
t8 = -pkin(4) * t14 + qJD(5) + t13;
t2 = qJ(5) * t14 + t4;
t1 = pkin(4) * t31 - qJ(5) * t15 + t3;
t5 = m(2) * (t36 ^ 2 + t37 ^ 2 + t55) / 0.2e1 + m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t55) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t30 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t20 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t37 * mrSges(2,1) - t36 * mrSges(2,2) + Ifges(2,3) * t46 / 0.2e1) * t46 + (t30 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,1) * t35 / 0.2e1) * t35 + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t33 / 0.2e1) * t33 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(1) * mrSges(2,2) - t37 * mrSges(2,3) + Ifges(2,5) * t46 + Ifges(2,1) * t39 / 0.2e1) * t39 + (-t30 * mrSges(3,1) + t23 * mrSges(3,3) + Ifges(3,4) * t35 + Ifges(3,2) * t34 / 0.2e1) * t34 + (t20 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t33 + Ifges(4,1) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t20 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,6) * t33 + Ifges(4,2) * t24 / 0.2e1) * t24 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t31) * t31 + (t13 * mrSges(5,2) + t8 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t15 + (Ifges(5,5) + Ifges(6,5)) * t31) * t15 + (V_base(1) * mrSges(2,1) + t22 * mrSges(3,1) - t23 * mrSges(3,2) - t36 * mrSges(2,3) - Ifges(2,4) * t39 + Ifges(3,5) * t35 - Ifges(2,6) * t46 + Ifges(3,6) * t34 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t40) * t40 + (-t13 * mrSges(5,1) - t8 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t14 + (Ifges(5,6) + Ifges(6,6)) * t31 + (Ifges(5,4) + Ifges(6,4)) * t15) * t14;
T = t5;
