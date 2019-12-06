% Calculate kinetic energy for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:35
% EndTime: 2019-12-05 15:06:35
% DurationCPUTime: 0.70s
% Computational Cost: add. (1303->128), mult. (2020->175), div. (0->0), fcn. (1528->8), ass. (0->42)
t48 = sin(pkin(7));
t50 = cos(pkin(7));
t39 = t48 * V_base(5) + t50 * V_base(4);
t47 = sin(pkin(8));
t49 = cos(pkin(8));
t33 = -t39 * t47 + t49 * V_base(6);
t34 = t39 * t49 + t47 * V_base(6);
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t24 = t33 * t54 - t52 * t34;
t25 = t52 * t33 + t34 * t54;
t43 = V_base(5) * qJ(1) + V_base(1);
t44 = -V_base(4) * qJ(1) + V_base(2);
t35 = -t48 * t43 + t44 * t50;
t31 = -V_base(6) * pkin(1) + qJD(2) - t35;
t26 = -pkin(2) * t33 + t31;
t13 = -pkin(3) * t24 - pkin(6) * t25 + t26;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t38 = t48 * V_base(4) - t50 * V_base(5);
t46 = V_base(3) + qJD(1);
t29 = pkin(1) * t38 - qJ(2) * t39 + t46;
t36 = t50 * t43 + t48 * t44;
t32 = V_base(6) * qJ(2) + t36;
t21 = t49 * t29 - t32 * t47;
t15 = pkin(2) * t38 - pkin(5) * t34 + t21;
t22 = t47 * t29 + t49 * t32;
t18 = pkin(5) * t33 + t22;
t10 = t52 * t15 + t54 * t18;
t37 = qJD(3) + t38;
t8 = pkin(6) * t37 + t10;
t4 = t51 * t13 + t53 * t8;
t3 = t53 * t13 - t51 * t8;
t9 = t15 * t54 - t52 * t18;
t7 = -t37 * pkin(3) - t9;
t23 = qJD(4) - t24;
t20 = t25 * t53 + t37 * t51;
t19 = -t25 * t51 + t37 * t53;
t5 = -t19 * pkin(4) + qJD(5) + t7;
t2 = qJ(5) * t19 + t4;
t1 = pkin(4) * t23 - qJ(5) * t20 + t3;
t6 = m(2) * (t35 ^ 2 + t36 ^ 2 + t46 ^ 2) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t26 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t46 * mrSges(2,2) - t35 * mrSges(2,3) + Ifges(2,1) * t39 / 0.2e1) * t39 + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + Ifges(4,3) * t37 / 0.2e1) * t37 + (t31 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,1) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t31 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t34 + Ifges(3,2) * t33 / 0.2e1) * t33 + (t26 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,5) * t37 + Ifges(4,1) * t25 / 0.2e1) * t25 + (-t26 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,6) * t37 + Ifges(4,2) * t24 / 0.2e1) * t24 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t23) * t23 + (t7 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t20 + (Ifges(5,5) + Ifges(6,5)) * t23) * t20 + (V_base(2) * mrSges(1,1) + t35 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t36 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t39 + Ifges(1,6) * V_base(5) + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t46 * mrSges(2,1) + t21 * mrSges(3,1) - t22 * mrSges(3,2) - t36 * mrSges(2,3) - Ifges(2,4) * t39 + Ifges(3,5) * t34 - Ifges(2,6) * V_base(6) + Ifges(3,6) * t33 + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t38) * t38 + (-t7 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t19 + (Ifges(5,6) + Ifges(6,6)) * t23 + (Ifges(5,4) + Ifges(6,4)) * t20) * t19;
T = t6;
