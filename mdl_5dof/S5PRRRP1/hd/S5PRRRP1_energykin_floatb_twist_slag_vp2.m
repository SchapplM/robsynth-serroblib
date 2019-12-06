% Calculate kinetic energy for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:42
% EndTime: 2019-12-05 16:39:43
% DurationCPUTime: 0.80s
% Computational Cost: add. (1549->128), mult. (2396->177), div. (0->0), fcn. (1856->8), ass. (0->43)
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t37 = -t47 * V_base(4) + t48 * V_base(5);
t38 = t47 * V_base(5) + t48 * V_base(4);
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t32 = t37 * t54 - t51 * t38;
t33 = t51 * t37 + t38 * t54;
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t24 = t32 * t53 - t33 * t50;
t25 = t32 * t50 + t33 * t53;
t46 = V_base(3) + qJD(1);
t36 = -pkin(1) * t37 + t46;
t26 = -pkin(2) * t32 + t36;
t13 = -pkin(3) * t24 - pkin(7) * t25 + t26;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t42 = V_base(5) * qJ(1) + V_base(1);
t43 = -V_base(4) * qJ(1) + V_base(2);
t34 = -t42 * t47 + t48 * t43;
t29 = V_base(6) * pkin(1) - pkin(5) * t38 + t34;
t35 = t48 * t42 + t47 * t43;
t31 = pkin(5) * t37 + t35;
t21 = t54 * t29 - t31 * t51;
t45 = V_base(6) + qJD(2);
t15 = pkin(2) * t45 - pkin(6) * t33 + t21;
t22 = t51 * t29 + t54 * t31;
t18 = pkin(6) * t32 + t22;
t10 = t50 * t15 + t53 * t18;
t44 = qJD(3) + t45;
t8 = pkin(7) * t44 + t10;
t4 = t49 * t13 + t52 * t8;
t3 = t52 * t13 - t49 * t8;
t9 = t15 * t53 - t50 * t18;
t7 = -pkin(3) * t44 - t9;
t23 = qJD(4) - t24;
t20 = t25 * t52 + t44 * t49;
t19 = -t25 * t49 + t44 * t52;
t5 = -pkin(4) * t19 + qJD(5) + t7;
t2 = qJ(5) * t19 + t4;
t1 = pkin(4) * t23 - qJ(5) * t20 + t3;
t6 = m(2) * (t34 ^ 2 + t35 ^ 2 + t46 ^ 2) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t26 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t21 * mrSges(3,1) - t22 * mrSges(3,2) + Ifges(3,3) * t45 / 0.2e1) * t45 + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + Ifges(4,3) * t44 / 0.2e1) * t44 + (t46 * mrSges(2,2) - t34 * mrSges(2,3) + Ifges(2,1) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t46 * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,2) * t37 / 0.2e1) * t37 + (t36 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,5) * t45 + Ifges(3,1) * t33 / 0.2e1) * t33 + (t26 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,5) * t44 + Ifges(4,1) * t25 / 0.2e1) * t25 + (-t36 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t33 + Ifges(3,6) * t45 + Ifges(3,2) * t32 / 0.2e1) * t32 + (-t26 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,6) * t44 + Ifges(4,2) * t24 / 0.2e1) * t24 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t23) * t23 + (t7 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t20 + (Ifges(5,5) + Ifges(6,5)) * t23) * t20 + (V_base(2) * mrSges(1,1) + t34 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t35 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t38 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t37 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t7 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t19 + (Ifges(5,6) + Ifges(6,6)) * t23 + (Ifges(5,4) + Ifges(6,4)) * t20) * t19;
T = t6;
