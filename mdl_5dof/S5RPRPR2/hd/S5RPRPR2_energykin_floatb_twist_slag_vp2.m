% Calculate kinetic energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:17
% EndTime: 2019-12-05 17:49:18
% DurationCPUTime: 0.86s
% Computational Cost: add. (2035->132), mult. (3094->191), div. (0->0), fcn. (2460->10), ass. (0->49)
t45 = V_base(6) * pkin(5) + V_base(2);
t46 = -V_base(5) * pkin(5) + V_base(3);
t55 = sin(qJ(1));
t57 = cos(qJ(1));
t40 = -t45 * t57 - t46 * t55;
t42 = -t55 * V_base(5) + t57 * V_base(6);
t48 = V_base(4) + qJD(1);
t33 = pkin(1) * t48 - qJ(2) * t42 + t40;
t39 = -t45 * t55 + t57 * t46;
t43 = -t55 * V_base(6) - t57 * V_base(5);
t36 = qJ(2) * t43 + t39;
t50 = sin(pkin(8));
t52 = cos(pkin(8));
t25 = t52 * t33 - t36 * t50;
t38 = t42 * t52 + t43 * t50;
t19 = pkin(2) * t48 - pkin(6) * t38 + t25;
t26 = t50 * t33 + t52 * t36;
t37 = -t42 * t50 + t43 * t52;
t22 = pkin(6) * t37 + t26;
t54 = sin(qJ(3));
t59 = cos(qJ(3));
t12 = t54 * t19 + t59 * t22;
t47 = qJD(3) + t48;
t10 = qJ(4) * t47 + t12;
t28 = -t59 * t37 + t38 * t54;
t29 = t54 * t37 + t59 * t38;
t41 = -pkin(1) * t43 + qJD(2) + V_base(1);
t30 = -pkin(2) * t37 + t41;
t15 = pkin(3) * t28 - qJ(4) * t29 + t30;
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t6 = t51 * t10 + t49 * t15;
t5 = -t10 * t49 + t51 * t15;
t11 = t59 * t19 - t54 * t22;
t9 = -t47 * pkin(3) + qJD(4) - t11;
t58 = V_base(1) ^ 2;
t56 = cos(qJ(5));
t53 = sin(qJ(5));
t27 = qJD(5) + t28;
t24 = t29 * t51 + t47 * t49;
t23 = -t29 * t49 + t47 * t51;
t17 = t23 * t53 + t24 * t56;
t16 = t23 * t56 - t24 * t53;
t7 = -t23 * pkin(4) + t9;
t4 = pkin(7) * t23 + t6;
t3 = pkin(4) * t28 - pkin(7) * t24 + t5;
t2 = t3 * t53 + t4 * t56;
t1 = t3 * t56 - t4 * t53;
t8 = m(2) * (t39 ^ 2 + t40 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t58) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t41 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t30 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t47 / 0.2e1) * t47 + (-V_base(1) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,2) * t43 / 0.2e1) * t43 + (t41 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,1) * t38 / 0.2e1) * t38 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t27 / 0.2e1) * t27 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(1) * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,4) * t43 + Ifges(2,1) * t42 / 0.2e1) * t42 + (-t41 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,2) * t37 / 0.2e1) * t37 + (t30 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t47 + Ifges(4,1) * t29 / 0.2e1) * t29 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,2) * t23 / 0.2e1) * t23 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t27 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t27 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t40 * mrSges(2,1) + t25 * mrSges(3,1) - t39 * mrSges(2,2) - t26 * mrSges(3,2) + Ifges(2,5) * t42 + Ifges(3,5) * t38 + Ifges(2,6) * t43 + Ifges(3,6) * t37 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t48) * t48 + (t30 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t12 * mrSges(4,3) - Ifges(4,4) * t29 + Ifges(5,5) * t24 - Ifges(4,6) * t47 + Ifges(5,6) * t23 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t28) * t28;
T = t8;
