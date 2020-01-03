% Calculate kinetic energy for
% S5RPRPR3
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:52
% EndTime: 2020-01-03 11:35:53
% DurationCPUTime: 0.91s
% Computational Cost: add. (1965->132), mult. (2994->191), div. (0->0), fcn. (2368->10), ass. (0->49)
t46 = V_base(6) * pkin(5) + V_base(2);
t47 = -V_base(5) * pkin(5) + V_base(3);
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t39 = t58 * t46 + t56 * t47;
t41 = t56 * V_base(5) - t58 * V_base(6);
t49 = V_base(4) + qJD(1);
t32 = pkin(1) * t49 - qJ(2) * t41 + t39;
t38 = t56 * t46 - t47 * t58;
t42 = t56 * V_base(6) + t58 * V_base(5);
t35 = qJ(2) * t42 + t38;
t51 = sin(pkin(8));
t53 = cos(pkin(8));
t25 = t53 * t32 - t35 * t51;
t37 = t41 * t53 + t42 * t51;
t18 = pkin(2) * t49 - pkin(6) * t37 + t25;
t26 = t51 * t32 + t53 * t35;
t36 = -t41 * t51 + t42 * t53;
t21 = pkin(6) * t36 + t26;
t55 = sin(qJ(3));
t60 = cos(qJ(3));
t12 = t55 * t18 + t60 * t21;
t48 = qJD(3) + t49;
t10 = qJ(4) * t48 + t12;
t27 = -t60 * t36 + t37 * t55;
t28 = t55 * t36 + t60 * t37;
t40 = -pkin(1) * t42 + qJD(2) + V_base(1);
t29 = -pkin(2) * t36 + t40;
t14 = pkin(3) * t27 - qJ(4) * t28 + t29;
t50 = sin(pkin(9));
t52 = cos(pkin(9));
t6 = t52 * t10 + t50 * t14;
t11 = t60 * t18 - t55 * t21;
t5 = -t10 * t50 + t14 * t52;
t23 = -t28 * t50 + t48 * t52;
t9 = -t48 * pkin(3) + qJD(4) - t11;
t59 = V_base(1) ^ 2;
t57 = cos(qJ(5));
t54 = sin(qJ(5));
t24 = t28 * t52 + t48 * t50;
t22 = qJD(5) - t23;
t16 = t24 * t57 + t27 * t54;
t15 = -t24 * t54 + t27 * t57;
t7 = -t23 * pkin(4) - t24 * pkin(7) + t9;
t4 = pkin(7) * t27 + t6;
t3 = -pkin(4) * t27 - t5;
t2 = t4 * t57 + t54 * t7;
t1 = -t4 * t54 + t57 * t7;
t8 = m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t59) / 0.2e1 + m(2) * (t38 ^ 2 + t39 ^ 2 + t59) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t29 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (-V_base(1) * mrSges(2,1) + t38 * mrSges(2,3) + Ifges(2,2) * t42 / 0.2e1) * t42 + (t40 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t24 / 0.2e1) * t24 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t22 / 0.2e1) * t22 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(1) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,4) * t42 + Ifges(2,1) * t41 / 0.2e1) * t41 + (-t40 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,2) * t36 / 0.2e1) * t36 + (t29 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t28 / 0.2e1) * t28 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,2) * t23 / 0.2e1) * t23 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t22 + Ifges(6,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t22 + Ifges(6,2) * t15 / 0.2e1) * t15 + (t39 * mrSges(2,1) + t25 * mrSges(3,1) - t38 * mrSges(2,2) - t26 * mrSges(3,2) + Ifges(2,5) * t41 + Ifges(3,5) * t37 + Ifges(2,6) * t42 + Ifges(3,6) * t36 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t49) * t49 + (t29 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t12 * mrSges(4,3) - Ifges(4,4) * t28 + Ifges(5,5) * t24 - Ifges(4,6) * t48 + Ifges(5,6) * t23 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t27) * t27;
T = t8;
