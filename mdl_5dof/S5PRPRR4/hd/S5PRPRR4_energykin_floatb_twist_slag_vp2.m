% Calculate kinetic energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:31
% EndTime: 2019-12-05 15:49:32
% DurationCPUTime: 1.08s
% Computational Cost: add. (2745->137), mult. (4920->202), div. (0->0), fcn. (4092->12), ass. (0->55)
t51 = V_base(5) * qJ(1) + V_base(1);
t52 = -V_base(4) * qJ(1) + V_base(2);
t55 = sin(pkin(9));
t58 = cos(pkin(9));
t44 = -t51 * t55 + t58 * t52;
t59 = cos(pkin(5));
t47 = t55 * V_base(5) + t58 * V_base(4);
t69 = pkin(6) * t47;
t38 = V_base(6) * pkin(1) - t59 * t69 + t44;
t46 = -t55 * V_base(4) + t58 * V_base(5);
t53 = V_base(3) + qJD(1);
t56 = sin(pkin(5));
t42 = -pkin(1) * t46 - t56 * t69 + t53;
t70 = t38 * t59 + t42 * t56;
t45 = t58 * t51 + t55 * t52;
t66 = t46 * t59 + t56 * V_base(6);
t35 = t66 * pkin(6) + t45;
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t22 = -t35 * t62 + t70 * t65;
t37 = t47 * t65 + t66 * t62;
t43 = -t46 * t56 + t59 * V_base(6) + qJD(2);
t18 = pkin(2) * t43 - qJ(3) * t37 + t22;
t23 = t65 * t35 + t70 * t62;
t36 = -t62 * t47 + t66 * t65;
t21 = qJ(3) * t36 + t23;
t54 = sin(pkin(10));
t57 = cos(pkin(10));
t12 = t54 * t18 + t57 * t21;
t10 = pkin(7) * t43 + t12;
t31 = -t38 * t56 + t59 * t42;
t24 = -pkin(2) * t36 + qJD(3) + t31;
t29 = t36 * t57 - t37 * t54;
t30 = t36 * t54 + t37 * t57;
t14 = -pkin(3) * t29 - pkin(7) * t30 + t24;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t6 = t64 * t10 + t61 * t14;
t11 = t18 * t57 - t54 * t21;
t5 = -t10 * t61 + t14 * t64;
t26 = -t30 * t61 + t43 * t64;
t9 = -pkin(3) * t43 - t11;
t63 = cos(qJ(5));
t60 = sin(qJ(5));
t28 = qJD(4) - t29;
t27 = t30 * t64 + t43 * t61;
t25 = qJD(5) - t26;
t16 = t27 * t63 + t28 * t60;
t15 = -t27 * t60 + t28 * t63;
t7 = -pkin(4) * t26 - pkin(8) * t27 + t9;
t4 = pkin(8) * t28 + t6;
t3 = -pkin(4) * t28 - t5;
t2 = t4 * t63 + t60 * t7;
t1 = -t4 * t60 + t63 * t7;
t8 = m(2) * (t44 ^ 2 + t45 ^ 2 + t53 ^ 2) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t53 * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,1) * t47 / 0.2e1) * t47 + (t31 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t24 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t53 * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t47 + Ifges(2,2) * t46 / 0.2e1) * t46 + (-t31 * mrSges(3,1) + t23 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t24 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t27 / 0.2e1) * t27 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t25 + Ifges(6,1) * t16 / 0.2e1) * t16 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t27 + Ifges(5,6) * t28 + Ifges(5,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t25 + Ifges(6,2) * t15 / 0.2e1) * t15 + (V_base(2) * mrSges(1,1) + t44 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t45 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t47 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t46 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t22 * mrSges(3,1) + t11 * mrSges(4,1) - t23 * mrSges(3,2) - t12 * mrSges(4,2) + Ifges(3,5) * t37 + Ifges(4,5) * t30 + Ifges(3,6) * t36 + Ifges(4,6) * t29 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t43) * t43;
T = t8;
