% Calculate kinetic energy for
% S5RRPRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:41
% EndTime: 2019-12-31 20:23:42
% DurationCPUTime: 1.10s
% Computational Cost: add. (3131->137), mult. (4920->203), div. (0->0), fcn. (4092->12), ass. (0->56)
t51 = V_base(5) * pkin(6) + V_base(1);
t52 = -V_base(4) * pkin(6) + V_base(2);
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t44 = -t51 * t61 + t65 * t52;
t53 = V_base(6) + qJD(1);
t57 = cos(pkin(5));
t47 = t61 * V_base(5) + t65 * V_base(4);
t70 = pkin(7) * t47;
t38 = pkin(1) * t53 - t57 * t70 + t44;
t46 = -t61 * V_base(4) + t65 * V_base(5);
t55 = sin(pkin(5));
t42 = -pkin(1) * t46 - t55 * t70 + V_base(3);
t71 = t38 * t57 + t42 * t55;
t45 = t65 * t51 + t61 * t52;
t67 = t46 * t57 + t53 * t55;
t35 = t67 * pkin(7) + t45;
t60 = sin(qJ(2));
t64 = cos(qJ(2));
t22 = -t35 * t60 + t71 * t64;
t37 = t47 * t64 + t67 * t60;
t43 = -t46 * t55 + t53 * t57 + qJD(2);
t18 = pkin(2) * t43 - qJ(3) * t37 + t22;
t23 = t64 * t35 + t71 * t60;
t36 = -t47 * t60 + t67 * t64;
t21 = qJ(3) * t36 + t23;
t54 = sin(pkin(10));
t56 = cos(pkin(10));
t12 = t54 * t18 + t56 * t21;
t10 = pkin(8) * t43 + t12;
t31 = -t38 * t55 + t57 * t42;
t25 = -pkin(2) * t36 + qJD(3) + t31;
t29 = t36 * t56 - t37 * t54;
t30 = t36 * t54 + t37 * t56;
t14 = -pkin(3) * t29 - pkin(8) * t30 + t25;
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t6 = t63 * t10 + t59 * t14;
t11 = t18 * t56 - t54 * t21;
t5 = -t10 * t59 + t14 * t63;
t26 = -t30 * t59 + t43 * t63;
t9 = -pkin(3) * t43 - t11;
t66 = V_base(3) ^ 2;
t62 = cos(qJ(5));
t58 = sin(qJ(5));
t28 = qJD(4) - t29;
t27 = t30 * t63 + t43 * t59;
t24 = qJD(5) - t26;
t16 = t27 * t62 + t28 * t58;
t15 = -t27 * t58 + t28 * t62;
t7 = -pkin(4) * t26 - pkin(9) * t27 + t9;
t4 = pkin(9) * t28 + t6;
t3 = -pkin(4) * t28 - t5;
t2 = t4 * t62 + t58 * t7;
t1 = -t4 * t58 + t62 * t7;
t8 = m(2) * (t44 ^ 2 + t45 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,3) * t53 / 0.2e1) * t53 + (t31 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t25 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,5) * t53 + Ifges(2,1) * t47 / 0.2e1) * t47 + (-t31 * mrSges(3,1) + t23 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t25 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t27 / 0.2e1) * t27 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t24 + Ifges(6,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t47 + Ifges(2,6) * t53 + Ifges(2,2) * t46 / 0.2e1) * t46 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t27 + Ifges(5,6) * t28 + Ifges(5,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t24 + Ifges(6,2) * t15 / 0.2e1) * t15 + (t22 * mrSges(3,1) + t11 * mrSges(4,1) - t23 * mrSges(3,2) - t12 * mrSges(4,2) + Ifges(3,5) * t37 + Ifges(4,5) * t30 + Ifges(3,6) * t36 + Ifges(4,6) * t29 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t43) * t43;
T = t8;
