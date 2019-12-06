% Calculate kinetic energy for
% S5PRPRR6
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:21
% EndTime: 2019-12-05 15:56:22
% DurationCPUTime: 1.06s
% Computational Cost: add. (2811->137), mult. (4866->202), div. (0->0), fcn. (4038->12), ass. (0->55)
t56 = sin(pkin(9));
t59 = cos(pkin(9));
t46 = -t56 * V_base(4) + t59 * V_base(5);
t57 = sin(pkin(5));
t60 = cos(pkin(5));
t67 = t46 * t60 + t57 * V_base(6);
t52 = V_base(5) * qJ(1) + V_base(1);
t53 = -V_base(4) * qJ(1) + V_base(2);
t43 = -t52 * t56 + t59 * t53;
t47 = t56 * V_base(5) + t59 * V_base(4);
t73 = pkin(6) * t47;
t38 = V_base(6) * pkin(1) - t60 * t73 + t43;
t54 = V_base(3) + qJD(1);
t41 = -pkin(1) * t46 - t57 * t73 + t54;
t74 = t38 * t60 + t41 * t57;
t44 = t59 * t52 + t56 * t53;
t34 = t67 * pkin(6) + t44;
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t26 = -t63 * t34 + t66 * t74;
t28 = -t38 * t57 + t60 * t41;
t36 = t47 * t63 - t66 * t67;
t37 = t47 * t66 + t67 * t63;
t19 = pkin(2) * t36 - qJ(3) * t37 + t28;
t27 = t66 * t34 + t63 * t74;
t42 = -t46 * t57 + t60 * V_base(6) + qJD(2);
t22 = qJ(3) * t42 + t27;
t55 = sin(pkin(10));
t58 = cos(pkin(10));
t13 = t55 * t19 + t58 * t22;
t29 = -t37 * t55 + t42 * t58;
t11 = pkin(7) * t29 + t13;
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t12 = t58 * t19 - t22 * t55;
t30 = t37 * t58 + t42 * t55;
t9 = pkin(3) * t36 - pkin(7) * t30 + t12;
t6 = t65 * t11 + t62 * t9;
t5 = -t11 * t62 + t65 * t9;
t24 = t29 * t65 - t30 * t62;
t21 = -t42 * pkin(2) + qJD(3) - t26;
t14 = -t29 * pkin(3) + t21;
t64 = cos(qJ(5));
t61 = sin(qJ(5));
t35 = qJD(4) + t36;
t25 = t29 * t62 + t30 * t65;
t23 = qJD(5) - t24;
t16 = t25 * t64 + t35 * t61;
t15 = -t25 * t61 + t35 * t64;
t7 = -t24 * pkin(4) - t25 * pkin(8) + t14;
t4 = pkin(8) * t35 + t6;
t3 = -pkin(4) * t35 - t5;
t2 = t4 * t64 + t61 * t7;
t1 = -t4 * t61 + t64 * t7;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t54 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t54 * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,1) * t47 / 0.2e1) * t47 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t42 / 0.2e1) * t42 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (t21 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t54 * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t47 + Ifges(2,2) * t46 / 0.2e1) * t46 + (t28 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t42 + Ifges(3,1) * t37 / 0.2e1) * t37 + (-t21 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (t14 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t35 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t16 / 0.2e1) * t16 + (-t14 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t35 + Ifges(5,2) * t24 / 0.2e1) * t24 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t23 + Ifges(6,2) * t15 / 0.2e1) * t15 + (V_base(2) * mrSges(1,1) + t43 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t44 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t47 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t46 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t28 * mrSges(3,1) + t12 * mrSges(4,1) - t13 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t37 + Ifges(4,5) * t30 - Ifges(3,6) * t42 + Ifges(4,6) * t29 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t36) * t36;
T = t8;
