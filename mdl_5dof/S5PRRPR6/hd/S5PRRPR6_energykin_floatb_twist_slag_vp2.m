% Calculate kinetic energy for
% S5PRRPR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:29:59
% EndTime: 2019-12-05 16:30:00
% DurationCPUTime: 1.07s
% Computational Cost: add. (2803->137), mult. (4814->202), div. (0->0), fcn. (3992->12), ass. (0->55)
t52 = V_base(5) * qJ(1) + V_base(1);
t53 = -V_base(4) * qJ(1) + V_base(2);
t56 = sin(pkin(9));
t59 = cos(pkin(9));
t45 = -t52 * t56 + t59 * t53;
t60 = cos(pkin(5));
t48 = t56 * V_base(5) + t59 * V_base(4);
t72 = pkin(6) * t48;
t39 = V_base(6) * pkin(1) - t60 * t72 + t45;
t47 = -t56 * V_base(4) + t59 * V_base(5);
t54 = V_base(3) + qJD(1);
t57 = sin(pkin(5));
t42 = -pkin(1) * t47 - t57 * t72 + t54;
t73 = t39 * t60 + t42 * t57;
t46 = t59 * t52 + t56 * t53;
t66 = t47 * t60 + t57 * V_base(6);
t35 = t66 * pkin(6) + t46;
t63 = sin(qJ(2));
t65 = cos(qJ(2));
t24 = -t63 * t35 + t73 * t65;
t37 = -t63 * t48 + t66 * t65;
t28 = -t39 * t57 + t60 * t42;
t38 = t48 * t65 + t66 * t63;
t19 = -pkin(2) * t37 - pkin(7) * t38 + t28;
t25 = t65 * t35 + t73 * t63;
t44 = -t47 * t57 + t60 * V_base(6) + qJD(2);
t23 = pkin(7) * t44 + t25;
t62 = sin(qJ(3));
t71 = cos(qJ(3));
t12 = t62 * t19 + t71 * t23;
t36 = qJD(3) - t37;
t10 = qJ(4) * t36 + t12;
t22 = -t44 * pkin(2) - t24;
t30 = t38 * t62 - t71 * t44;
t31 = t71 * t38 + t62 * t44;
t15 = t30 * pkin(3) - t31 * qJ(4) + t22;
t55 = sin(pkin(10));
t58 = cos(pkin(10));
t6 = t58 * t10 + t55 * t15;
t5 = -t10 * t55 + t58 * t15;
t11 = t71 * t19 - t62 * t23;
t9 = -t36 * pkin(3) + qJD(4) - t11;
t64 = cos(qJ(5));
t61 = sin(qJ(5));
t29 = qJD(5) + t30;
t27 = t31 * t58 + t36 * t55;
t26 = -t31 * t55 + t36 * t58;
t17 = t26 * t61 + t27 * t64;
t16 = t26 * t64 - t27 * t61;
t7 = -t26 * pkin(4) + t9;
t4 = pkin(8) * t26 + t6;
t3 = pkin(4) * t30 - pkin(8) * t27 + t5;
t2 = t3 * t61 + t4 * t64;
t1 = t3 * t64 - t4 * t61;
t8 = m(2) * (t45 ^ 2 + t46 ^ 2 + t54 ^ 2) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t22 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t54 * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t36 / 0.2e1) * t36 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t54 * mrSges(2,1) + t46 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (t28 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t44 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t22 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t36 + Ifges(4,1) * t31 / 0.2e1) * t31 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t27 + Ifges(5,2) * t26 / 0.2e1) * t26 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + t29 * Ifges(6,5) + Ifges(6,1) * t17 / 0.2e1) * t17 + (-t28 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t44 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t29 + Ifges(6,2) * t16 / 0.2e1) * t16 + (V_base(2) * mrSges(1,1) + t45 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t46 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t48 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t47 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t22 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t12 * mrSges(4,3) - Ifges(4,4) * t31 + Ifges(5,5) * t27 - Ifges(4,6) * t36 + Ifges(5,6) * t26 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t30) * t30;
T = t8;
