% Calculate kinetic energy for
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:42
% EndTime: 2019-12-05 17:14:43
% DurationCPUTime: 1.14s
% Computational Cost: add. (2843->137), mult. (4866->204), div. (0->0), fcn. (4038->12), ass. (0->56)
t52 = V_base(5) * qJ(1) + V_base(1);
t53 = -V_base(4) * qJ(1) + V_base(2);
t55 = sin(pkin(10));
t57 = cos(pkin(10));
t44 = -t52 * t55 + t57 * t53;
t58 = cos(pkin(5));
t48 = t55 * V_base(5) + t57 * V_base(4);
t71 = pkin(6) * t48;
t39 = V_base(6) * pkin(1) - t58 * t71 + t44;
t47 = -t55 * V_base(4) + t57 * V_base(5);
t54 = V_base(3) + qJD(1);
t56 = sin(pkin(5));
t42 = -pkin(1) * t47 - t56 * t71 + t54;
t72 = t39 * t58 + t42 * t56;
t45 = t57 * t52 + t55 * t53;
t67 = t47 * t58 + t56 * V_base(6);
t35 = t67 * pkin(6) + t45;
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t26 = -t62 * t35 + t72 * t66;
t37 = -t62 * t48 + t67 * t66;
t28 = -t39 * t56 + t58 * t42;
t38 = t48 * t66 + t67 * t62;
t19 = -pkin(2) * t37 - pkin(7) * t38 + t28;
t27 = t66 * t35 + t72 * t62;
t43 = -t47 * t56 + t58 * V_base(6) + qJD(2);
t22 = pkin(7) * t43 + t27;
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t13 = t61 * t19 + t65 * t22;
t29 = -t38 * t61 + t43 * t65;
t11 = pkin(8) * t29 + t13;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t12 = t65 * t19 - t22 * t61;
t30 = t38 * t65 + t43 * t61;
t36 = qJD(3) - t37;
t9 = pkin(3) * t36 - pkin(8) * t30 + t12;
t6 = t64 * t11 + t60 * t9;
t5 = -t11 * t60 + t64 * t9;
t24 = t29 * t64 - t30 * t60;
t21 = -t43 * pkin(2) - t26;
t14 = -t29 * pkin(3) + t21;
t63 = cos(qJ(5));
t59 = sin(qJ(5));
t34 = qJD(4) + t36;
t25 = t29 * t60 + t30 * t64;
t23 = qJD(5) - t24;
t16 = t25 * t63 + t34 * t59;
t15 = -t25 * t59 + t34 * t63;
t7 = -t24 * pkin(4) - t25 * pkin(9) + t14;
t4 = pkin(9) * t34 + t6;
t3 = -pkin(4) * t34 - t5;
t2 = t4 * t63 + t59 * t7;
t1 = -t4 * t59 + t63 * t7;
t8 = m(2) * (t44 ^ 2 + t45 ^ 2 + t54 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t54 * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t43 / 0.2e1) * t43 + (t12 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,3) * t36 / 0.2e1) * t36 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t34 / 0.2e1) * t34 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t54 * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (t28 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t43 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t21 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,5) * t36 + Ifges(4,1) * t30 / 0.2e1) * t30 + (t14 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t34 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t16 / 0.2e1) * t16 + (-t28 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t43 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t21 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,6) * t36 + Ifges(4,2) * t29 / 0.2e1) * t29 + (-t14 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t34 + Ifges(5,2) * t24 / 0.2e1) * t24 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + t16 * Ifges(6,4) + t23 * Ifges(6,6) + Ifges(6,2) * t15 / 0.2e1) * t15 + (V_base(2) * mrSges(1,1) + t44 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t45 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t48 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t47 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t8;
