% Calculate kinetic energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:46
% EndTime: 2019-03-08 20:13:47
% DurationCPUTime: 0.97s
% Computational Cost: add. (2683->156), mult. (4587->209), div. (0->0), fcn. (3756->10), ass. (0->56)
t72 = pkin(2) + pkin(8);
t56 = sin(pkin(10));
t58 = cos(pkin(10));
t48 = t56 * V_base(5) + t58 * V_base(4);
t71 = pkin(7) * t48;
t53 = V_base(5) * qJ(1) + V_base(1);
t54 = -V_base(4) * qJ(1) + V_base(2);
t45 = t58 * t53 + t56 * t54;
t47 = -t56 * V_base(4) + t58 * V_base(5);
t57 = sin(pkin(6));
t59 = cos(pkin(6));
t66 = t47 * t59 + t57 * V_base(6);
t34 = pkin(7) * t66 + t45;
t44 = -t53 * t56 + t58 * t54;
t38 = V_base(6) * pkin(1) - t59 * t71 + t44;
t55 = V_base(3) + qJD(1);
t41 = -pkin(1) * t47 - t57 * t71 + t55;
t62 = sin(qJ(2));
t70 = cos(qJ(2));
t24 = t70 * t34 + (t38 * t59 + t41 * t57) * t62;
t43 = -t47 * t57 + t59 * V_base(6) + qJD(2);
t22 = -t43 * qJ(3) - t24;
t68 = t59 * t70;
t69 = t57 * t70;
t36 = -t47 * t68 + t48 * t62 - t69 * V_base(6);
t19 = -pkin(3) * t36 - t22;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t29 = t36 * t64 - t61 * t43;
t30 = t61 * t36 + t43 * t64;
t13 = -pkin(4) * t29 - pkin(9) * t30 + t19;
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t37 = t48 * t70 + t62 * t66;
t23 = -t62 * t34 + t38 * t68 + t41 * t69;
t65 = qJD(3) - t23;
t15 = t37 * pkin(3) - t43 * t72 + t65;
t27 = -t38 * t57 + t59 * t41;
t67 = -qJ(3) * t37 + t27;
t18 = t36 * t72 + t67;
t10 = t61 * t15 + t64 * t18;
t35 = qJD(4) + t37;
t8 = pkin(9) * t35 + t10;
t4 = t60 * t13 + t63 * t8;
t3 = t63 * t13 - t60 * t8;
t9 = t15 * t64 - t61 * t18;
t7 = -t35 * pkin(4) - t9;
t28 = qJD(5) - t29;
t26 = t30 * t63 + t35 * t60;
t25 = -t30 * t60 + t35 * t63;
t21 = -t43 * pkin(2) + t65;
t20 = pkin(2) * t36 + t67;
t5 = -t25 * pkin(5) + qJD(6) + t7;
t2 = qJ(6) * t25 + t4;
t1 = pkin(5) * t28 - qJ(6) * t26 + t3;
t6 = m(2) * (t44 ^ 2 + t45 ^ 2 + t55 ^ 2) / 0.2e1 + m(3) * (t23 ^ 2 + t24 ^ 2 + t27 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t19 ^ 2 + t9 ^ 2) / 0.2e1 + m(4) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t55 * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t55 * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (t19 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t35 + Ifges(5,1) * t30 / 0.2e1) * t30 + (-t19 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t30 + Ifges(5,6) * t35 + Ifges(5,2) * t29 / 0.2e1) * t29 + (t23 * mrSges(3,1) - t24 * mrSges(3,2) + t21 * mrSges(4,2) - t22 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t43) * t43 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t28) * t28 + (t21 * mrSges(4,1) + t27 * mrSges(3,2) - t23 * mrSges(3,3) - t20 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t37 + (-Ifges(4,4) + Ifges(3,5)) * t43) * t37 + (t7 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t26 + (Ifges(6,5) + Ifges(7,5)) * t28) * t26 + (V_base(2) * mrSges(1,1) + t44 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t45 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t48 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t47 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t27 * mrSges(3,1) + t22 * mrSges(4,1) - t20 * mrSges(4,2) - t24 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t36 + (Ifges(4,5) - Ifges(3,6)) * t43 + (-Ifges(3,4) - Ifges(4,6)) * t37) * t36 + (-t7 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t25 + (Ifges(6,6) + Ifges(7,6)) * t28 + (Ifges(6,4) + Ifges(7,4)) * t26) * t25;
T  = t6;
