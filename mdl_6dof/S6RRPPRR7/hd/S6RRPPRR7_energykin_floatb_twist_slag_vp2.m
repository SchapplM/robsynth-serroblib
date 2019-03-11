% Calculate kinetic energy for
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:03
% EndTime: 2019-03-09 09:17:04
% DurationCPUTime: 1.14s
% Computational Cost: add. (2495->158), mult. (3707->209), div. (0->0), fcn. (2974->10), ass. (0->59)
t76 = -pkin(3) - pkin(9);
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t48 = t60 * V_base(5) + t63 * V_base(4);
t75 = pkin(8) * t48;
t47 = -t60 * V_base(4) + t63 * V_base(5);
t56 = sin(pkin(6));
t71 = V_base(6) + qJD(1);
t72 = cos(pkin(6));
t43 = t47 * t56 - t71 * t72 - qJD(2);
t59 = sin(qJ(2));
t68 = t56 * t71;
t67 = t47 * t72 + t68;
t74 = cos(qJ(2));
t37 = t48 * t74 + t59 * t67;
t54 = V_base(5) * pkin(7) + V_base(1);
t55 = -V_base(4) * pkin(7) + V_base(2);
t45 = t63 * t54 + t60 * t55;
t35 = pkin(8) * t67 + t45;
t44 = -t60 * t54 + t63 * t55;
t38 = pkin(1) * t71 - t72 * t75 + t44;
t70 = t72 * t74;
t41 = -pkin(1) * t47 - t56 * t75 + V_base(3);
t73 = t41 * t56;
t19 = -t59 * t35 + t38 * t70 + t73 * t74;
t66 = qJD(3) - t19;
t65 = -t37 * qJ(4) + t66;
t11 = (pkin(2) - t76) * t43 + t65;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t36 = -t47 * t70 + t48 * t59 - t68 * t74;
t23 = -t56 * t38 + t72 * t41;
t16 = t36 * pkin(2) - t37 * qJ(3) + t23;
t69 = qJD(4) - t16;
t9 = pkin(4) * t37 + t36 * t76 + t69;
t6 = t62 * t11 + t58 * t9;
t20 = t74 * t35 + (t38 * t72 + t73) * t59;
t18 = -t43 * qJ(3) + t20;
t15 = -t36 * qJ(4) - t18;
t5 = -t11 * t58 + t62 * t9;
t25 = -t36 * t58 + t43 * t62;
t13 = -pkin(4) * t43 - t15;
t64 = V_base(3) ^ 2;
t61 = cos(qJ(6));
t57 = sin(qJ(6));
t34 = qJD(5) + t37;
t26 = t36 * t62 + t43 * t58;
t24 = qJD(6) - t25;
t22 = t26 * t61 + t34 * t57;
t21 = -t26 * t57 + t34 * t61;
t17 = t43 * pkin(2) + t66;
t14 = -pkin(3) * t36 + t69;
t12 = (pkin(2) + pkin(3)) * t43 + t65;
t7 = -pkin(5) * t25 - pkin(10) * t26 + t13;
t4 = pkin(10) * t34 + t6;
t3 = -pkin(5) * t34 - t5;
t2 = t4 * t61 + t57 * t7;
t1 = -t4 * t57 + t61 * t7;
t8 = (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (-t13 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t34 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t45 * mrSges(2,3) - V_base(3) * mrSges(2,1) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t22 / 0.2e1) * t22 + m(2) * (t44 ^ 2 + t45 ^ 2 + t64) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t64) / 0.2e1 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t22 + Ifges(7,6) * t24 + Ifges(7,2) * t21 / 0.2e1) * t21 + m(3) * (t19 ^ 2 + t20 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + (t23 * mrSges(3,1) + t16 * mrSges(4,1) - t18 * mrSges(4,2) + t14 * mrSges(5,2) - t20 * mrSges(3,3) - t15 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t36 + (Ifges(5,5) + Ifges(3,6) - Ifges(4,6)) * t43 + (-Ifges(3,4) - Ifges(5,4) + Ifges(4,5)) * t37) * t36 + (t14 * mrSges(5,1) + t23 * mrSges(3,2) + t17 * mrSges(4,2) - t19 * mrSges(3,3) - t16 * mrSges(4,3) - t12 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t37 + (-Ifges(4,4) - Ifges(3,5) - Ifges(5,6)) * t43) * t37 + (t13 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t19 * mrSges(3,1) + t17 * mrSges(4,1) + t15 * mrSges(5,1) + t20 * mrSges(3,2) - t12 * mrSges(5,2) - t18 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t43) * t43 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,5) * t48 + Ifges(2,6) * t47 + Ifges(2,3) * t71 / 0.2e1) * t71;
T  = t8;
