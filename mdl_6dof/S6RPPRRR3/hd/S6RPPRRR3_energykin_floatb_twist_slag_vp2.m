% Calculate kinetic energy for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:54
% EndTime: 2019-03-09 02:22:54
% DurationCPUTime: 0.96s
% Computational Cost: add. (2319->153), mult. (3351->211), div. (0->0), fcn. (2612->10), ass. (0->56)
t70 = pkin(2) + pkin(7);
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t47 = -t61 * V_base(4) + t65 * V_base(5);
t48 = t61 * V_base(5) + t65 * V_base(4);
t57 = sin(pkin(10));
t69 = cos(pkin(10));
t42 = t57 * t47 + t69 * t48;
t56 = V_base(6) + qJD(1);
t53 = V_base(5) * pkin(6) + V_base(1);
t54 = -V_base(4) * pkin(6) + V_base(2);
t43 = -t53 * t61 + t65 * t54;
t36 = pkin(1) * t56 - qJ(2) * t48 + t43;
t44 = t65 * t53 + t61 * t54;
t39 = qJ(2) * t47 + t44;
t29 = t69 * t36 - t57 * t39;
t68 = qJD(3) - t29;
t19 = t42 * pkin(3) - t70 * t56 + t68;
t41 = -t69 * t47 + t48 * t57;
t45 = -pkin(1) * t47 + qJD(2) + V_base(3);
t67 = -qJ(3) * t42 + t45;
t23 = t70 * t41 + t67;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t12 = t60 * t19 + t64 * t23;
t40 = qJD(4) + t42;
t10 = pkin(8) * t40 + t12;
t30 = t57 * t36 + t69 * t39;
t25 = -t56 * qJ(3) - t30;
t20 = -pkin(3) * t41 - t25;
t34 = t41 * t64 - t60 * t56;
t35 = t41 * t60 + t56 * t64;
t15 = -pkin(4) * t34 - pkin(8) * t35 + t20;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t6 = t63 * t10 + t59 * t15;
t5 = -t10 * t59 + t63 * t15;
t11 = t19 * t64 - t60 * t23;
t32 = qJD(5) - t34;
t9 = -pkin(4) * t40 - t11;
t66 = V_base(3) ^ 2;
t62 = cos(qJ(6));
t58 = sin(qJ(6));
t31 = qJD(6) + t32;
t28 = t35 * t63 + t40 * t59;
t27 = -t35 * t59 + t40 * t63;
t26 = pkin(2) * t41 + t67;
t24 = -t56 * pkin(2) + t68;
t17 = t27 * t58 + t28 * t62;
t16 = t27 * t62 - t28 * t58;
t7 = -pkin(5) * t27 + t9;
t4 = pkin(9) * t27 + t6;
t3 = pkin(5) * t32 - pkin(9) * t28 + t5;
t2 = t3 * t58 + t4 * t62;
t1 = t3 * t62 - t4 * t58;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t20 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t40 / 0.2e1) * t40 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t32 / 0.2e1) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (t20 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t40 + Ifges(5,1) * t35 / 0.2e1) * t35 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t32 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t31 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t20 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,6) * t40 + Ifges(5,2) * t34 / 0.2e1) * t34 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t32 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t31 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t24 * mrSges(4,1) + t45 * mrSges(3,2) - t29 * mrSges(3,3) - t26 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t42) * t42 + (t45 * mrSges(3,1) + t25 * mrSges(4,1) - t26 * mrSges(4,2) - t30 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t41 + (-Ifges(3,4) - Ifges(4,6)) * t42) * t41 + (t43 * mrSges(2,1) + t29 * mrSges(3,1) - t44 * mrSges(2,2) - t30 * mrSges(3,2) + t24 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(2,5) * t48 + Ifges(2,6) * t47 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t56 + (-Ifges(4,4) + Ifges(3,5)) * t42 + (Ifges(4,5) - Ifges(3,6)) * t41) * t56;
T  = t8;
