% Calculate kinetic energy for
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:13
% EndTime: 2019-03-09 01:33:14
% DurationCPUTime: 0.93s
% Computational Cost: add. (2197->153), mult. (3035->209), div. (0->0), fcn. (2276->10), ass. (0->53)
t63 = sin(qJ(1));
t68 = cos(qJ(1));
t47 = t63 * V_base(4) - t68 * V_base(5);
t48 = t63 * V_base(5) + t68 * V_base(4);
t39 = t47 * pkin(1) - t48 * qJ(2) + V_base(3);
t30 = -pkin(2) * t47 + qJD(3) - t39;
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t37 = -t60 * t47 + t48 * t58;
t38 = t47 * t58 + t48 * t60;
t19 = pkin(3) * t37 - qJ(4) * t38 + t30;
t56 = V_base(6) + qJD(1);
t52 = V_base(5) * pkin(6) + V_base(1);
t53 = -V_base(4) * pkin(6) + V_base(2);
t42 = -t63 * t52 + t53 * t68;
t67 = qJD(2) - t42;
t29 = -t48 * qJ(3) + (-pkin(1) - pkin(2)) * t56 + t67;
t43 = t68 * t52 + t63 * t53;
t41 = t56 * qJ(2) + t43;
t35 = qJ(3) * t47 + t41;
t25 = t58 * t29 + t60 * t35;
t22 = -qJ(4) * t56 + t25;
t57 = sin(pkin(10));
t59 = cos(pkin(10));
t13 = t57 * t19 + t59 * t22;
t31 = -t38 * t57 - t56 * t59;
t11 = pkin(7) * t31 + t13;
t62 = sin(qJ(5));
t65 = cos(qJ(5));
t12 = t59 * t19 - t22 * t57;
t32 = t38 * t59 - t56 * t57;
t9 = pkin(4) * t37 - pkin(7) * t32 + t12;
t6 = t65 * t11 + t62 * t9;
t24 = t29 * t60 - t58 * t35;
t5 = -t11 * t62 + t65 * t9;
t26 = t31 * t65 - t32 * t62;
t21 = pkin(3) * t56 + qJD(4) - t24;
t14 = -pkin(4) * t31 + t21;
t66 = V_base(3) ^ 2;
t64 = cos(qJ(6));
t61 = sin(qJ(6));
t40 = -t56 * pkin(1) + t67;
t36 = qJD(5) + t37;
t27 = t31 * t62 + t32 * t65;
t23 = qJD(6) - t26;
t16 = t27 * t64 + t36 * t61;
t15 = -t27 * t61 + t36 * t64;
t7 = -pkin(5) * t26 - pkin(8) * t27 + t14;
t4 = pkin(8) * t36 + t6;
t3 = -pkin(5) * t36 - t5;
t2 = t4 * t64 + t61 * t7;
t1 = -t4 * t61 + t64 * t7;
t8 = m(2) * (t42 ^ 2 + t43 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t30 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t21 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t30 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,1) * t38 / 0.2e1) * t38 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t36 / 0.2e1) * t36 + (t21 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t32 / 0.2e1) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t21 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t32 + Ifges(5,2) * t31 / 0.2e1) * t31 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t36 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t23 + Ifges(7,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t36 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t23 + Ifges(7,2) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(2,2) + t40 * mrSges(3,2) - t42 * mrSges(2,3) - t39 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t48) * t48 + (V_base(3) * mrSges(2,1) + t39 * mrSges(3,1) - t41 * mrSges(3,2) - t43 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t47 + (-Ifges(2,4) + Ifges(3,5)) * t48) * t47 + (t30 * mrSges(4,1) + t12 * mrSges(5,1) - t13 * mrSges(5,2) - t25 * mrSges(4,3) - Ifges(4,4) * t38 + Ifges(5,5) * t32 + Ifges(5,6) * t31 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t37) * t37 + (t42 * mrSges(2,1) - t40 * mrSges(3,1) - t24 * mrSges(4,1) - t43 * mrSges(2,2) + t25 * mrSges(4,2) + t41 * mrSges(3,3) - Ifges(4,5) * t38 + Ifges(4,6) * t37 + (Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t56 + (Ifges(3,4) + Ifges(2,5)) * t48 + (-Ifges(2,6) + Ifges(3,6)) * t47) * t56;
T  = t8;
