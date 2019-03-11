% Calculate kinetic energy for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:35
% EndTime: 2019-03-09 04:38:36
% DurationCPUTime: 0.95s
% Computational Cost: add. (3111->152), mult. (4145->209), div. (0->0), fcn. (3328->10), ass. (0->52)
t57 = sin(pkin(10));
t67 = cos(pkin(10));
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t48 = t62 * V_base(4) - t65 * V_base(5);
t49 = t62 * V_base(5) + t65 * V_base(4);
t39 = pkin(1) * t48 - qJ(2) * t49 + V_base(3);
t53 = V_base(5) * pkin(6) + V_base(1);
t54 = -V_base(4) * pkin(6) + V_base(2);
t46 = t65 * t53 + t62 * t54;
t56 = V_base(6) + qJD(1);
t42 = qJ(2) * t56 + t46;
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t31 = t59 * t39 - t42 * t58;
t44 = t49 * t59 + t56 * t58;
t24 = pkin(2) * t48 - pkin(7) * t44 + t31;
t32 = t58 * t39 + t59 * t42;
t43 = -t49 * t58 + t56 * t59;
t27 = pkin(7) * t43 + t32;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t17 = t61 * t24 + t64 * t27;
t47 = qJD(3) + t48;
t15 = pkin(8) * t47 + t17;
t34 = t43 * t64 - t44 * t61;
t35 = t43 * t61 + t44 * t64;
t45 = -t62 * t53 + t54 * t65;
t41 = -pkin(1) * t56 + qJD(2) - t45;
t36 = -pkin(2) * t43 + t41;
t20 = -pkin(3) * t34 - pkin(8) * t35 + t36;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t10 = -t15 * t60 + t63 * t20;
t30 = t35 * t63 + t47 * t60;
t33 = qJD(4) - t34;
t7 = pkin(4) * t33 - qJ(5) * t30 + t10;
t11 = t63 * t15 + t60 * t20;
t29 = -t35 * t60 + t47 * t63;
t9 = qJ(5) * t29 + t11;
t4 = t57 * t7 + t67 * t9;
t16 = t24 * t64 - t61 * t27;
t14 = -pkin(3) * t47 - t16;
t3 = -t57 * t9 + t67 * t7;
t12 = -pkin(4) * t29 + qJD(5) + t14;
t66 = V_base(3) ^ 2;
t22 = t57 * t29 + t30 * t67;
t21 = -t29 * t67 + t30 * t57;
t5 = pkin(5) * t21 - qJ(6) * t22 + t12;
t2 = qJ(6) * t33 + t4;
t1 = -t33 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t45 ^ 2 + t46 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t36 ^ 2) / 0.2e1 + m(3) * (t31 ^ 2 + t32 ^ 2 + t41 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t45 * mrSges(2,1) - t46 * mrSges(2,2) + Ifges(2,3) * t56 / 0.2e1) * t56 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t47 / 0.2e1) * t47 + (t41 * mrSges(3,2) - t31 * mrSges(3,3) + Ifges(3,1) * t44 / 0.2e1) * t44 + (t14 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t30 / 0.2e1) * t30 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,5) * t56 + Ifges(2,1) * t49 / 0.2e1) * t49 + (-t41 * mrSges(3,1) + t32 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t36 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t47 + Ifges(4,1) * t35 / 0.2e1) * t35 + (-t14 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t30 + Ifges(5,2) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t36 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t35 + Ifges(4,6) * t47 + Ifges(4,2) * t34 / 0.2e1) * t34 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t22) * t22 + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21 + (V_base(3) * mrSges(2,1) + t31 * mrSges(3,1) - t32 * mrSges(3,2) - t46 * mrSges(2,3) - Ifges(2,4) * t49 + Ifges(3,5) * t44 - Ifges(2,6) * t56 + Ifges(3,6) * t43 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t48) * t48 + (t10 * mrSges(5,1) + t3 * mrSges(6,1) - t1 * mrSges(7,1) - t11 * mrSges(5,2) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + Ifges(5,5) * t30 + Ifges(5,6) * t29 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t33 + (Ifges(7,4) + Ifges(6,5)) * t22 + (-Ifges(6,6) + Ifges(7,6)) * t21) * t33;
T  = t6;
