% Calculate kinetic energy for
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:05
% EndTime: 2019-03-08 21:40:06
% DurationCPUTime: 1.01s
% Computational Cost: add. (2843->154), mult. (4849->205), div. (0->0), fcn. (4002->10), ass. (0->55)
t52 = V_base(5) * qJ(1) + V_base(1);
t53 = -V_base(4) * qJ(1) + V_base(2);
t55 = sin(pkin(10));
t57 = cos(pkin(10));
t45 = -t52 * t55 + t57 * t53;
t58 = cos(pkin(6));
t48 = t55 * V_base(5) + t57 * V_base(4);
t72 = pkin(7) * t48;
t39 = V_base(6) * pkin(1) - t58 * t72 + t45;
t47 = -t55 * V_base(4) + t57 * V_base(5);
t54 = V_base(3) + qJD(1);
t56 = sin(pkin(6));
t42 = -pkin(1) * t47 - t56 * t72 + t54;
t74 = t39 * t58 + t42 * t56;
t46 = t57 * t52 + t55 * t53;
t65 = t47 * t58 + t56 * V_base(6);
t35 = pkin(7) * t65 + t46;
t61 = sin(qJ(2));
t63 = cos(qJ(2));
t23 = -t61 * t35 + t63 * t74;
t37 = -t61 * t48 + t63 * t65;
t73 = pkin(3) + pkin(9);
t38 = t48 * t63 + t61 * t65;
t44 = -t47 * t56 + t58 * V_base(6) + qJD(2);
t60 = sin(qJ(3));
t71 = cos(qJ(3));
t29 = t38 * t60 - t44 * t71;
t21 = -t44 * pkin(2) - t23;
t30 = t38 * t71 + t60 * t44;
t64 = -t30 * qJ(4) + t21;
t12 = t29 * t73 + t64;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t36 = qJD(3) - t37;
t27 = -t39 * t56 + t58 * t42;
t18 = -pkin(2) * t37 - pkin(8) * t38 + t27;
t24 = t63 * t35 + t74 * t61;
t22 = pkin(8) * t44 + t24;
t14 = t18 * t71 - t60 * t22;
t66 = qJD(4) - t14;
t8 = t30 * pkin(4) - t36 * t73 + t66;
t4 = t62 * t12 + t59 * t8;
t15 = t60 * t18 + t71 * t22;
t13 = -t36 * qJ(4) - t15;
t3 = -t12 * t59 + t62 * t8;
t9 = -pkin(4) * t29 - t13;
t28 = qJD(5) + t30;
t26 = t29 * t59 + t36 * t62;
t25 = t29 * t62 - t36 * t59;
t16 = t29 * pkin(3) + t64;
t11 = -t36 * pkin(3) + t66;
t5 = -pkin(5) * t25 + qJD(6) + t9;
t2 = qJ(6) * t25 + t4;
t1 = pkin(5) * t28 - qJ(6) * t26 + t3;
t6 = m(2) * (t45 ^ 2 + t46 ^ 2 + t54 ^ 2) / 0.2e1 + m(3) * (t23 ^ 2 + t24 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t54 * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t23 * mrSges(3,1) - t24 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t54 * mrSges(2,1) + t46 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (t27 * mrSges(3,2) - t23 * mrSges(3,3) + Ifges(3,5) * t44 + Ifges(3,1) * t38 / 0.2e1) * t38 + (-t27 * mrSges(3,1) + t24 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t44 + Ifges(3,2) * t37 / 0.2e1) * t37 + (t14 * mrSges(4,1) - t15 * mrSges(4,2) + t11 * mrSges(5,2) - t13 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t36) * t36 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t28) * t28 + (t11 * mrSges(5,1) + t21 * mrSges(4,2) - t14 * mrSges(4,3) - t16 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t30 + (-Ifges(5,4) + Ifges(4,5)) * t36) * t30 + (t9 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t26 + (Ifges(6,5) + Ifges(7,5)) * t28) * t26 + (V_base(2) * mrSges(1,1) + t45 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t46 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t48 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t47 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t21 * mrSges(4,1) + t13 * mrSges(5,1) - t16 * mrSges(5,2) - t15 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29 + (Ifges(5,5) - Ifges(4,6)) * t36 + (-Ifges(4,4) - Ifges(5,6)) * t30) * t29 + (-t9 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t25 + (Ifges(6,6) + Ifges(7,6)) * t28 + (Ifges(6,4) + Ifges(7,4)) * t26) * t25;
T  = t6;
