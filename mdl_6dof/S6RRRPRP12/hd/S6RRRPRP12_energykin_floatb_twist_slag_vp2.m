% Calculate kinetic energy for
% S6RRRPRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:31
% EndTime: 2019-03-09 17:51:32
% DurationCPUTime: 1.06s
% Computational Cost: add. (3243->154), mult. (4827->206), div. (0->0), fcn. (3980->10), ass. (0->56)
t52 = V_base(5) * pkin(7) + V_base(1);
t53 = -V_base(4) * pkin(7) + V_base(2);
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t45 = -t52 * t60 + t62 * t53;
t54 = V_base(6) + qJD(1);
t56 = cos(pkin(6));
t48 = t60 * V_base(5) + t62 * V_base(4);
t73 = pkin(8) * t48;
t39 = pkin(1) * t54 - t56 * t73 + t45;
t47 = -t60 * V_base(4) + t62 * V_base(5);
t55 = sin(pkin(6));
t42 = -pkin(1) * t47 - t55 * t73 + V_base(3);
t75 = t39 * t56 + t42 * t55;
t46 = t62 * t52 + t60 * t53;
t66 = t47 * t56 + t54 * t55;
t36 = t66 * pkin(8) + t46;
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t22 = -t59 * t36 + t75 * t61;
t37 = -t48 * t59 + t66 * t61;
t74 = pkin(3) + pkin(10);
t38 = t48 * t61 + t66 * t59;
t44 = -t47 * t55 + t54 * t56 + qJD(2);
t58 = sin(qJ(3));
t72 = cos(qJ(3));
t29 = t38 * t58 - t72 * t44;
t20 = -pkin(2) * t44 - t22;
t30 = t72 * t38 + t58 * t44;
t64 = -qJ(4) * t30 + t20;
t10 = t74 * t29 + t64;
t57 = sin(qJ(5));
t35 = qJD(3) - t37;
t26 = -t39 * t55 + t56 * t42;
t17 = -pkin(2) * t37 - pkin(9) * t38 + t26;
t23 = t61 * t36 + t75 * t59;
t21 = pkin(9) * t44 + t23;
t13 = t72 * t17 - t58 * t21;
t65 = qJD(4) - t13;
t7 = t30 * pkin(4) - t74 * t35 + t65;
t71 = cos(qJ(5));
t4 = t71 * t10 + t57 * t7;
t14 = t58 * t17 + t72 * t21;
t12 = -t35 * qJ(4) - t14;
t8 = -pkin(4) * t29 - t12;
t3 = -t57 * t10 + t71 * t7;
t63 = V_base(3) ^ 2;
t28 = qJD(5) + t30;
t25 = t57 * t29 + t71 * t35;
t24 = -t71 * t29 + t35 * t57;
t15 = pkin(3) * t29 + t64;
t11 = -t35 * pkin(3) + t65;
t5 = pkin(5) * t24 - qJ(6) * t25 + t8;
t2 = qJ(6) * t28 + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t45 ^ 2 + t46 ^ 2 + t63) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t63) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t15 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t45 * mrSges(2,1) - t46 * mrSges(2,2) + Ifges(2,3) * t54 / 0.2e1) * t54 + (t22 * mrSges(3,1) - t23 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,5) * t54 + Ifges(2,1) * t48 / 0.2e1) * t48 + (t26 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,5) * t44 + Ifges(3,1) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t46 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t54 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t26 * mrSges(3,1) + t23 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t44 + Ifges(3,2) * t37 / 0.2e1) * t37 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + t11 * mrSges(5,2) - t12 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t35) * t35 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t28) * t28 + (t11 * mrSges(5,1) + t20 * mrSges(4,2) - t13 * mrSges(4,3) - t15 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t30 + (-Ifges(5,4) + Ifges(4,5)) * t35) * t30 + (t8 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t25 + (Ifges(7,4) + Ifges(6,5)) * t28) * t25 + (t20 * mrSges(4,1) + t12 * mrSges(5,1) - t15 * mrSges(5,2) - t14 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29 + (Ifges(5,5) - Ifges(4,6)) * t35 + (-Ifges(4,4) - Ifges(5,6)) * t30) * t29 + (t8 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t24 + (-Ifges(6,6) + Ifges(7,6)) * t28 + (-Ifges(6,4) + Ifges(7,5)) * t25) * t24;
T  = t6;
