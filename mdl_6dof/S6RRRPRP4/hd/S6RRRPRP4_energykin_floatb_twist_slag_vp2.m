% Calculate kinetic energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:28
% EndTime: 2019-03-09 16:43:29
% DurationCPUTime: 0.93s
% Computational Cost: add. (2085->149), mult. (2657->196), div. (0->0), fcn. (2040->8), ass. (0->50)
t64 = pkin(3) + pkin(9);
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t45 = t56 * V_base(5) + t58 * V_base(4);
t52 = V_base(6) + qJD(1);
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t37 = -t45 * t55 + t52 * t57;
t38 = t45 * t57 + t52 * t55;
t54 = sin(qJ(3));
t63 = cos(qJ(3));
t27 = -t37 * t63 + t38 * t54;
t28 = t54 * t37 + t38 * t63;
t49 = V_base(5) * pkin(6) + V_base(1);
t50 = -V_base(4) * pkin(6) + V_base(2);
t39 = -t56 * t49 + t50 * t58;
t35 = -pkin(1) * t52 - t39;
t29 = -pkin(2) * t37 + t35;
t60 = -qJ(4) * t28 + t29;
t10 = t27 * t64 + t60;
t53 = sin(qJ(5));
t62 = cos(qJ(5));
t44 = -t56 * V_base(4) + t58 * V_base(5);
t43 = qJD(2) - t44;
t42 = qJD(3) + t43;
t32 = -pkin(1) * t44 - pkin(7) * t45 + V_base(3);
t40 = t58 * t49 + t56 * t50;
t36 = pkin(7) * t52 + t40;
t23 = t57 * t32 - t36 * t55;
t17 = pkin(2) * t43 - pkin(8) * t38 + t23;
t24 = t55 * t32 + t57 * t36;
t20 = pkin(8) * t37 + t24;
t13 = t17 * t63 - t54 * t20;
t61 = qJD(4) - t13;
t7 = t28 * pkin(4) - t42 * t64 + t61;
t4 = t62 * t10 + t53 * t7;
t14 = t54 * t17 + t63 * t20;
t12 = -t42 * qJ(4) - t14;
t8 = -pkin(4) * t27 - t12;
t3 = -t53 * t10 + t62 * t7;
t59 = V_base(3) ^ 2;
t26 = qJD(5) + t28;
t22 = t53 * t27 + t42 * t62;
t21 = -t27 * t62 + t42 * t53;
t15 = pkin(3) * t27 + t60;
t11 = -t42 * pkin(3) + t61;
t5 = pkin(5) * t21 - qJ(6) * t22 + t8;
t2 = qJ(6) * t26 + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t39 ^ 2 + t40 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t23 ^ 2 + t24 ^ 2 + t35 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t15 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t52 / 0.2e1) * t52 + (t23 * mrSges(3,1) - t24 * mrSges(3,2) + Ifges(3,3) * t43 / 0.2e1) * t43 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,5) * t52 + Ifges(2,1) * t45 / 0.2e1) * t45 + (t35 * mrSges(3,2) - t23 * mrSges(3,3) + Ifges(3,5) * t43 + Ifges(3,1) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,4) * t45 + Ifges(2,6) * t52 + Ifges(2,2) * t44 / 0.2e1) * t44 + (-t35 * mrSges(3,1) + t24 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t43 + Ifges(3,2) * t37 / 0.2e1) * t37 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + t11 * mrSges(5,2) - t12 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t42) * t42 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t26) * t26 + (t11 * mrSges(5,1) + t29 * mrSges(4,2) - t13 * mrSges(4,3) - t15 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t28 + (-Ifges(5,4) + Ifges(4,5)) * t42) * t28 + (t8 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5)) * t26) * t22 + (t29 * mrSges(4,1) + t12 * mrSges(5,1) - t15 * mrSges(5,2) - t14 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t27 + (Ifges(5,5) - Ifges(4,6)) * t42 + (-Ifges(4,4) - Ifges(5,6)) * t28) * t27 + (t8 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,6) + Ifges(7,6)) * t26 + (-Ifges(6,4) + Ifges(7,5)) * t22) * t21;
T  = t6;
