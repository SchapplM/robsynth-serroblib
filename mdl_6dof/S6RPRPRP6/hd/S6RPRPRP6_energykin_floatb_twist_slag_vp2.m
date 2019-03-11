% Calculate kinetic energy for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:17:52
% EndTime: 2019-03-09 03:17:53
% DurationCPUTime: 0.84s
% Computational Cost: add. (2013->149), mult. (2665->194), div. (0->0), fcn. (2048->8), ass. (0->49)
t63 = pkin(3) + pkin(8);
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t44 = t56 * V_base(5) + t58 * V_base(4);
t51 = V_base(6) + qJD(1);
t52 = sin(pkin(9));
t53 = cos(pkin(9));
t37 = -t44 * t52 + t51 * t53;
t38 = t44 * t53 + t51 * t52;
t55 = sin(qJ(3));
t62 = cos(qJ(3));
t27 = -t62 * t37 + t38 * t55;
t28 = t55 * t37 + t62 * t38;
t48 = V_base(5) * pkin(6) + V_base(1);
t49 = -V_base(4) * pkin(6) + V_base(2);
t39 = -t56 * t48 + t49 * t58;
t34 = -pkin(1) * t51 + qJD(2) - t39;
t29 = -pkin(2) * t37 + t34;
t60 = -qJ(4) * t28 + t29;
t11 = t63 * t27 + t60;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t43 = t56 * V_base(4) - t58 * V_base(5);
t42 = qJD(3) + t43;
t32 = pkin(1) * t43 - qJ(2) * t44 + V_base(3);
t40 = t58 * t48 + t56 * t49;
t36 = qJ(2) * t51 + t40;
t24 = t53 * t32 - t36 * t52;
t18 = pkin(2) * t43 - pkin(7) * t38 + t24;
t25 = t52 * t32 + t53 * t36;
t21 = pkin(7) * t37 + t25;
t14 = t62 * t18 - t55 * t21;
t61 = qJD(4) - t14;
t8 = t28 * pkin(4) - t63 * t42 + t61;
t4 = t57 * t11 + t54 * t8;
t15 = t55 * t18 + t62 * t21;
t13 = -t42 * qJ(4) - t15;
t3 = -t11 * t54 + t57 * t8;
t9 = -pkin(4) * t27 - t13;
t59 = V_base(3) ^ 2;
t26 = qJD(5) + t28;
t23 = t27 * t54 + t42 * t57;
t22 = t27 * t57 - t42 * t54;
t16 = pkin(3) * t27 + t60;
t12 = -t42 * pkin(3) + t61;
t5 = -pkin(5) * t22 + qJD(6) + t9;
t2 = qJ(6) * t22 + t4;
t1 = pkin(5) * t26 - qJ(6) * t23 + t3;
t6 = m(2) * (t39 ^ 2 + t40 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t34 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t51 / 0.2e1) * t51 + (t34 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,1) * t38 / 0.2e1) * t38 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t44 / 0.2e1) * t44 + (-t34 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,2) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t14 * mrSges(4,1) - t15 * mrSges(4,2) + t12 * mrSges(5,2) - t13 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t42) * t42 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t26) * t26 + (t12 * mrSges(5,1) + t29 * mrSges(4,2) - t14 * mrSges(4,3) - t16 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t28 + (-Ifges(5,4) + Ifges(4,5)) * t42) * t28 + (t9 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t23 + (Ifges(6,5) + Ifges(7,5)) * t26) * t23 + (V_base(3) * mrSges(2,1) + t24 * mrSges(3,1) - t25 * mrSges(3,2) - t40 * mrSges(2,3) - Ifges(2,4) * t44 + Ifges(3,5) * t38 - Ifges(2,6) * t51 + Ifges(3,6) * t37 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t43) * t43 + (t29 * mrSges(4,1) + t13 * mrSges(5,1) - t16 * mrSges(5,2) - t15 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t27 + (Ifges(5,5) - Ifges(4,6)) * t42 + (-Ifges(4,4) - Ifges(5,6)) * t28) * t27 + (-t9 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t22 + (Ifges(6,6) + Ifges(7,6)) * t26 + (Ifges(6,4) + Ifges(7,4)) * t23) * t22;
T  = t6;
