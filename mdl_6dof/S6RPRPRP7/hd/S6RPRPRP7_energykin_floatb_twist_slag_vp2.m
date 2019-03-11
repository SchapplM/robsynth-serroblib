% Calculate kinetic energy for
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:07
% EndTime: 2019-03-09 03:21:07
% DurationCPUTime: 0.85s
% Computational Cost: add. (1879->149), mult. (2393->194), div. (0->0), fcn. (1736->8), ass. (0->49)
t62 = pkin(1) + pkin(7);
t46 = V_base(5) * pkin(6) + V_base(1);
t47 = -V_base(4) * pkin(6) + V_base(2);
t55 = sin(qJ(1));
t61 = cos(qJ(1));
t39 = t61 * t46 + t55 * t47;
t50 = V_base(6) + qJD(1);
t35 = -t50 * qJ(2) - t39;
t41 = t55 * V_base(4) - t61 * V_base(5);
t32 = -pkin(2) * t41 - t35;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t36 = t41 * t57 - t50 * t54;
t23 = -pkin(3) * t36 + qJD(4) + t32;
t37 = t41 * t54 + t50 * t57;
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t25 = t36 * t52 - t37 * t51;
t26 = t36 * t51 + t37 * t52;
t13 = -pkin(4) * t25 - pkin(8) * t26 + t23;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t42 = t55 * V_base(5) + t61 * V_base(4);
t38 = -t55 * t46 + t61 * t47;
t59 = qJD(2) - t38;
t29 = t42 * pkin(2) - t62 * t50 + t59;
t60 = -qJ(2) * t42 + V_base(3);
t31 = t62 * t41 + t60;
t19 = t57 * t29 - t31 * t54;
t40 = qJD(3) + t42;
t15 = pkin(3) * t40 - qJ(4) * t37 + t19;
t20 = t54 * t29 + t57 * t31;
t18 = qJ(4) * t36 + t20;
t10 = t51 * t15 + t52 * t18;
t8 = pkin(8) * t40 + t10;
t4 = t53 * t13 + t56 * t8;
t3 = t56 * t13 - t53 * t8;
t9 = t15 * t52 - t51 * t18;
t7 = -pkin(4) * t40 - t9;
t58 = V_base(3) ^ 2;
t34 = -t50 * pkin(1) + t59;
t33 = pkin(1) * t41 + t60;
t24 = qJD(5) - t25;
t22 = t26 * t56 + t40 * t53;
t21 = -t26 * t53 + t40 * t56;
t5 = -pkin(5) * t21 + qJD(6) + t7;
t2 = qJ(6) * t21 + t4;
t1 = pkin(5) * t24 - qJ(6) * t22 + t3;
t6 = m(2) * (t38 ^ 2 + t39 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t58) / 0.2e1 + m(4) * (t19 ^ 2 + t20 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t23 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t32 * mrSges(4,2) - t19 * mrSges(4,3) + Ifges(4,1) * t37 / 0.2e1) * t37 + (t23 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,1) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t32 * mrSges(4,1) + t20 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,2) * t36 / 0.2e1) * t36 + (-t23 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t26 + Ifges(5,2) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + t34 * mrSges(3,2) - t35 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t50) * t50 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t24) * t24 + (t34 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) - t33 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t42 + (-Ifges(3,4) + Ifges(2,5)) * t50) * t42 + (t7 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t22 + (Ifges(6,5) + Ifges(7,5)) * t24) * t22 + (V_base(3) * mrSges(2,1) + t35 * mrSges(3,1) - t33 * mrSges(3,2) - t39 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t41 + (Ifges(3,5) - Ifges(2,6)) * t50 + (-Ifges(2,4) - Ifges(3,6)) * t42) * t41 + (t19 * mrSges(4,1) + t9 * mrSges(5,1) - t20 * mrSges(4,2) - t10 * mrSges(5,2) + Ifges(4,5) * t37 + Ifges(5,5) * t26 + Ifges(4,6) * t36 + Ifges(5,6) * t25 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t40) * t40 + (-t7 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t21 + (Ifges(6,6) + Ifges(7,6)) * t24 + (Ifges(6,4) + Ifges(7,4)) * t22) * t21;
T  = t6;
