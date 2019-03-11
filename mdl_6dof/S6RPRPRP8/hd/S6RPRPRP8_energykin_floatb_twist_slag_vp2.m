% Calculate kinetic energy for
% S6RPRPRP8
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:24
% EndTime: 2019-03-09 03:24:25
% DurationCPUTime: 0.81s
% Computational Cost: add. (1875->149), mult. (2385->194), div. (0->0), fcn. (1728->8), ass. (0->49)
t62 = pkin(1) + pkin(7);
t46 = V_base(5) * pkin(6) + V_base(1);
t47 = -V_base(4) * pkin(6) + V_base(2);
t55 = sin(qJ(1));
t61 = cos(qJ(1));
t38 = t61 * t46 + t55 * t47;
t50 = V_base(6) + qJD(1);
t34 = -t50 * qJ(2) - t38;
t41 = t55 * V_base(4) - t61 * V_base(5);
t31 = -pkin(2) * t41 - t34;
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t35 = t41 * t56 - t50 * t54;
t22 = -pkin(3) * t35 + qJD(4) + t31;
t36 = t41 * t54 + t50 * t56;
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t24 = t35 * t52 - t36 * t51;
t25 = t35 * t51 + t36 * t52;
t12 = -pkin(4) * t24 - pkin(8) * t25 + t22;
t53 = sin(qJ(5));
t60 = cos(qJ(5));
t42 = t55 * V_base(5) + t61 * V_base(4);
t37 = -t55 * t46 + t61 * t47;
t58 = qJD(2) - t37;
t28 = t42 * pkin(2) - t62 * t50 + t58;
t59 = -qJ(2) * t42 + V_base(3);
t30 = t62 * t41 + t59;
t18 = t56 * t28 - t30 * t54;
t40 = qJD(3) + t42;
t14 = pkin(3) * t40 - qJ(4) * t36 + t18;
t19 = t54 * t28 + t56 * t30;
t17 = qJ(4) * t35 + t19;
t10 = t51 * t14 + t52 * t17;
t8 = pkin(8) * t40 + t10;
t4 = t53 * t12 + t60 * t8;
t9 = t14 * t52 - t51 * t17;
t7 = -pkin(4) * t40 - t9;
t3 = t60 * t12 - t53 * t8;
t57 = V_base(3) ^ 2;
t33 = -t50 * pkin(1) + t58;
t32 = pkin(1) * t41 + t59;
t23 = qJD(5) - t24;
t21 = t60 * t25 + t53 * t40;
t20 = t25 * t53 - t60 * t40;
t5 = pkin(5) * t20 - qJ(6) * t21 + t7;
t2 = qJ(6) * t23 + t4;
t1 = -t23 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t37 ^ 2 + t38 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t31 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t22 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t31 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t36 / 0.2e1) * t36 + (t22 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,1) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t31 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t36 + Ifges(4,2) * t35 / 0.2e1) * t35 + (-t22 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,2) * t24 / 0.2e1) * t24 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t37 * mrSges(2,1) - t38 * mrSges(2,2) + t33 * mrSges(3,2) - t34 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t50) * t50 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t23) * t23 + (t33 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t37 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t42 + (-Ifges(3,4) + Ifges(2,5)) * t50) * t42 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t21 + (Ifges(7,4) + Ifges(6,5)) * t23) * t21 + (V_base(3) * mrSges(2,1) + t34 * mrSges(3,1) - t32 * mrSges(3,2) - t38 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t41 + (Ifges(3,5) - Ifges(2,6)) * t50 + (-Ifges(2,4) - Ifges(3,6)) * t42) * t41 + (t18 * mrSges(4,1) + t9 * mrSges(5,1) - t19 * mrSges(4,2) - t10 * mrSges(5,2) + Ifges(4,5) * t36 + Ifges(5,5) * t25 + Ifges(4,6) * t35 + Ifges(5,6) * t24 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t40) * t40 + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t20 + (-Ifges(6,6) + Ifges(7,6)) * t23 + (-Ifges(6,4) + Ifges(7,5)) * t21) * t20;
T  = t6;
