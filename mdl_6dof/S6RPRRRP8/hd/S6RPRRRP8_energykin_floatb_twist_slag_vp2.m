% Calculate kinetic energy for
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:55
% EndTime: 2019-03-09 06:22:56
% DurationCPUTime: 0.88s
% Computational Cost: add. (1911->149), mult. (2385->196), div. (0->0), fcn. (1728->8), ass. (0->50)
t63 = pkin(1) + pkin(7);
t47 = V_base(5) * pkin(6) + V_base(1);
t48 = -V_base(4) * pkin(6) + V_base(2);
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t38 = t58 * t47 + t55 * t48;
t51 = V_base(6) + qJD(1);
t34 = -t51 * qJ(2) - t38;
t42 = t55 * V_base(4) - t58 * V_base(5);
t31 = -pkin(2) * t42 - t34;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t35 = t42 * t57 - t51 * t54;
t22 = -pkin(3) * t35 + t31;
t36 = t42 * t54 + t51 * t57;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t24 = t35 * t56 - t36 * t53;
t25 = t35 * t53 + t36 * t56;
t12 = -pkin(4) * t24 - pkin(9) * t25 + t22;
t52 = sin(qJ(5));
t62 = cos(qJ(5));
t43 = t55 * V_base(5) + t58 * V_base(4);
t37 = -t55 * t47 + t48 * t58;
t60 = qJD(2) - t37;
t28 = pkin(2) * t43 - t51 * t63 + t60;
t61 = -qJ(2) * t43 + V_base(3);
t30 = t42 * t63 + t61;
t18 = t57 * t28 - t30 * t54;
t41 = qJD(3) + t43;
t14 = pkin(3) * t41 - pkin(8) * t36 + t18;
t19 = t54 * t28 + t57 * t30;
t17 = pkin(8) * t35 + t19;
t10 = t53 * t14 + t56 * t17;
t40 = qJD(4) + t41;
t8 = pkin(9) * t40 + t10;
t4 = t52 * t12 + t62 * t8;
t9 = t14 * t56 - t53 * t17;
t7 = -pkin(4) * t40 - t9;
t3 = t12 * t62 - t52 * t8;
t59 = V_base(3) ^ 2;
t33 = -pkin(1) * t51 + t60;
t32 = pkin(1) * t42 + t61;
t23 = qJD(5) - t24;
t21 = t25 * t62 + t52 * t40;
t20 = t25 * t52 - t40 * t62;
t5 = pkin(5) * t20 - qJ(6) * t21 + t7;
t2 = qJ(6) * t23 + t4;
t1 = -t23 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t37 ^ 2 + t38 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t31 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t22 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t18 * mrSges(4,1) - t19 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t40 / 0.2e1) * t40 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t31 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t36 / 0.2e1) * t36 + (t22 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t40 + Ifges(5,1) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t31 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t36 + Ifges(4,6) * t41 + Ifges(4,2) * t35 / 0.2e1) * t35 + (-t22 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t40 + Ifges(5,2) * t24 / 0.2e1) * t24 + (t37 * mrSges(2,1) - t38 * mrSges(2,2) + t33 * mrSges(3,2) - t34 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t51) * t51 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t23) * t23 + (t33 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t37 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t43 + (-Ifges(3,4) + Ifges(2,5)) * t51) * t43 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t21 + (Ifges(7,4) + Ifges(6,5)) * t23) * t21 + (V_base(3) * mrSges(2,1) + t34 * mrSges(3,1) - t32 * mrSges(3,2) - t38 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t42 + (Ifges(3,5) - Ifges(2,6)) * t51 + (-Ifges(2,4) - Ifges(3,6)) * t43) * t42 + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t20 + (-Ifges(6,6) + Ifges(7,6)) * t23 + (-Ifges(6,4) + Ifges(7,5)) * t21) * t20;
T  = t6;
