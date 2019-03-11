% Calculate kinetic energy for
% S6RPRRRP9
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:03
% EndTime: 2019-03-09 06:26:04
% DurationCPUTime: 0.92s
% Computational Cost: add. (1927->149), mult. (2413->196), div. (0->0), fcn. (1756->8), ass. (0->50)
t64 = pkin(1) + pkin(7);
t56 = sin(qJ(1));
t63 = cos(qJ(1));
t43 = t56 * V_base(5) + t63 * V_base(4);
t52 = V_base(6) + qJD(1);
t48 = V_base(5) * pkin(6) + V_base(1);
t49 = -V_base(4) * pkin(6) + V_base(2);
t39 = -t56 * t48 + t49 * t63;
t61 = qJD(2) - t39;
t25 = t43 * pkin(2) - t52 * t64 + t61;
t42 = t56 * V_base(4) - t63 * V_base(5);
t62 = -qJ(2) * t43 + V_base(3);
t30 = t42 * t64 + t62;
t55 = sin(qJ(3));
t59 = cos(qJ(3));
t18 = t55 * t25 + t59 * t30;
t41 = qJD(3) + t43;
t16 = pkin(8) * t41 + t18;
t40 = t63 * t48 + t56 * t49;
t34 = -t52 * qJ(2) - t40;
t31 = -pkin(2) * t42 - t34;
t37 = t42 * t59 - t55 * t52;
t38 = t42 * t55 + t52 * t59;
t23 = -pkin(3) * t37 - pkin(8) * t38 + t31;
t54 = sin(qJ(4));
t58 = cos(qJ(4));
t12 = t58 * t16 + t54 * t23;
t28 = -t38 * t54 + t41 * t58;
t10 = pkin(9) * t28 + t12;
t53 = sin(qJ(5));
t57 = cos(qJ(5));
t11 = -t16 * t54 + t58 * t23;
t29 = t38 * t58 + t41 * t54;
t36 = qJD(4) - t37;
t7 = pkin(4) * t36 - pkin(9) * t29 + t11;
t4 = t57 * t10 + t53 * t7;
t3 = -t10 * t53 + t57 * t7;
t17 = t25 * t59 - t55 * t30;
t15 = -pkin(3) * t41 - t17;
t13 = -pkin(4) * t28 + t15;
t60 = V_base(3) ^ 2;
t35 = qJD(5) + t36;
t33 = -t52 * pkin(1) + t61;
t32 = pkin(1) * t42 + t62;
t20 = t28 * t53 + t29 * t57;
t19 = t28 * t57 - t29 * t53;
t9 = -pkin(5) * t19 + qJD(6) + t13;
t2 = qJ(6) * t19 + t4;
t1 = pkin(5) * t35 - qJ(6) * t20 + t3;
t5 = m(2) * (t39 ^ 2 + t40 ^ 2 + t60) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t60) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t31 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t15 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t17 * mrSges(4,1) - t18 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t36 / 0.2e1) * t36 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t31 * mrSges(4,2) - t17 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t38 / 0.2e1) * t38 + (t15 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t36 + Ifges(5,1) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t31 * mrSges(4,1) + t18 * mrSges(4,3) + Ifges(4,4) * t38 + Ifges(4,6) * t41 + Ifges(4,2) * t37 / 0.2e1) * t37 + (-t15 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,6) * t36 + Ifges(5,2) * t28 / 0.2e1) * t28 + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + t33 * mrSges(3,2) - t34 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t52) * t52 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (t33 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t43 + (-Ifges(3,4) + Ifges(2,5)) * t52) * t43 + (t13 * mrSges(6,2) + t9 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t20 + (Ifges(6,5) + Ifges(7,5)) * t35) * t20 + (V_base(3) * mrSges(2,1) + t34 * mrSges(3,1) - t32 * mrSges(3,2) - t40 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t42 + (Ifges(3,5) - Ifges(2,6)) * t52 + (-Ifges(2,4) - Ifges(3,6)) * t43) * t42 + (-t13 * mrSges(6,1) - t9 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t19 + (Ifges(6,6) + Ifges(7,6)) * t35 + (Ifges(6,4) + Ifges(7,4)) * t20) * t19;
T  = t5;
