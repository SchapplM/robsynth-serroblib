% Calculate kinetic energy for
% S6RPRPRP9
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:28
% EndTime: 2019-03-09 03:27:29
% DurationCPUTime: 0.85s
% Computational Cost: add. (1887->149), mult. (2385->194), div. (0->0), fcn. (1728->8), ass. (0->49)
t63 = pkin(1) + pkin(7);
t54 = sin(qJ(5));
t61 = cos(qJ(5));
t56 = sin(qJ(1));
t62 = cos(qJ(1));
t43 = t56 * V_base(5) + t62 * V_base(4);
t51 = V_base(6) + qJD(1);
t47 = V_base(5) * pkin(6) + V_base(1);
t48 = -V_base(4) * pkin(6) + V_base(2);
t38 = -t56 * t47 + t62 * t48;
t59 = qJD(2) - t38;
t25 = t43 * pkin(2) - t63 * t51 + t59;
t42 = t56 * V_base(4) - t62 * V_base(5);
t60 = -qJ(2) * t43 + V_base(3);
t30 = t63 * t42 + t60;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t17 = t55 * t25 + t57 * t30;
t41 = qJD(3) + t43;
t15 = qJ(4) * t41 + t17;
t39 = t62 * t47 + t56 * t48;
t34 = -t51 * qJ(2) - t39;
t31 = -pkin(2) * t42 - t34;
t36 = -t57 * t42 + t51 * t55;
t37 = t42 * t55 + t51 * t57;
t22 = pkin(3) * t36 - qJ(4) * t37 + t31;
t52 = sin(pkin(9));
t53 = cos(pkin(9));
t10 = -t15 * t52 + t53 * t22;
t29 = t37 * t53 + t41 * t52;
t7 = pkin(4) * t36 - pkin(8) * t29 + t10;
t11 = t53 * t15 + t52 * t22;
t28 = -t37 * t52 + t41 * t53;
t9 = pkin(8) * t28 + t11;
t4 = t54 * t7 + t61 * t9;
t16 = t25 * t57 - t55 * t30;
t3 = -t54 * t9 + t61 * t7;
t14 = -pkin(3) * t41 + qJD(4) - t16;
t12 = -pkin(4) * t28 + t14;
t58 = V_base(3) ^ 2;
t35 = qJD(5) + t36;
t33 = -t51 * pkin(1) + t59;
t32 = pkin(1) * t42 + t60;
t19 = t54 * t28 + t61 * t29;
t18 = -t61 * t28 + t29 * t54;
t5 = pkin(5) * t18 - qJ(6) * t19 + t12;
t2 = qJ(6) * t35 + t4;
t1 = -t35 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t38 ^ 2 + t39 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t58) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t31 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (t14 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t31 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t37 / 0.2e1) * t37 + (-t14 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,2) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + t33 * mrSges(3,2) - t34 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t51) * t51 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (t33 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t43 + (-Ifges(3,4) + Ifges(2,5)) * t51) * t43 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t19 + (Ifges(7,4) + Ifges(6,5)) * t35) * t19 + (V_base(3) * mrSges(2,1) + t34 * mrSges(3,1) - t32 * mrSges(3,2) - t39 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t42 + (Ifges(3,5) - Ifges(2,6)) * t51 + (-Ifges(2,4) - Ifges(3,6)) * t43) * t42 + (t31 * mrSges(4,1) + t10 * mrSges(5,1) - t11 * mrSges(5,2) - t17 * mrSges(4,3) - Ifges(4,4) * t37 + Ifges(5,5) * t29 - Ifges(4,6) * t41 + Ifges(5,6) * t28 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t36) * t36 + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t18 + (-Ifges(6,6) + Ifges(7,6)) * t35 + (-Ifges(6,4) + Ifges(7,5)) * t19) * t18;
T  = t6;
