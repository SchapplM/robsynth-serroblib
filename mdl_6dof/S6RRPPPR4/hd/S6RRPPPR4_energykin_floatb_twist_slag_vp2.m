% Calculate kinetic energy for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:10
% EndTime: 2019-03-09 08:17:11
% DurationCPUTime: 0.98s
% Computational Cost: add. (1673->150), mult. (2111->194), div. (0->0), fcn. (1540->8), ass. (0->52)
t65 = -pkin(4) - pkin(5);
t64 = pkin(2) + qJ(4);
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t42 = t53 * V_base(5) + t56 * V_base(4);
t49 = V_base(6) + qJD(1);
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t36 = t42 * t55 + t49 * t52;
t41 = -t53 * V_base(4) + t56 * V_base(5);
t40 = qJD(2) - t41;
t27 = -pkin(1) * t41 - pkin(7) * t42 + V_base(3);
t47 = V_base(5) * pkin(6) + V_base(1);
t48 = -V_base(4) * pkin(6) + V_base(2);
t38 = t56 * t47 + t53 * t48;
t32 = pkin(7) * t49 + t38;
t22 = t27 * t55 - t52 * t32;
t62 = qJD(3) - t22;
t12 = pkin(3) * t36 - t64 * t40 + t62;
t35 = t42 * t52 - t55 * t49;
t37 = -t53 * t47 + t48 * t56;
t31 = -pkin(1) * t49 - t37;
t59 = -qJ(3) * t36 + t31;
t16 = t64 * t35 + t59;
t50 = sin(pkin(9));
t63 = cos(pkin(9));
t9 = t50 * t12 + t63 * t16;
t23 = t52 * t27 + t55 * t32;
t6 = t36 * qJ(5) + t9;
t20 = -t40 * qJ(3) - t23;
t8 = t63 * t12 - t50 * t16;
t61 = qJD(5) - t8;
t60 = pkin(3) * t35 - qJD(4) + t20;
t25 = t50 * t35 + t63 * t40;
t58 = qJ(5) * t25 + t60;
t57 = V_base(3) ^ 2;
t54 = cos(qJ(6));
t51 = sin(qJ(6));
t34 = qJD(6) - t36;
t24 = -t63 * t35 + t40 * t50;
t21 = pkin(2) * t35 + t59;
t19 = -pkin(2) * t40 + t62;
t18 = t24 * t51 + t25 * t54;
t17 = t24 * t54 - t25 * t51;
t10 = pkin(4) * t24 - t58;
t7 = t65 * t24 + t58;
t5 = -t36 * pkin(4) + t61;
t4 = pkin(8) * t24 + t6;
t3 = -t25 * pkin(8) + t65 * t36 + t61;
t2 = t3 * t51 + t4 * t54;
t1 = t3 * t54 - t4 * t51;
t11 = m(2) * (t37 ^ 2 + t38 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t20 ^ 2 + t21 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t37 * mrSges(2,1) - t38 * mrSges(2,2) + Ifges(2,3) * t49 / 0.2e1) * t49 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t37 * mrSges(2,3) + Ifges(2,5) * t49 + Ifges(2,1) * t42 / 0.2e1) * t42 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t18 / 0.2e1) * t18 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t38 * mrSges(2,3) + Ifges(2,4) * t42 + Ifges(2,6) * t49 + Ifges(2,2) * t41 / 0.2e1) * t41 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t18 + Ifges(7,6) * t34 + Ifges(7,2) * t17 / 0.2e1) * t17 + (t22 * mrSges(3,1) - t23 * mrSges(3,2) + t19 * mrSges(4,2) - t20 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t40) * t40 + (-t60 * mrSges(5,2) + t5 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t25) * t25 + (t31 * mrSges(3,1) + t20 * mrSges(4,1) - t21 * mrSges(4,2) - t23 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t35 + (Ifges(4,5) - Ifges(3,6)) * t40) * t35 + (-t60 * mrSges(5,1) + t10 * mrSges(6,1) - t6 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t24 + (-Ifges(5,4) + Ifges(6,5)) * t25) * t24 + (t19 * mrSges(4,1) + t8 * mrSges(5,1) - t5 * mrSges(6,1) + t31 * mrSges(3,2) - t9 * mrSges(5,2) - t22 * mrSges(3,3) - t21 * mrSges(4,3) + t6 * mrSges(6,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t36 + (-Ifges(4,4) + Ifges(3,5)) * t40 + (-Ifges(3,4) - Ifges(4,6)) * t35 + (Ifges(6,4) + Ifges(5,5)) * t25 + (-Ifges(5,6) + Ifges(6,6)) * t24) * t36;
T  = t11;
