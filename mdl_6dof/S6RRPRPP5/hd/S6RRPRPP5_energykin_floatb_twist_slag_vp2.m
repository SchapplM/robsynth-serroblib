% Calculate kinetic energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:24
% EndTime: 2019-03-09 10:03:25
% DurationCPUTime: 0.87s
% Computational Cost: add. (1393->146), mult. (1733->179), div. (0->0), fcn. (1228->6), ass. (0->46)
t59 = pkin(2) + pkin(8);
t58 = -pkin(4) - pkin(5);
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t38 = t48 * V_base(5) + t50 * V_base(4);
t45 = V_base(6) + qJD(1);
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t32 = t38 * t49 + t45 * t47;
t37 = -t48 * V_base(4) + t50 * V_base(5);
t36 = qJD(2) - t37;
t23 = -pkin(1) * t37 - pkin(7) * t38 + V_base(3);
t43 = V_base(5) * pkin(6) + V_base(1);
t44 = -V_base(4) * pkin(6) + V_base(2);
t34 = t50 * t43 + t48 * t44;
t29 = pkin(7) * t45 + t34;
t18 = t23 * t49 - t47 * t29;
t56 = qJD(3) - t18;
t10 = pkin(3) * t32 - t36 * t59 + t56;
t31 = t38 * t47 - t49 * t45;
t33 = -t48 * t43 + t44 * t50;
t28 = -pkin(1) * t45 - t33;
t52 = -qJ(3) * t32 + t28;
t14 = t31 * t59 + t52;
t46 = sin(qJ(4));
t57 = cos(qJ(4));
t7 = t46 * t10 + t57 * t14;
t19 = t47 * t23 + t49 * t29;
t30 = qJD(4) + t32;
t4 = t30 * qJ(5) + t7;
t16 = -t36 * qJ(3) - t19;
t6 = t10 * t57 - t46 * t14;
t55 = pkin(3) * t31 + t16;
t54 = qJD(5) - t6;
t21 = t46 * t31 + t36 * t57;
t53 = qJ(5) * t21 + t55;
t51 = V_base(3) ^ 2;
t20 = -t31 * t57 + t36 * t46;
t17 = pkin(2) * t31 + t52;
t15 = -pkin(2) * t36 + t56;
t8 = pkin(4) * t20 - t53;
t5 = t20 * t58 + qJD(6) + t53;
t3 = -t30 * pkin(4) + t54;
t2 = qJ(6) * t20 + t4;
t1 = -t21 * qJ(6) + t30 * t58 + t54;
t9 = m(2) * (t33 ^ 2 + t34 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(5) * (t55 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t45 / 0.2e1) * t45 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,5) * t45 + Ifges(2,1) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,6) * t45 + Ifges(2,2) * t37 / 0.2e1) * t37 + (t18 * mrSges(3,1) - t19 * mrSges(3,2) + t15 * mrSges(4,2) - t16 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t36) * t36 + (t15 * mrSges(4,1) + t28 * mrSges(3,2) - t18 * mrSges(3,3) - t17 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t32 + (-Ifges(4,4) + Ifges(3,5)) * t36) * t32 + (t28 * mrSges(3,1) + t16 * mrSges(4,1) - t17 * mrSges(4,2) - t19 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t31 + (Ifges(4,5) - Ifges(3,6)) * t36 + (-Ifges(3,4) - Ifges(4,6)) * t32) * t31 + (t6 * mrSges(5,1) - t3 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t30) * t30 + (-t55 * mrSges(5,2) + t3 * mrSges(6,2) + t5 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t21 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t30) * t21 + (-t55 * mrSges(5,1) + t8 * mrSges(6,1) - t5 * mrSges(7,1) - t4 * mrSges(6,2) - t7 * mrSges(5,3) + t2 * mrSges(7,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t20 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t30 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t21) * t20;
T  = t9;
