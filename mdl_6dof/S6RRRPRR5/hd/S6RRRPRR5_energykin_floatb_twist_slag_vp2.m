% Calculate kinetic energy for
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:19
% EndTime: 2019-03-09 18:20:21
% DurationCPUTime: 1.08s
% Computational Cost: add. (2641->153), mult. (3359->213), div. (0->0), fcn. (2644->10), ass. (0->57)
t71 = pkin(3) + pkin(9);
t62 = sin(qJ(1));
t67 = cos(qJ(1));
t50 = t62 * V_base(5) + t67 * V_base(4);
t57 = V_base(6) + qJD(1);
t61 = sin(qJ(2));
t66 = cos(qJ(2));
t42 = -t50 * t61 + t57 * t66;
t43 = t50 * t66 + t57 * t61;
t60 = sin(qJ(3));
t65 = cos(qJ(3));
t33 = t42 * t60 + t43 * t65;
t49 = -t62 * V_base(4) + t67 * V_base(5);
t48 = qJD(2) - t49;
t47 = qJD(3) + t48;
t37 = -pkin(1) * t49 - pkin(7) * t50 + V_base(3);
t54 = V_base(5) * pkin(6) + V_base(1);
t55 = -V_base(4) * pkin(6) + V_base(2);
t45 = t67 * t54 + t62 * t55;
t41 = pkin(7) * t57 + t45;
t28 = t66 * t37 - t41 * t61;
t22 = pkin(2) * t48 - pkin(8) * t43 + t28;
t29 = t61 * t37 + t66 * t41;
t25 = pkin(8) * t42 + t29;
t16 = t22 * t65 - t60 * t25;
t70 = qJD(4) - t16;
t10 = pkin(4) * t33 - t47 * t71 + t70;
t32 = -t65 * t42 + t43 * t60;
t44 = -t62 * t54 + t55 * t67;
t40 = -pkin(1) * t57 - t44;
t34 = -pkin(2) * t42 + t40;
t69 = -qJ(4) * t33 + t34;
t13 = t32 * t71 + t69;
t59 = sin(qJ(5));
t64 = cos(qJ(5));
t6 = t59 * t10 + t64 * t13;
t17 = t60 * t22 + t65 * t25;
t15 = -t47 * qJ(4) - t17;
t5 = t64 * t10 - t13 * t59;
t11 = -pkin(4) * t32 - t15;
t31 = qJD(5) + t33;
t68 = V_base(3) ^ 2;
t63 = cos(qJ(6));
t58 = sin(qJ(6));
t30 = qJD(6) + t31;
t27 = t32 * t59 + t47 * t64;
t26 = t32 * t64 - t47 * t59;
t20 = t26 * t58 + t27 * t63;
t19 = t26 * t63 - t27 * t58;
t18 = pkin(3) * t32 + t69;
t14 = -pkin(3) * t47 + t70;
t7 = -pkin(5) * t26 + t11;
t4 = pkin(10) * t26 + t6;
t3 = pkin(5) * t31 - pkin(10) * t27 + t5;
t2 = t3 * t58 + t4 * t63;
t1 = t3 * t63 - t4 * t58;
t8 = (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + t20 * Ifges(7,4) + t30 * Ifges(7,6) + Ifges(7,2) * t19 / 0.2e1) * t19 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + t14 * mrSges(5,2) - t15 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t47) * t47 + (-t40 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t43 + Ifges(3,6) * t48 + Ifges(3,2) * t42 / 0.2e1) * t42 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t30 + Ifges(7,1) * t20 / 0.2e1) * t20 + (t28 * mrSges(3,1) - t29 * mrSges(3,2) + Ifges(3,3) * t48 / 0.2e1) * t48 + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t50 / 0.2e1) * t50 + (t40 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,5) * t48 + Ifges(3,1) * t43 / 0.2e1) * t43 + (t11 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t31 + Ifges(6,1) * t27 / 0.2e1) * t27 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,6) * t57 + Ifges(2,2) * t49 / 0.2e1) * t49 + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(2) * (t44 ^ 2 + t45 ^ 2 + t68) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t68) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t34 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t18 ^ 2) / 0.2e1 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t30 / 0.2e1) * t30 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t11 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t31 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t34 * mrSges(4,1) + t15 * mrSges(5,1) - t18 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t32 + (Ifges(5,5) - Ifges(4,6)) * t47 + (-Ifges(4,4) - Ifges(5,6)) * t33) * t32 + (t14 * mrSges(5,1) + t34 * mrSges(4,2) - t16 * mrSges(4,3) - t18 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t33 + (-Ifges(5,4) + Ifges(4,5)) * t47) * t33;
T  = t8;
