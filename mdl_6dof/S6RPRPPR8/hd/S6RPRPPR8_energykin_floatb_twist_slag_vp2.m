% Calculate kinetic energy for
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:41
% EndTime: 2019-03-09 02:58:41
% DurationCPUTime: 0.71s
% Computational Cost: add. (1075->148), mult. (1325->179), div. (0->0), fcn. (848->6), ass. (0->47)
t59 = pkin(1) + pkin(7);
t58 = -pkin(4) - pkin(8);
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t37 = t48 * V_base(5) + t51 * V_base(4);
t45 = V_base(6) + qJD(1);
t41 = V_base(5) * pkin(6) + V_base(1);
t42 = -V_base(4) * pkin(6) + V_base(2);
t30 = -t48 * t41 + t42 * t51;
t56 = qJD(2) - t30;
t15 = pkin(2) * t37 - t45 * t59 + t56;
t36 = t48 * V_base(4) - t51 * V_base(5);
t57 = -qJ(2) * t37 + V_base(3);
t20 = t36 * t59 + t57;
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t12 = t47 * t15 + t50 * t20;
t31 = t51 * t41 + t48 * t42;
t35 = qJD(3) + t37;
t10 = t35 * qJ(4) + t12;
t24 = -t45 * qJ(2) - t31;
t11 = t15 * t50 - t47 * t20;
t21 = -t36 * pkin(2) - t24;
t55 = qJD(4) - t11;
t29 = t36 * t47 + t45 * t50;
t28 = -t50 * t36 + t45 * t47;
t7 = -qJ(5) * t28 - t10;
t13 = t28 * pkin(3) - t29 * qJ(4) + t21;
t54 = qJD(5) - t13;
t53 = -qJ(5) * t29 + t55;
t52 = V_base(3) ^ 2;
t49 = cos(qJ(6));
t46 = sin(qJ(6));
t27 = qJD(6) + t29;
t23 = -pkin(1) * t45 + t56;
t22 = pkin(1) * t36 + t57;
t19 = t28 * t49 - t35 * t46;
t18 = -t28 * t46 - t35 * t49;
t9 = -pkin(3) * t35 + t55;
t8 = -pkin(4) * t28 + t54;
t6 = pkin(5) * t35 - t7;
t5 = (-pkin(3) - pkin(4)) * t35 + t53;
t4 = pkin(5) * t29 + t28 * t58 + t54;
t3 = (-pkin(3) + t58) * t35 + t53;
t2 = t3 * t49 + t4 * t46;
t1 = -t3 * t46 + t4 * t49;
t14 = m(2) * (t30 ^ 2 + t31 ^ 2 + t52) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t52) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t21 ^ 2) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t13 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t6 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t6 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t27 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t30 * mrSges(2,1) - t31 * mrSges(2,2) + t23 * mrSges(3,2) - t24 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t45) * t45 + (t23 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t30 * mrSges(2,3) - t22 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t37 + (-Ifges(3,4) + Ifges(2,5)) * t45) * t37 + (V_base(3) * mrSges(2,1) + t24 * mrSges(3,1) - t22 * mrSges(3,2) - t31 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t36 + (Ifges(3,5) - Ifges(2,6)) * t45 + (-Ifges(2,4) - Ifges(3,6)) * t37) * t36 + (t11 * mrSges(4,1) - t9 * mrSges(5,1) - t7 * mrSges(6,1) - t12 * mrSges(4,2) + t5 * mrSges(6,2) + t10 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t35) * t35 + (t8 * mrSges(6,1) + t21 * mrSges(4,2) + t9 * mrSges(5,2) - t11 * mrSges(4,3) - t13 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t29 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t35) * t29 + (t21 * mrSges(4,1) + t13 * mrSges(5,1) - t10 * mrSges(5,2) + t8 * mrSges(6,2) - t12 * mrSges(4,3) - t7 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t28 + (-Ifges(6,5) - Ifges(4,6) + Ifges(5,6)) * t35 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t29) * t28;
T  = t14;
