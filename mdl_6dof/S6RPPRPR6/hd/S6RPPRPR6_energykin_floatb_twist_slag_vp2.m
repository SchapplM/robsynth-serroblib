% Calculate kinetic energy for
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:32
% EndTime: 2019-03-09 01:50:33
% DurationCPUTime: 0.69s
% Computational Cost: add. (991->148), mult. (1235->179), div. (0->0), fcn. (760->6), ass. (0->46)
t58 = pkin(4) + pkin(8);
t57 = cos(qJ(1));
t56 = cos(qJ(4));
t49 = sin(qJ(1));
t34 = t49 * V_base(4) - t57 * V_base(5);
t46 = V_base(6) + qJD(1);
t40 = V_base(5) * pkin(6) + V_base(1);
t41 = -V_base(4) * pkin(6) + V_base(2);
t29 = t57 * t40 + t49 * t41;
t24 = -t46 * qJ(2) - t29;
t55 = qJD(3) - t24;
t13 = -pkin(7) * t46 + (-pkin(2) - pkin(3)) * t34 + t55;
t31 = t34 * qJ(3);
t35 = t49 * V_base(5) + t57 * V_base(4);
t54 = pkin(1) * t34 + V_base(3);
t16 = t31 + (pkin(7) - qJ(2)) * t35 + t54;
t48 = sin(qJ(4));
t9 = t48 * t13 + t56 * t16;
t28 = -t49 * t40 + t57 * t41;
t23 = -t46 * pkin(1) + qJD(2) - t28;
t33 = qJD(4) - t34;
t6 = -qJ(5) * t33 - t9;
t8 = t13 * t56 - t48 * t16;
t17 = t35 * pkin(2) - t46 * qJ(3) + t23;
t53 = qJD(5) - t8;
t27 = t48 * t35 + t46 * t56;
t22 = -qJ(2) * t35 + t54;
t12 = -pkin(3) * t35 - t17;
t52 = -qJ(5) * t27 + t12;
t51 = V_base(3) ^ 2;
t50 = cos(qJ(6));
t47 = sin(qJ(6));
t26 = -t35 * t56 + t46 * t48;
t25 = qJD(6) + t27;
t21 = -pkin(2) * t34 + t55;
t20 = t22 + t31;
t19 = t26 * t47 + t33 * t50;
t18 = t26 * t50 - t33 * t47;
t10 = pkin(4) * t26 + t52;
t7 = t26 * t58 + t52;
t5 = -t33 * pkin(4) + t53;
t4 = -pkin(5) * t26 - t6;
t3 = t27 * pkin(5) - t33 * t58 + t53;
t2 = t3 * t47 + t50 * t7;
t1 = t3 * t50 - t47 * t7;
t11 = m(2) * (t28 ^ 2 + t29 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(4) * (t17 ^ 2 + t20 ^ 2 + t21 ^ 2) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t25 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t8 * mrSges(5,1) - t9 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t33) * t33 + (t5 * mrSges(6,1) + t12 * mrSges(5,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t27 + (-Ifges(6,4) + Ifges(5,5)) * t33) * t27 + (t12 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t26 + (Ifges(6,5) - Ifges(5,6)) * t33 + (-Ifges(5,4) - Ifges(6,6)) * t27) * t26 + (t28 * mrSges(2,1) - t29 * mrSges(2,2) + t23 * mrSges(3,2) + t21 * mrSges(4,2) - t24 * mrSges(3,3) - t17 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t46) * t46 + (t23 * mrSges(3,1) + t17 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t20 * mrSges(4,2) - t28 * mrSges(2,3) - t22 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t35 + (-Ifges(3,4) + Ifges(2,5) + Ifges(4,5)) * t46) * t35 + (V_base(3) * mrSges(2,1) + t24 * mrSges(3,1) - t21 * mrSges(4,1) - t22 * mrSges(3,2) - t29 * mrSges(2,3) + t20 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t34 + (Ifges(4,4) + Ifges(3,5) - Ifges(2,6)) * t46 + (-Ifges(2,4) - Ifges(3,6) + Ifges(4,6)) * t35) * t34;
T  = t11;
