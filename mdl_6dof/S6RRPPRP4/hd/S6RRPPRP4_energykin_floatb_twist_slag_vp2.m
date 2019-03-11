% Calculate kinetic energy for
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:36:55
% EndTime: 2019-03-09 08:36:56
% DurationCPUTime: 0.84s
% Computational Cost: add. (1985->149), mult. (2545->194), div. (0->0), fcn. (1928->8), ass. (0->47)
t57 = sin(qJ(5));
t64 = cos(qJ(5));
t59 = sin(qJ(1));
t60 = cos(qJ(1));
t48 = t59 * V_base(5) + t60 * V_base(4);
t55 = V_base(6) + qJD(1);
t58 = sin(qJ(2));
t65 = cos(qJ(2));
t41 = t48 * t65 + t55 * t58;
t47 = -t59 * V_base(4) + t60 * V_base(5);
t46 = qJD(2) - t47;
t56 = sin(pkin(9));
t63 = cos(pkin(9));
t30 = t41 * t63 + t46 * t56;
t40 = t48 * t58 - t55 * t65;
t33 = -pkin(1) * t47 - pkin(7) * t48 + V_base(3);
t53 = pkin(6) * V_base(5) + V_base(1);
t54 = -pkin(6) * V_base(4) + V_base(2);
t43 = t53 * t60 + t54 * t59;
t37 = pkin(7) * t55 + t43;
t25 = t33 * t58 + t37 * t65;
t22 = qJ(3) * t46 + t25;
t42 = -t53 * t59 + t54 * t60;
t36 = -pkin(1) * t55 - t42;
t23 = pkin(2) * t40 - qJ(3) * t41 + t36;
t14 = -t22 * t56 + t23 * t63;
t62 = qJD(4) - t14;
t7 = -t30 * pkin(8) + (-pkin(3) - pkin(4)) * t40 + t62;
t15 = t22 * t63 + t23 * t56;
t12 = qJ(4) * t40 + t15;
t29 = t41 * t56 - t46 * t63;
t9 = pkin(8) * t29 + t12;
t4 = t57 * t7 + t64 * t9;
t24 = t33 * t65 - t37 * t58;
t21 = -pkin(2) * t46 + qJD(3) - t24;
t13 = pkin(3) * t29 - qJ(4) * t30 + t21;
t3 = -t57 * t9 + t64 * t7;
t10 = -pkin(4) * t29 - t13;
t61 = V_base(3) ^ 2;
t39 = qJD(5) - t40;
t17 = t29 * t57 + t30 * t64;
t16 = -t29 * t64 + t30 * t57;
t11 = -pkin(3) * t40 + t62;
t5 = pkin(5) * t16 - qJ(6) * t17 + t10;
t2 = qJ(6) * t39 + t4;
t1 = -pkin(5) * t39 + qJD(6) - t3;
t6 = m(2) * (t42 ^ 2 + t43 ^ 2 + t61) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t61) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t21 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t13 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t42 * mrSges(2,1) - t43 * mrSges(2,2) + Ifges(2,3) * t55 / 0.2e1) * t55 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t42 * mrSges(2,3) + Ifges(2,5) * t55 + Ifges(2,1) * t48 / 0.2e1) * t48 + (t36 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t41 / 0.2e1) * t41 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t43 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t55 + Ifges(2,2) * t47 / 0.2e1) * t47 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t39) * t39 + (t21 * mrSges(4,2) + t11 * mrSges(5,2) - t14 * mrSges(4,3) - t13 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t30) * t30 + (t21 * mrSges(4,1) + t13 * mrSges(5,1) - t12 * mrSges(5,2) - t15 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29 + (-Ifges(4,4) + Ifges(5,5)) * t30) * t29 + (t10 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t17 + (Ifges(7,4) + Ifges(6,5)) * t39) * t17 + (t10 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t16 + (-Ifges(6,6) + Ifges(7,6)) * t39 + (-Ifges(6,4) + Ifges(7,5)) * t17) * t16 + (t36 * mrSges(3,1) + t14 * mrSges(4,1) - t11 * mrSges(5,1) - t15 * mrSges(4,2) - t25 * mrSges(3,3) + t12 * mrSges(5,3) - Ifges(3,4) * t41 - Ifges(3,6) * t46 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t40 + (Ifges(5,4) + Ifges(4,5)) * t30 + (-Ifges(4,6) + Ifges(5,6)) * t29) * t40;
T  = t6;
