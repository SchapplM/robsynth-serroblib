% Calculate kinetic energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:15
% EndTime: 2019-12-31 20:03:16
% DurationCPUTime: 0.65s
% Computational Cost: add. (1005->125), mult. (1282->161), div. (0->0), fcn. (888->6), ass. (0->38)
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t34 = -t48 * V_base(4) + t50 * V_base(5);
t35 = t48 * V_base(5) + t50 * V_base(4);
t20 = -pkin(1) * t34 - pkin(6) * t35 + V_base(3);
t41 = V_base(5) * pkin(5) + V_base(1);
t42 = -V_base(4) * pkin(5) + V_base(2);
t30 = t50 * t41 + t48 * t42;
t45 = V_base(6) + qJD(1);
t24 = pkin(6) * t45 + t30;
t47 = sin(qJ(2));
t53 = cos(qJ(2));
t16 = t47 * t20 + t53 * t24;
t33 = qJD(2) - t34;
t13 = t33 * qJ(3) + t16;
t27 = t35 * t47 - t53 * t45;
t10 = pkin(7) * t27 + t13;
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t28 = t53 * t35 + t47 * t45;
t15 = t53 * t20 - t47 * t24;
t52 = qJD(3) - t15;
t8 = -t28 * pkin(7) + (-pkin(2) - pkin(3)) * t33 + t52;
t4 = t49 * t10 + t46 * t8;
t29 = -t48 * t41 + t50 * t42;
t23 = -t45 * pkin(1) - t29;
t3 = -t10 * t46 + t49 * t8;
t14 = t27 * pkin(2) - t28 * qJ(3) + t23;
t11 = -pkin(3) * t27 - t14;
t51 = V_base(3) ^ 2;
t32 = qJD(4) - t33;
t18 = t27 * t46 + t28 * t49;
t17 = t27 * t49 - t28 * t46;
t12 = -t33 * pkin(2) + t52;
t5 = -pkin(4) * t17 + qJD(5) + t11;
t2 = qJ(5) * t17 + t4;
t1 = pkin(4) * t32 - qJ(5) * t18 + t3;
t6 = m(2) * (t29 ^ 2 + t30 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(3) * (t15 ^ 2 + t16 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t29 * mrSges(2,1) - t30 * mrSges(2,2) + Ifges(2,3) * t45 / 0.2e1) * t45 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t29 * mrSges(2,3) + Ifges(2,5) * t45 + Ifges(2,1) * t35 / 0.2e1) * t35 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t30 * mrSges(2,3) + Ifges(2,4) * t35 + Ifges(2,6) * t45 + Ifges(2,2) * t34 / 0.2e1) * t34 + (t15 * mrSges(3,1) - t12 * mrSges(4,1) - t16 * mrSges(3,2) + t13 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t33) * t33 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t32) * t32 + (t23 * mrSges(3,2) + t12 * mrSges(4,2) - t15 * mrSges(3,3) - t14 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t28 + (Ifges(4,4) + Ifges(3,5)) * t33) * t28 + (t11 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t18 + (Ifges(5,5) + Ifges(6,5)) * t32) * t18 + (t23 * mrSges(3,1) + t14 * mrSges(4,1) - t13 * mrSges(4,2) - t16 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t27 + (-Ifges(3,6) + Ifges(4,6)) * t33 + (-Ifges(3,4) + Ifges(4,5)) * t28) * t27 + (-t11 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t17 + (Ifges(5,6) + Ifges(6,6)) * t32 + (Ifges(5,4) + Ifges(6,4)) * t18) * t17;
T = t6;
