% Calculate kinetic energy for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:17
% DurationCPUTime: 0.60s
% Computational Cost: add. (949->125), mult. (1266->159), div. (0->0), fcn. (872->6), ass. (0->37)
t47 = sin(qJ(1));
t52 = cos(qJ(1));
t33 = t47 * V_base(4) - t52 * V_base(5);
t34 = t47 * V_base(5) + t52 * V_base(4);
t19 = pkin(1) * t33 - qJ(2) * t34 + V_base(3);
t40 = V_base(5) * pkin(5) + V_base(1);
t41 = -V_base(4) * pkin(5) + V_base(2);
t30 = t52 * t40 + t47 * t41;
t44 = V_base(6) + qJD(1);
t24 = qJ(2) * t44 + t30;
t45 = sin(pkin(7));
t50 = cos(pkin(7));
t15 = t45 * t19 + t50 * t24;
t13 = t33 * qJ(3) + t15;
t27 = t34 * t45 - t50 * t44;
t10 = pkin(6) * t27 + t13;
t46 = sin(qJ(4));
t51 = cos(qJ(4));
t28 = t50 * t34 + t45 * t44;
t14 = t50 * t19 - t45 * t24;
t49 = qJD(3) - t14;
t7 = -t28 * pkin(6) + (-pkin(2) - pkin(3)) * t33 + t49;
t4 = t51 * t10 + t46 * t7;
t29 = -t47 * t40 + t52 * t41;
t22 = -t44 * pkin(1) + qJD(2) - t29;
t11 = t27 * pkin(2) - t28 * qJ(3) + t22;
t3 = -t46 * t10 + t51 * t7;
t9 = -pkin(3) * t27 - t11;
t48 = V_base(3) ^ 2;
t32 = qJD(4) - t33;
t17 = t46 * t27 + t51 * t28;
t16 = -t51 * t27 + t28 * t46;
t12 = -t33 * pkin(2) + t49;
t5 = pkin(4) * t16 - qJ(5) * t17 + t9;
t2 = qJ(5) * t32 + t4;
t1 = -t32 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t29 ^ 2 + t30 ^ 2 + t48) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t48) / 0.2e1 + m(3) * (t14 ^ 2 + t15 ^ 2 + t22 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t13 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t29 * mrSges(2,1) - t30 * mrSges(2,2) + Ifges(2,3) * t44 / 0.2e1) * t44 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t29 * mrSges(2,3) + Ifges(2,5) * t44 + Ifges(2,1) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t3 * mrSges(5,1) - t1 * mrSges(6,1) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t32) * t32 + (t22 * mrSges(3,2) + t12 * mrSges(4,2) - t14 * mrSges(3,3) - t11 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t28) * t28 + (t22 * mrSges(3,1) + t11 * mrSges(4,1) - t13 * mrSges(4,2) - t15 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t27 + (-Ifges(3,4) + Ifges(4,5)) * t28) * t27 + (t9 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t17 + (Ifges(6,4) + Ifges(5,5)) * t32) * t17 + (t9 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t16 + (-Ifges(5,6) + Ifges(6,6)) * t32 + (-Ifges(5,4) + Ifges(6,5)) * t17) * t16 + (V_base(3) * mrSges(2,1) + t14 * mrSges(3,1) - t12 * mrSges(4,1) - t15 * mrSges(3,2) - t30 * mrSges(2,3) + t13 * mrSges(4,3) - Ifges(2,4) * t34 - Ifges(2,6) * t44 + (Ifges(2,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t33 + (Ifges(4,4) + Ifges(3,5)) * t28 + (-Ifges(3,6) + Ifges(4,6)) * t27) * t33;
T = t6;
