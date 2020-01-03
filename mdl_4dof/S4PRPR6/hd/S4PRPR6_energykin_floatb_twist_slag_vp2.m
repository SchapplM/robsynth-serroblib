% Calculate kinetic energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:16
% EndTime: 2019-12-31 16:24:17
% DurationCPUTime: 0.52s
% Computational Cost: add. (841->108), mult. (1245->157), div. (0->0), fcn. (896->8), ass. (0->39)
t48 = cos(qJ(2));
t42 = sin(pkin(6));
t44 = cos(pkin(6));
t32 = -t42 * V_base(4) + t44 * V_base(5);
t33 = t42 * V_base(5) + t44 * V_base(4);
t40 = V_base(3) + qJD(1);
t21 = -pkin(1) * t32 - pkin(4) * t33 + t40;
t37 = V_base(5) * qJ(1) + V_base(1);
t38 = -V_base(4) * qJ(1) + V_base(2);
t30 = t44 * t37 + t42 * t38;
t25 = V_base(6) * pkin(4) + t30;
t46 = sin(qJ(2));
t17 = t46 * t21 + t48 * t25;
t31 = qJD(2) - t32;
t12 = qJ(3) * t31 + t17;
t29 = -t42 * t37 + t38 * t44;
t24 = -V_base(6) * pkin(1) - t29;
t27 = t33 * t46 - t48 * V_base(6);
t28 = t48 * t33 + t46 * V_base(6);
t15 = pkin(2) * t27 - qJ(3) * t28 + t24;
t41 = sin(pkin(7));
t43 = cos(pkin(7));
t6 = t43 * t12 + t41 * t15;
t5 = -t12 * t41 + t43 * t15;
t16 = t48 * t21 - t46 * t25;
t11 = -t31 * pkin(2) + qJD(3) - t16;
t47 = cos(qJ(4));
t45 = sin(qJ(4));
t26 = qJD(4) + t27;
t19 = t28 * t43 + t31 * t41;
t18 = -t28 * t41 + t31 * t43;
t9 = t45 * t18 + t19 * t47;
t8 = t18 * t47 - t45 * t19;
t7 = -t18 * pkin(3) + t11;
t4 = pkin(5) * t18 + t6;
t3 = pkin(3) * t27 - pkin(5) * t19 + t5;
t2 = t45 * t3 + t4 * t47;
t1 = t3 * t47 - t45 * t4;
t10 = m(2) * (t29 ^ 2 + t30 ^ 2 + t40 ^ 2) / 0.2e1 + m(3) * (t16 ^ 2 + t17 ^ 2 + t24 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t7 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t9 / 0.2e1) * t9 + (t40 * mrSges(2,2) - t29 * mrSges(2,3) + Ifges(2,1) * t33 / 0.2e1) * t33 + (t16 * mrSges(3,1) - t17 * mrSges(3,2) + Ifges(3,3) * t31 / 0.2e1) * t31 + (t11 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,2) * t8 / 0.2e1) * t8 + (-t40 * mrSges(2,1) + t30 * mrSges(2,3) + Ifges(2,4) * t33 + Ifges(2,2) * t32 / 0.2e1) * t32 + (t24 * mrSges(3,2) - t16 * mrSges(3,3) + Ifges(3,5) * t31 + Ifges(3,1) * t28 / 0.2e1) * t28 + (-t11 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t19 + Ifges(4,2) * t18 / 0.2e1) * t18 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t8 + Ifges(5,3) * t26 / 0.2e1) * t26 + (V_base(2) * mrSges(1,1) + t29 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t30 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t33 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t32 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t24 * mrSges(3,1) + t5 * mrSges(4,1) - t6 * mrSges(4,2) - t17 * mrSges(3,3) - Ifges(3,4) * t28 + Ifges(4,5) * t19 - Ifges(3,6) * t31 + Ifges(4,6) * t18 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t27) * t27;
T = t10;
