% Calculate kinetic energy for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:00
% EndTime: 2019-12-31 17:52:00
% DurationCPUTime: 0.60s
% Computational Cost: add. (897->125), mult. (1208->159), div. (0->0), fcn. (800->6), ass. (0->37)
t45 = sin(qJ(1));
t49 = cos(qJ(1));
t33 = t45 * V_base(5) + t49 * V_base(4);
t41 = V_base(6) + qJD(1);
t37 = V_base(5) * pkin(5) + V_base(1);
t38 = -V_base(4) * pkin(5) + V_base(2);
t28 = -t45 * t37 + t49 * t38;
t48 = qJD(2) - t28;
t15 = -t33 * qJ(3) + (-pkin(1) - pkin(2)) * t41 + t48;
t29 = t49 * t37 + t45 * t38;
t27 = t41 * qJ(2) + t29;
t32 = t45 * V_base(4) - t49 * V_base(5);
t21 = qJ(3) * t32 + t27;
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t13 = t42 * t15 + t43 * t21;
t11 = -pkin(6) * t41 + t13;
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t25 = t32 * pkin(1) - t33 * qJ(2) + V_base(3);
t16 = -pkin(2) * t32 + qJD(3) - t25;
t23 = t32 * t43 - t33 * t42;
t24 = t32 * t42 + t33 * t43;
t8 = -pkin(3) * t23 - pkin(6) * t24 + t16;
t4 = t46 * t11 + t44 * t8;
t3 = -t11 * t44 + t46 * t8;
t12 = t15 * t43 - t42 * t21;
t10 = pkin(3) * t41 - t12;
t47 = V_base(3) ^ 2;
t26 = -t41 * pkin(1) + t48;
t22 = qJD(4) - t23;
t20 = t24 * t46 - t41 * t44;
t19 = -t24 * t44 - t41 * t46;
t5 = -pkin(4) * t19 + qJD(5) + t10;
t2 = qJ(5) * t19 + t4;
t1 = pkin(4) * t22 - qJ(5) * t20 + t3;
t6 = m(2) * (t28 ^ 2 + t29 ^ 2 + t47) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t47) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t16 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,1) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t16 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t24 + Ifges(4,2) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (V_base(3) * mrSges(2,2) + t26 * mrSges(3,2) - t28 * mrSges(2,3) - t25 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t33) * t33 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t22) * t22 + (V_base(3) * mrSges(2,1) + t25 * mrSges(3,1) - t27 * mrSges(3,2) - t29 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t32 + (-Ifges(2,4) + Ifges(3,5)) * t33) * t32 + (t10 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t20 + (Ifges(5,5) + Ifges(6,5)) * t22) * t20 + (-t10 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t19 + (Ifges(5,6) + Ifges(6,6)) * t22 + (Ifges(5,4) + Ifges(6,4)) * t20) * t19 + (t28 * mrSges(2,1) - t26 * mrSges(3,1) - t12 * mrSges(4,1) - t29 * mrSges(2,2) + t13 * mrSges(4,2) + t27 * mrSges(3,3) - Ifges(4,5) * t24 - Ifges(4,6) * t23 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * t41 + (Ifges(3,4) + Ifges(2,5)) * t33 + (-Ifges(2,6) + Ifges(3,6)) * t32) * t41;
T = t6;
