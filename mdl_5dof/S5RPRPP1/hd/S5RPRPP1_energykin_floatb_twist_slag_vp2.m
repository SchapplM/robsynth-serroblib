% Calculate kinetic energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:35
% EndTime: 2019-12-31 18:08:35
% DurationCPUTime: 0.68s
% Computational Cost: add. (1505->128), mult. (2200->174), div. (0->0), fcn. (1684->8), ass. (0->42)
t45 = sin(pkin(8));
t53 = cos(pkin(8));
t42 = V_base(5) * pkin(5) + V_base(1);
t43 = -V_base(4) * pkin(5) + V_base(2);
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t34 = -t42 * t49 + t51 * t43;
t38 = t49 * V_base(5) + t51 * V_base(4);
t44 = V_base(6) + qJD(1);
t27 = pkin(1) * t44 - qJ(2) * t38 + t34;
t35 = t51 * t42 + t49 * t43;
t37 = -t49 * V_base(4) + t51 * V_base(5);
t30 = qJ(2) * t37 + t35;
t46 = sin(pkin(7));
t47 = cos(pkin(7));
t22 = t46 * t27 + t47 * t30;
t17 = pkin(6) * t44 + t22;
t32 = t37 * t47 - t38 * t46;
t33 = t37 * t46 + t38 * t47;
t36 = -pkin(1) * t37 + qJD(2) + V_base(3);
t20 = -pkin(2) * t32 - pkin(6) * t33 + t36;
t48 = sin(qJ(3));
t50 = cos(qJ(3));
t10 = -t17 * t48 + t50 * t20;
t26 = t33 * t50 + t44 * t48;
t31 = qJD(3) - t32;
t7 = pkin(3) * t31 - qJ(4) * t26 + t10;
t11 = t50 * t17 + t48 * t20;
t25 = -t33 * t48 + t44 * t50;
t9 = qJ(4) * t25 + t11;
t4 = t45 * t7 + t53 * t9;
t21 = t27 * t47 - t46 * t30;
t16 = -pkin(2) * t44 - t21;
t3 = -t45 * t9 + t53 * t7;
t12 = -pkin(3) * t25 + qJD(4) + t16;
t52 = V_base(3) ^ 2;
t14 = t45 * t25 + t53 * t26;
t13 = -t53 * t25 + t26 * t45;
t5 = pkin(4) * t13 - qJ(5) * t14 + t12;
t2 = qJ(5) * t31 + t4;
t1 = -t31 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t34 ^ 2 + t35 ^ 2 + t52) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t52) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t11 ^ 2 + t16 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t34 * mrSges(2,3) + Ifges(2,1) * t38 / 0.2e1) * t38 + (t36 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,1) * t33 / 0.2e1) * t33 + (t16 * mrSges(4,2) - t10 * mrSges(4,3) + Ifges(4,1) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,2) * t37 / 0.2e1) * t37 + (-t36 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t33 + Ifges(3,2) * t32 / 0.2e1) * t32 + (-t16 * mrSges(4,1) + t11 * mrSges(4,3) + Ifges(4,4) * t26 + Ifges(4,2) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t12 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t14) * t14 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t13 + (-Ifges(5,4) + Ifges(6,5)) * t14) * t13 + (t34 * mrSges(2,1) + t21 * mrSges(3,1) - t35 * mrSges(2,2) - t22 * mrSges(3,2) + Ifges(2,5) * t38 + Ifges(3,5) * t33 + Ifges(2,6) * t37 + Ifges(3,6) * t32 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * t44) * t44 + (t10 * mrSges(4,1) + t3 * mrSges(5,1) - t1 * mrSges(6,1) - t11 * mrSges(4,2) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + Ifges(4,5) * t26 + Ifges(4,6) * t25 + (Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t31 + (Ifges(6,4) + Ifges(5,5)) * t14 + (-Ifges(5,6) + Ifges(6,6)) * t13) * t31;
T = t6;
