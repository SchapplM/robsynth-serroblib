% Calculate kinetic energy for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:20
% EndTime: 2019-12-31 20:49:20
% DurationCPUTime: 0.71s
% Computational Cost: add. (1603->128), mult. (2200->176), div. (0->0), fcn. (1684->8), ass. (0->43)
t46 = sin(pkin(8));
t54 = cos(pkin(8));
t42 = V_base(5) * pkin(5) + V_base(1);
t43 = -V_base(4) * pkin(5) + V_base(2);
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t34 = -t42 * t49 + t52 * t43;
t38 = t49 * V_base(5) + t52 * V_base(4);
t45 = V_base(6) + qJD(1);
t27 = pkin(1) * t45 - pkin(6) * t38 + t34;
t35 = t52 * t42 + t49 * t43;
t37 = -t49 * V_base(4) + t52 * V_base(5);
t30 = pkin(6) * t37 + t35;
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t22 = t48 * t27 + t51 * t30;
t44 = qJD(2) + t45;
t17 = pkin(7) * t44 + t22;
t32 = t37 * t51 - t38 * t48;
t33 = t37 * t48 + t38 * t51;
t36 = -pkin(1) * t37 + V_base(3);
t20 = -pkin(2) * t32 - pkin(7) * t33 + t36;
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t10 = -t17 * t47 + t50 * t20;
t25 = t33 * t50 + t44 * t47;
t31 = qJD(3) - t32;
t7 = pkin(3) * t31 - qJ(4) * t25 + t10;
t11 = t50 * t17 + t47 * t20;
t24 = -t33 * t47 + t44 * t50;
t9 = qJ(4) * t24 + t11;
t4 = t46 * t7 + t54 * t9;
t21 = t27 * t51 - t48 * t30;
t16 = -pkin(2) * t44 - t21;
t3 = -t46 * t9 + t54 * t7;
t12 = -pkin(3) * t24 + qJD(4) + t16;
t53 = V_base(3) ^ 2;
t14 = t46 * t24 + t54 * t25;
t13 = -t54 * t24 + t25 * t46;
t5 = pkin(4) * t13 - qJ(5) * t14 + t12;
t2 = qJ(5) * t31 + t4;
t1 = -t31 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t34 ^ 2 + t35 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t11 ^ 2 + t16 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t34 * mrSges(2,1) - t35 * mrSges(2,2) + Ifges(2,3) * t45 / 0.2e1) * t45 + (t21 * mrSges(3,1) - t22 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (t16 * mrSges(4,2) - t10 * mrSges(4,3) + Ifges(4,1) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t34 * mrSges(2,3) + Ifges(2,5) * t45 + Ifges(2,1) * t38 / 0.2e1) * t38 + (t36 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,5) * t44 + Ifges(3,1) * t33 / 0.2e1) * t33 + (-t16 * mrSges(4,1) + t11 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,2) * t24 / 0.2e1) * t24 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,6) * t45 + Ifges(2,2) * t37 / 0.2e1) * t37 + (-t36 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t33 + Ifges(3,6) * t44 + Ifges(3,2) * t32 / 0.2e1) * t32 + (t12 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t14) * t14 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t13 + (-Ifges(5,4) + Ifges(6,5)) * t14) * t13 + (t10 * mrSges(4,1) + t3 * mrSges(5,1) - t1 * mrSges(6,1) - t11 * mrSges(4,2) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + Ifges(4,5) * t25 + Ifges(4,6) * t24 + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t31 + (Ifges(6,4) + Ifges(5,5)) * t14 + (-Ifges(5,6) + Ifges(6,6)) * t13) * t31;
T = t6;
