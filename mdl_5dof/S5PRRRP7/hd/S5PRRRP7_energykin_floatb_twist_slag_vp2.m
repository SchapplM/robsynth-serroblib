% Calculate kinetic energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:01
% EndTime: 2019-12-05 16:54:02
% DurationCPUTime: 0.89s
% Computational Cost: add. (2171->133), mult. (3706->187), div. (0->0), fcn. (3022->10), ass. (0->49)
t47 = V_base(5) * qJ(1) + V_base(1);
t48 = -V_base(4) * qJ(1) + V_base(2);
t50 = sin(pkin(9));
t52 = cos(pkin(9));
t40 = -t47 * t50 + t52 * t48;
t53 = cos(pkin(5));
t43 = t50 * V_base(5) + t52 * V_base(4);
t65 = pkin(6) * t43;
t35 = V_base(6) * pkin(1) - t53 * t65 + t40;
t42 = -t50 * V_base(4) + t52 * V_base(5);
t49 = V_base(3) + qJD(1);
t51 = sin(pkin(5));
t38 = -pkin(1) * t42 - t51 * t65 + t49;
t66 = t35 * t53 + t38 * t51;
t41 = t52 * t47 + t50 * t48;
t60 = t42 * t53 + t51 * V_base(6);
t31 = t60 * pkin(6) + t41;
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t20 = -t56 * t31 + t66 * t59;
t33 = -t56 * t43 + t60 * t59;
t39 = -t42 * t51 + t53 * V_base(6) + qJD(2);
t18 = -t39 * pkin(2) - t20;
t34 = t43 * t59 + t60 * t56;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t26 = -t34 * t55 + t39 * t58;
t27 = t34 * t58 + t39 * t55;
t13 = -t26 * pkin(3) - t27 * pkin(8) + t18;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t24 = -t35 * t51 + t53 * t38;
t15 = -pkin(2) * t33 - pkin(7) * t34 + t24;
t21 = t59 * t31 + t66 * t56;
t19 = pkin(7) * t39 + t21;
t10 = t55 * t15 + t58 * t19;
t32 = qJD(3) - t33;
t8 = pkin(8) * t32 + t10;
t4 = t54 * t13 + t57 * t8;
t3 = t57 * t13 - t54 * t8;
t9 = t15 * t58 - t55 * t19;
t7 = -pkin(3) * t32 - t9;
t25 = qJD(4) - t26;
t23 = t27 * t57 + t32 * t54;
t22 = -t27 * t54 + t32 * t57;
t5 = -pkin(4) * t22 + qJD(5) + t7;
t2 = qJ(5) * t22 + t4;
t1 = pkin(4) * t25 - qJ(5) * t23 + t3;
t6 = m(2) * (t40 ^ 2 + t41 ^ 2 + t49 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t24 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t18 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t49 * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,1) * t43 / 0.2e1) * t43 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + Ifges(3,3) * t39 / 0.2e1) * t39 + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + Ifges(4,3) * t32 / 0.2e1) * t32 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t49 * mrSges(2,1) + t41 * mrSges(2,3) + Ifges(2,4) * t43 + Ifges(2,2) * t42 / 0.2e1) * t42 + (t24 * mrSges(3,2) - t20 * mrSges(3,3) + Ifges(3,5) * t39 + Ifges(3,1) * t34 / 0.2e1) * t34 + (t18 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,5) * t32 + Ifges(4,1) * t27 / 0.2e1) * t27 + (-t24 * mrSges(3,1) + t21 * mrSges(3,3) + Ifges(3,4) * t34 + Ifges(3,6) * t39 + Ifges(3,2) * t33 / 0.2e1) * t33 + (-t18 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t27 + Ifges(4,6) * t32 + Ifges(4,2) * t26 / 0.2e1) * t26 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t25) * t25 + (t7 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t23 + (Ifges(5,5) + Ifges(6,5)) * t25) * t23 + (V_base(2) * mrSges(1,1) + t40 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t41 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t43 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t42 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t7 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t22 + (Ifges(5,6) + Ifges(6,6)) * t25 + (Ifges(5,4) + Ifges(6,4)) * t23) * t22;
T = t6;
