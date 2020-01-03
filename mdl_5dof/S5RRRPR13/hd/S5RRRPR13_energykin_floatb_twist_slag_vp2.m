% Calculate kinetic energy for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR13_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:12
% EndTime: 2019-12-31 21:43:13
% DurationCPUTime: 0.94s
% Computational Cost: add. (2157->134), mult. (3228->188), div. (0->0), fcn. (2604->10), ass. (0->53)
t45 = pkin(6) * V_base(5) + V_base(1);
t46 = -pkin(6) * V_base(4) + V_base(2);
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t38 = -t45 * t53 + t56 * t46;
t47 = V_base(6) + qJD(1);
t49 = cos(pkin(5));
t41 = t53 * V_base(5) + t56 * V_base(4);
t66 = pkin(7) * t41;
t32 = pkin(1) * t47 - t49 * t66 + t38;
t40 = -t53 * V_base(4) + t56 * V_base(5);
t48 = sin(pkin(5));
t35 = -pkin(1) * t40 - t48 * t66 + V_base(3);
t68 = t32 * t49 + t35 * t48;
t39 = t56 * t45 + t53 * t46;
t60 = t40 * t49 + t47 * t48;
t29 = pkin(7) * t60 + t39;
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t17 = -t52 * t29 + t55 * t68;
t30 = -t41 * t52 + t55 * t60;
t67 = pkin(3) + pkin(9);
t65 = cos(qJ(3));
t21 = -t32 * t48 + t49 * t35;
t31 = t41 * t55 + t52 * t60;
t12 = -pkin(2) * t30 - pkin(8) * t31 + t21;
t18 = t55 * t29 + t52 * t68;
t37 = -t40 * t48 + t47 * t49 + qJD(2);
t16 = pkin(8) * t37 + t18;
t51 = sin(qJ(3));
t9 = t51 * t12 + t65 * t16;
t28 = qJD(3) - t30;
t7 = -qJ(4) * t28 - t9;
t8 = t12 * t65 - t51 * t16;
t59 = qJD(4) - t8;
t24 = t31 * t65 + t51 * t37;
t15 = -pkin(2) * t37 - t17;
t58 = -qJ(4) * t24 + t15;
t57 = V_base(3) ^ 2;
t54 = cos(qJ(5));
t50 = sin(qJ(5));
t23 = t31 * t51 - t37 * t65;
t22 = qJD(5) + t24;
t20 = t23 * t50 + t28 * t54;
t19 = t23 * t54 - t28 * t50;
t10 = pkin(3) * t23 + t58;
t6 = -t28 * pkin(3) + t59;
t5 = t23 * t67 + t58;
t4 = -pkin(4) * t23 - t7;
t3 = t24 * pkin(4) - t28 * t67 + t59;
t2 = t3 * t50 + t5 * t54;
t1 = t3 * t54 - t5 * t50;
t11 = m(2) * (t38 ^ 2 + t39 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(4) * (t15 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(3) * (t17 ^ 2 + t18 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + Ifges(2,3) * t47 / 0.2e1) * t47 + (t17 * mrSges(3,1) - t18 * mrSges(3,2) + Ifges(3,3) * t37 / 0.2e1) * t37 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t22 / 0.2e1) * t22 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,5) * t47 + Ifges(2,1) * t41 / 0.2e1) * t41 + (t21 * mrSges(3,2) - t17 * mrSges(3,3) + Ifges(3,5) * t37 + Ifges(3,1) * t31 / 0.2e1) * t31 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t22 + Ifges(6,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t41 + Ifges(2,6) * t47 + Ifges(2,2) * t40 / 0.2e1) * t40 + (-t21 * mrSges(3,1) + t18 * mrSges(3,3) + Ifges(3,4) * t31 + Ifges(3,6) * t37 + Ifges(3,2) * t30 / 0.2e1) * t30 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t20 + Ifges(6,6) * t22 + Ifges(6,2) * t19 / 0.2e1) * t19 + (t8 * mrSges(4,1) - t9 * mrSges(4,2) + t6 * mrSges(5,2) - t7 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t28) * t28 + (t6 * mrSges(5,1) + t15 * mrSges(4,2) - t8 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t24 + (-Ifges(5,4) + Ifges(4,5)) * t28) * t24 + (t15 * mrSges(4,1) + t7 * mrSges(5,1) - t10 * mrSges(5,2) - t9 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t23 + (Ifges(5,5) - Ifges(4,6)) * t28 + (-Ifges(4,4) - Ifges(5,6)) * t24) * t23;
T = t11;
