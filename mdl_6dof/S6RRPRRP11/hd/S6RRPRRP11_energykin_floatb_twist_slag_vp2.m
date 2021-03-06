% Calculate kinetic energy for
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP11_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:18
% EndTime: 2019-03-09 12:45:19
% DurationCPUTime: 0.92s
% Computational Cost: add. (2017->149), mult. (2525->196), div. (0->0), fcn. (1896->8), ass. (0->50)
t63 = pkin(2) + pkin(8);
t52 = sin(qJ(5));
t56 = cos(qJ(5));
t55 = sin(qJ(1));
t59 = cos(qJ(1));
t44 = t55 * V_base(5) + t59 * V_base(4);
t51 = V_base(6) + qJD(1);
t54 = sin(qJ(2));
t58 = cos(qJ(2));
t38 = t44 * t58 + t51 * t54;
t43 = -t55 * V_base(4) + t59 * V_base(5);
t42 = qJD(2) - t43;
t30 = -pkin(1) * t43 - pkin(7) * t44 + V_base(3);
t49 = V_base(5) * pkin(6) + V_base(1);
t50 = -V_base(4) * pkin(6) + V_base(2);
t40 = t59 * t49 + t55 * t50;
t35 = pkin(7) * t51 + t40;
t25 = t30 * t58 - t54 * t35;
t62 = qJD(3) - t25;
t16 = pkin(3) * t38 - t42 * t63 + t62;
t37 = t44 * t54 - t58 * t51;
t39 = -t55 * t49 + t50 * t59;
t34 = -pkin(1) * t51 - t39;
t61 = -qJ(3) * t38 + t34;
t19 = t37 * t63 + t61;
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t11 = t57 * t16 - t19 * t53;
t28 = t37 * t53 + t42 * t57;
t36 = qJD(4) + t38;
t7 = pkin(4) * t36 - pkin(9) * t28 + t11;
t12 = t53 * t16 + t57 * t19;
t27 = t37 * t57 - t42 * t53;
t9 = pkin(9) * t27 + t12;
t4 = t52 * t7 + t56 * t9;
t26 = t54 * t30 + t58 * t35;
t23 = -t42 * qJ(3) - t26;
t3 = -t52 * t9 + t56 * t7;
t18 = -pkin(3) * t37 - t23;
t13 = -pkin(4) * t27 + t18;
t60 = V_base(3) ^ 2;
t33 = qJD(5) + t36;
t24 = pkin(2) * t37 + t61;
t22 = -pkin(2) * t42 + t62;
t21 = t27 * t52 + t28 * t56;
t20 = t27 * t56 - t28 * t52;
t10 = -pkin(5) * t20 + qJD(6) + t13;
t2 = qJ(6) * t20 + t4;
t1 = pkin(5) * t33 - qJ(6) * t21 + t3;
t5 = m(2) * (t39 ^ 2 + t40 ^ 2 + t60) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t60) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t18 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t51 / 0.2e1) * t51 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t36 / 0.2e1) * t36 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t44 / 0.2e1) * t44 + (t18 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t36 + Ifges(5,1) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t51 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t18 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t28 + Ifges(5,6) * t36 + Ifges(5,2) * t27 / 0.2e1) * t27 + (t25 * mrSges(3,1) - t26 * mrSges(3,2) + t22 * mrSges(4,2) - t23 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t42) * t42 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t33) * t33 + (t22 * mrSges(4,1) + t34 * mrSges(3,2) - t25 * mrSges(3,3) - t24 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t38 + (-Ifges(4,4) + Ifges(3,5)) * t42) * t38 + (t13 * mrSges(6,2) + t10 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t21 + (Ifges(6,5) + Ifges(7,5)) * t33) * t21 + (t34 * mrSges(3,1) + t23 * mrSges(4,1) - t24 * mrSges(4,2) - t26 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t37 + (Ifges(4,5) - Ifges(3,6)) * t42 + (-Ifges(3,4) - Ifges(4,6)) * t38) * t37 + (-t13 * mrSges(6,1) - t10 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t20 + (Ifges(6,6) + Ifges(7,6)) * t33 + (Ifges(6,4) + Ifges(7,4)) * t21) * t20;
T  = t5;
