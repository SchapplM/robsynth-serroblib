% Calculate kinetic energy for
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:27
% EndTime: 2019-03-09 05:34:28
% DurationCPUTime: 1.01s
% Computational Cost: add. (1623->150), mult. (2011->196), div. (0->0), fcn. (1416->8), ass. (0->53)
t68 = pkin(1) + pkin(7);
t67 = -pkin(4) - pkin(5);
t66 = cos(qJ(1));
t65 = cos(qJ(4));
t56 = sin(qJ(1));
t44 = t56 * V_base(5) + t66 * V_base(4);
t52 = V_base(6) + qJD(1);
t48 = V_base(5) * pkin(6) + V_base(1);
t49 = -V_base(4) * pkin(6) + V_base(2);
t38 = -t56 * t48 + t49 * t66;
t62 = qJD(2) - t38;
t23 = t44 * pkin(2) - t52 * t68 + t62;
t43 = t56 * V_base(4) - t66 * V_base(5);
t63 = -qJ(2) * t44 + V_base(3);
t28 = t43 * t68 + t63;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t16 = t55 * t23 + t58 * t28;
t42 = qJD(3) + t44;
t14 = pkin(8) * t42 + t16;
t39 = t66 * t48 + t56 * t49;
t33 = -t52 * qJ(2) - t39;
t29 = -pkin(2) * t43 - t33;
t36 = t58 * t43 - t52 * t55;
t37 = t43 * t55 + t52 * t58;
t20 = -pkin(3) * t36 - pkin(8) * t37 + t29;
t54 = sin(qJ(4));
t9 = t65 * t14 + t54 * t20;
t15 = t58 * t23 - t55 * t28;
t35 = qJD(4) - t36;
t7 = t35 * qJ(5) + t9;
t64 = pkin(3) * t42 + t15;
t8 = -t54 * t14 + t20 * t65;
t61 = qJD(5) - t8;
t27 = t37 * t65 + t54 * t42;
t60 = qJ(5) * t27 + t64;
t59 = V_base(3) ^ 2;
t57 = cos(qJ(6));
t53 = sin(qJ(6));
t34 = qJD(6) - t35;
t31 = -t52 * pkin(1) + t62;
t30 = pkin(1) * t43 + t63;
t26 = t37 * t54 - t42 * t65;
t18 = t26 * t53 + t27 * t57;
t17 = t26 * t57 - t27 * t53;
t10 = pkin(4) * t26 - t60;
t6 = -t35 * pkin(4) + t61;
t5 = t26 * t67 + t60;
t4 = pkin(9) * t26 + t7;
t3 = -t27 * pkin(9) + t35 * t67 + t61;
t2 = t3 * t53 + t4 * t57;
t1 = t3 * t57 - t4 * t53;
t11 = m(2) * (t38 ^ 2 + t39 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t30 ^ 2 + t31 ^ 2 + t33 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t42 / 0.2e1) * t42 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t29 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t42 + Ifges(4,1) * t37 / 0.2e1) * t37 + (t5 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t18 / 0.2e1) * t18 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t29 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,6) * t42 + Ifges(4,2) * t36 / 0.2e1) * t36 + (-t5 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t18 + Ifges(7,6) * t34 + Ifges(7,2) * t17 / 0.2e1) * t17 + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + t31 * mrSges(3,2) - t33 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t52) * t52 + (t8 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) + t7 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t35) * t35 + (t31 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) - t30 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t44 + (-Ifges(3,4) + Ifges(2,5)) * t52) * t44 + (-t64 * mrSges(5,2) + t6 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t27 + (Ifges(6,4) + Ifges(5,5)) * t35) * t27 + (V_base(3) * mrSges(2,1) + t33 * mrSges(3,1) - t30 * mrSges(3,2) - t39 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t43 + (Ifges(3,5) - Ifges(2,6)) * t52 + (-Ifges(2,4) - Ifges(3,6)) * t44) * t43 + (-t64 * mrSges(5,1) + t10 * mrSges(6,1) - t7 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t26 + (-Ifges(5,6) + Ifges(6,6)) * t35 + (-Ifges(5,4) + Ifges(6,5)) * t27) * t26;
T  = t11;
