% Calculate kinetic energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:20
% EndTime: 2019-03-08 19:21:21
% DurationCPUTime: 1.06s
% Computational Cost: add. (3267->160), mult. (5741->224), div. (0->0), fcn. (4760->12), ass. (0->60)
t64 = sin(pkin(10));
t67 = cos(pkin(10));
t54 = t64 * V_base(5) + t67 * V_base(4);
t79 = pkin(7) * t54;
t70 = sin(qJ(2));
t53 = -t64 * V_base(4) + t67 * V_base(5);
t65 = sin(pkin(6));
t77 = cos(pkin(6));
t74 = t53 * t77 + t65 * V_base(6);
t78 = cos(qJ(2));
t43 = t54 * t78 + t70 * t74;
t49 = t53 * t65 - t77 * V_base(6) - qJD(2);
t59 = qJ(1) * V_base(5) + V_base(1);
t60 = -qJ(1) * V_base(4) + V_base(2);
t51 = t59 * t67 + t60 * t64;
t40 = pkin(7) * t74 + t51;
t50 = -t59 * t64 + t60 * t67;
t44 = V_base(6) * pkin(1) - t77 * t79 + t50;
t62 = V_base(3) + qJD(1);
t47 = -pkin(1) * t53 - t65 * t79 + t62;
t75 = t77 * t78;
t76 = t65 * t78;
t26 = -t40 * t70 + t44 * t75 + t47 * t76;
t73 = qJD(3) - t26;
t16 = -t43 * qJ(4) + (pkin(2) + pkin(3)) * t49 + t73;
t27 = t78 * t40 + (t44 * t77 + t47 * t65) * t70;
t25 = -qJ(3) * t49 + t27;
t42 = -t53 * t75 + t54 * t70 - t76 * V_base(6);
t22 = qJ(4) * t42 + t25;
t63 = sin(pkin(11));
t66 = cos(pkin(11));
t12 = t16 * t63 + t22 * t66;
t10 = pkin(8) * t49 + t12;
t34 = -t44 * t65 + t47 * t77;
t23 = pkin(2) * t42 - qJ(3) * t43 + t34;
t19 = -pkin(3) * t42 + qJD(4) - t23;
t32 = t42 * t66 - t43 * t63;
t33 = t42 * t63 + t43 * t66;
t14 = -pkin(4) * t32 - pkin(8) * t33 + t19;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t6 = t10 * t72 + t14 * t69;
t11 = t16 * t66 - t22 * t63;
t5 = -t10 * t69 + t14 * t72;
t29 = -t33 * t69 + t49 * t72;
t9 = -pkin(4) * t49 - t11;
t71 = cos(qJ(6));
t68 = sin(qJ(6));
t31 = qJD(5) - t32;
t30 = t33 * t72 + t49 * t69;
t28 = qJD(6) - t29;
t24 = pkin(2) * t49 + t73;
t18 = t30 * t71 + t31 * t68;
t17 = -t30 * t68 + t31 * t71;
t7 = -pkin(5) * t29 - pkin(9) * t30 + t9;
t4 = pkin(9) * t31 + t6;
t3 = -pkin(5) * t31 - t5;
t2 = t4 * t71 + t68 * t7;
t1 = -t4 * t68 + t7 * t71;
t8 = m(2) * (t50 ^ 2 + t51 ^ 2 + t62 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t19 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t62 * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,1) * t54 / 0.2e1) * t54 + (t19 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t33 / 0.2e1) * t33 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t31 / 0.2e1) * t31 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t62 * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,2) * t53 / 0.2e1) * t53 + (-t19 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,2) * t32 / 0.2e1) * t32 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t31 + Ifges(6,1) * t30 / 0.2e1) * t30 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t18 / 0.2e1) * t18 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t30 + Ifges(6,6) * t31 + Ifges(6,2) * t29 / 0.2e1) * t29 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t18 + Ifges(7,6) * t28 + Ifges(7,2) * t17 / 0.2e1) * t17 + (t34 * mrSges(3,2) + t24 * mrSges(4,2) - t26 * mrSges(3,3) - t23 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t43) * t43 + (t34 * mrSges(3,1) + t23 * mrSges(4,1) - t25 * mrSges(4,2) - t27 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t42 + (-Ifges(3,4) + Ifges(4,5)) * t43) * t42 + (V_base(2) * mrSges(1,1) + t50 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t51 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t54 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t53 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t26 * mrSges(3,1) + t24 * mrSges(4,1) + t11 * mrSges(5,1) + t27 * mrSges(3,2) - t12 * mrSges(5,2) - t25 * mrSges(4,3) + Ifges(5,5) * t33 + Ifges(5,6) * t32 + (Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t49 + (-Ifges(4,4) - Ifges(3,5)) * t43 + (Ifges(3,6) - Ifges(4,6)) * t42) * t49;
T  = t8;
