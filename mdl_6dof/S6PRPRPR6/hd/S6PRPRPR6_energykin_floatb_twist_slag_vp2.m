% Calculate kinetic energy for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:33
% EndTime: 2019-03-08 19:46:34
% DurationCPUTime: 1.12s
% Computational Cost: add. (3395->160), mult. (5835->224), div. (0->0), fcn. (4850->12), ass. (0->62)
t79 = pkin(2) + pkin(8);
t62 = sin(pkin(10));
t65 = cos(pkin(10));
t53 = t62 * V_base(5) + t65 * V_base(4);
t78 = pkin(7) * t53;
t69 = sin(qJ(2));
t52 = -t62 * V_base(4) + t65 * V_base(5);
t63 = sin(pkin(6));
t66 = cos(pkin(6));
t72 = t52 * t66 + t63 * V_base(6);
t77 = cos(qJ(2));
t42 = t77 * t53 + t72 * t69;
t48 = -t52 * t63 + t66 * V_base(6) + qJD(2);
t58 = V_base(5) * qJ(1) + V_base(1);
t59 = -V_base(4) * qJ(1) + V_base(2);
t50 = t65 * t58 + t62 * t59;
t39 = t72 * pkin(7) + t50;
t49 = -t58 * t62 + t65 * t59;
t43 = V_base(6) * pkin(1) - t66 * t78 + t49;
t60 = V_base(3) + qJD(1);
t46 = -pkin(1) * t52 - t63 * t78 + t60;
t74 = t66 * t77;
t75 = t63 * t77;
t27 = -t69 * t39 + t43 * t74 + t46 * t75;
t71 = qJD(3) - t27;
t19 = t42 * pkin(3) - t79 * t48 + t71;
t41 = -t52 * t74 + t53 * t69 - V_base(6) * t75;
t31 = -t43 * t63 + t66 * t46;
t73 = -qJ(3) * t42 + t31;
t22 = t79 * t41 + t73;
t68 = sin(qJ(4));
t76 = cos(qJ(4));
t12 = t68 * t19 + t76 * t22;
t40 = qJD(4) + t42;
t10 = qJ(5) * t40 + t12;
t28 = t77 * t39 + (t43 * t66 + t46 * t63) * t69;
t26 = -t48 * qJ(3) - t28;
t23 = -pkin(3) * t41 - t26;
t33 = -t76 * t41 + t48 * t68;
t34 = t68 * t41 + t76 * t48;
t15 = pkin(4) * t33 - qJ(5) * t34 + t23;
t61 = sin(pkin(11));
t64 = cos(pkin(11));
t6 = t64 * t10 + t61 * t15;
t5 = -t10 * t61 + t64 * t15;
t11 = t76 * t19 - t68 * t22;
t9 = -t40 * pkin(4) + qJD(5) - t11;
t70 = cos(qJ(6));
t67 = sin(qJ(6));
t32 = qJD(6) + t33;
t30 = t34 * t64 + t40 * t61;
t29 = -t34 * t61 + t40 * t64;
t25 = -t48 * pkin(2) + t71;
t24 = pkin(2) * t41 + t73;
t17 = t67 * t29 + t30 * t70;
t16 = t29 * t70 - t67 * t30;
t7 = -t29 * pkin(5) + t9;
t4 = pkin(9) * t29 + t6;
t3 = pkin(5) * t33 - pkin(9) * t30 + t5;
t2 = t67 * t3 + t4 * t70;
t1 = t3 * t70 - t67 * t4;
t8 = m(2) * (t49 ^ 2 + t50 ^ 2 + t60 ^ 2) / 0.2e1 + m(3) * (t27 ^ 2 + t28 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t23 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t60 * mrSges(2,2) - t49 * mrSges(2,3) + Ifges(2,1) * t53 / 0.2e1) * t53 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t40 / 0.2e1) * t40 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t30 / 0.2e1) * t30 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t60 * mrSges(2,1) + t50 * mrSges(2,3) + Ifges(2,4) * t53 + Ifges(2,2) * t52 / 0.2e1) * t52 + (t23 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t40 + Ifges(5,1) * t34 / 0.2e1) * t34 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t30 + Ifges(6,2) * t29 / 0.2e1) * t29 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t17 / 0.2e1) * t17 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t32 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t27 * mrSges(3,1) - t28 * mrSges(3,2) + t25 * mrSges(4,2) - t26 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t48) * t48 + (t25 * mrSges(4,1) + t31 * mrSges(3,2) - t27 * mrSges(3,3) - t24 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t42 + (-Ifges(4,4) + Ifges(3,5)) * t48) * t42 + (V_base(2) * mrSges(1,1) + t49 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t50 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t53 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t52 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t31 * mrSges(3,1) + t26 * mrSges(4,1) - t24 * mrSges(4,2) - t28 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t41 + (Ifges(4,5) - Ifges(3,6)) * t48 + (-Ifges(3,4) - Ifges(4,6)) * t42) * t41 + (t23 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t34 + Ifges(6,5) * t30 - Ifges(5,6) * t40 + Ifges(6,6) * t29 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t33) * t33;
T  = t8;
