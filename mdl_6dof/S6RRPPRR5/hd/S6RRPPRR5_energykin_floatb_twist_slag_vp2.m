% Calculate kinetic energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:17
% EndTime: 2019-03-09 09:07:18
% DurationCPUTime: 1.07s
% Computational Cost: add. (2471->156), mult. (3687->209), div. (0->0), fcn. (2958->10), ass. (0->56)
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t53 = t65 * V_base(5) + t68 * V_base(4);
t78 = pkin(8) * t53;
t59 = V_base(5) * pkin(7) + V_base(1);
t60 = -V_base(4) * pkin(7) + V_base(2);
t50 = t68 * t59 + t65 * t60;
t52 = -t65 * V_base(4) + t68 * V_base(5);
t61 = sin(pkin(6));
t74 = V_base(6) + qJD(1);
t71 = t61 * t74;
t75 = cos(pkin(6));
t70 = t52 * t75 + t71;
t37 = pkin(8) * t70 + t50;
t49 = -t65 * t59 + t68 * t60;
t40 = pkin(1) * t74 - t75 * t78 + t49;
t64 = sin(qJ(2));
t44 = -pkin(1) * t52 - t61 * t78 + V_base(3);
t76 = t44 * t61;
t77 = cos(qJ(2));
t20 = t77 * t37 + (t40 * t75 + t76) * t64;
t48 = t52 * t61 - t74 * t75 - qJD(2);
t18 = -t48 * qJ(3) + t20;
t73 = t75 * t77;
t38 = -t52 * t73 + t53 * t64 - t77 * t71;
t15 = t38 * qJ(4) + t18;
t13 = pkin(9) * t48 + t15;
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t39 = t53 * t77 + t64 * t70;
t23 = -t61 * t40 + t75 * t44;
t16 = t38 * pkin(2) - t39 * qJ(3) + t23;
t72 = qJD(4) - t16;
t9 = -pkin(9) * t39 + (-pkin(3) - pkin(4)) * t38 + t72;
t6 = t67 * t13 + t63 * t9;
t19 = -t64 * t37 + t40 * t73 + t77 * t76;
t5 = -t13 * t63 + t67 * t9;
t17 = t48 * pkin(2) + qJD(3) - t19;
t25 = -t39 * t63 + t48 * t67;
t12 = t48 * pkin(3) - t39 * qJ(4) + t17;
t10 = -pkin(4) * t48 - t12;
t69 = V_base(3) ^ 2;
t66 = cos(qJ(6));
t62 = sin(qJ(6));
t36 = qJD(5) - t38;
t26 = t39 * t67 + t48 * t63;
t24 = qJD(6) - t25;
t22 = t26 * t66 + t36 * t62;
t21 = -t26 * t62 + t36 * t66;
t14 = -pkin(3) * t38 + t72;
t7 = -pkin(5) * t25 - pkin(10) * t26 + t10;
t4 = pkin(10) * t36 + t6;
t3 = -pkin(5) * t36 - t5;
t2 = t4 * t66 + t62 * t7;
t1 = -t4 * t62 + t66 * t7;
t8 = (t10 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t36 + Ifges(6,1) * t26 / 0.2e1) * t26 + (-t10 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t36 + Ifges(6,2) * t25 / 0.2e1) * t25 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t22 / 0.2e1) * t22 + (t50 * mrSges(2,3) - V_base(3) * mrSges(2,1) + Ifges(2,4) * t53 + Ifges(2,2) * t52 / 0.2e1) * t52 + (-t49 * mrSges(2,3) + V_base(3) * mrSges(2,2) + Ifges(2,1) * t53 / 0.2e1) * t53 + m(2) * (t49 ^ 2 + t50 ^ 2 + t69) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t69) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(3) * (t19 ^ 2 + t20 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t22 + Ifges(7,6) * t24 + Ifges(7,2) * t21 / 0.2e1) * t21 + (t23 * mrSges(3,1) + t16 * mrSges(4,1) - t14 * mrSges(5,1) - t18 * mrSges(4,2) - t20 * mrSges(3,3) + t15 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t38 + (Ifges(3,6) - Ifges(4,6) + Ifges(5,6)) * t48 + (-Ifges(3,4) + Ifges(5,4) + Ifges(4,5)) * t39) * t38 + (t23 * mrSges(3,2) + t17 * mrSges(4,2) + t14 * mrSges(5,2) - t19 * mrSges(3,3) - t16 * mrSges(4,3) - t12 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t39 + (-Ifges(4,4) - Ifges(3,5) + Ifges(5,5)) * t48) * t39 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t19 * mrSges(3,1) + t17 * mrSges(4,1) + t12 * mrSges(5,1) + t20 * mrSges(3,2) - t15 * mrSges(5,2) - t18 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t48) * t48 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t36 / 0.2e1) * t36 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t49 * mrSges(2,1) - t50 * mrSges(2,2) + Ifges(2,5) * t53 + Ifges(2,6) * t52 + Ifges(2,3) * t74 / 0.2e1) * t74;
T  = t8;
