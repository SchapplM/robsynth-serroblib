% Calculate kinetic energy for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:38:42
% EndTime: 2019-03-08 23:38:43
% DurationCPUTime: 1.65s
% Computational Cost: add. (8953->166), mult. (16221->247), div. (0->0), fcn. (14174->16), ass. (0->71)
t67 = qJ(1) * V_base(5) + V_base(1);
t68 = -qJ(1) * V_base(4) + V_base(2);
t71 = sin(pkin(12));
t75 = cos(pkin(12));
t60 = -t67 * t71 + t75 * t68;
t77 = cos(pkin(6));
t63 = t71 * V_base(5) + t75 * V_base(4);
t95 = pkin(8) * t63;
t54 = V_base(6) * pkin(1) - t77 * t95 + t60;
t62 = -t71 * V_base(4) + t75 * V_base(5);
t69 = V_base(3) + qJD(1);
t73 = sin(pkin(6));
t58 = -pkin(1) * t62 - t73 * t95 + t69;
t97 = t54 * t77 + t58 * t73;
t61 = t75 * t67 + t71 * t68;
t85 = t62 * t77 + t73 * V_base(6);
t51 = t85 * pkin(8) + t61;
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t40 = -t51 * t81 + t97 * t84;
t59 = -t62 * t73 + t77 * V_base(6) + qJD(2);
t76 = cos(pkin(7));
t53 = t63 * t84 + t85 * t81;
t94 = pkin(9) * t53;
t34 = pkin(2) * t59 - t76 * t94 + t40;
t45 = -t54 * t73 + t77 * t58;
t52 = -t81 * t63 + t85 * t84;
t72 = sin(pkin(7));
t39 = -pkin(2) * t52 - t72 * t94 + t45;
t96 = t34 * t76 + t39 * t72;
t41 = t84 * t51 + t97 * t81;
t86 = t52 * t76 + t59 * t72;
t32 = t86 * pkin(9) + t41;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t24 = -t80 * t32 + t96 * t83;
t43 = -t53 * t80 + t86 * t83;
t25 = t83 * t32 + t96 * t80;
t47 = -t52 * t72 + t59 * t76 + qJD(3);
t21 = pkin(10) * t47 + t25;
t26 = -t34 * t72 + t76 * t39;
t44 = t53 * t83 + t86 * t80;
t23 = -pkin(3) * t43 - pkin(10) * t44 + t26;
t79 = sin(qJ(4));
t93 = cos(qJ(4));
t12 = t93 * t21 + t79 * t23;
t42 = qJD(4) - t43;
t10 = qJ(5) * t42 + t12;
t20 = -pkin(3) * t47 - t24;
t35 = t44 * t79 - t93 * t47;
t36 = t93 * t44 + t79 * t47;
t15 = pkin(4) * t35 - qJ(5) * t36 + t20;
t70 = sin(pkin(13));
t74 = cos(pkin(13));
t6 = t74 * t10 + t70 * t15;
t5 = -t10 * t70 + t74 * t15;
t11 = -t79 * t21 + t93 * t23;
t9 = -t42 * pkin(4) + qJD(5) - t11;
t82 = cos(qJ(6));
t78 = sin(qJ(6));
t33 = qJD(6) + t35;
t28 = t36 * t74 + t42 * t70;
t27 = -t36 * t70 + t42 * t74;
t19 = t27 * t78 + t28 * t82;
t18 = t27 * t82 - t28 * t78;
t7 = -t27 * pkin(5) + t9;
t4 = pkin(11) * t27 + t6;
t3 = pkin(5) * t35 - pkin(11) * t28 + t5;
t2 = t3 * t78 + t4 * t82;
t1 = t3 * t82 - t4 * t78;
t8 = m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t33 / 0.2e1) * t33 + (-t26 * mrSges(4,1) + t25 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,6) * t47 + Ifges(4,2) * t43 / 0.2e1) * t43 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t28 / 0.2e1) * t28 + (t20 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t36 + Ifges(6,5) * t28 - Ifges(5,6) * t42 + Ifges(6,6) * t27 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (-t69 * mrSges(2,1) + t61 * mrSges(2,3) + Ifges(2,4) * t63 + Ifges(2,2) * t62 / 0.2e1) * t62 + (t20 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t42 + Ifges(5,1) * t36 / 0.2e1) * t36 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t33 + Ifges(7,1) * t19 / 0.2e1) * t19 + (t69 * mrSges(2,2) - t60 * mrSges(2,3) + Ifges(2,1) * t63 / 0.2e1) * t63 + (t26 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,5) * t47 + Ifges(4,1) * t44 / 0.2e1) * t44 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t33 + Ifges(7,2) * t18 / 0.2e1) * t18 + (-t45 * mrSges(3,1) + t41 * mrSges(3,3) + Ifges(3,4) * t53 + Ifges(3,6) * t59 + Ifges(3,2) * t52 / 0.2e1) * t52 + (t45 * mrSges(3,2) - t40 * mrSges(3,3) + Ifges(3,5) * t59 + Ifges(3,1) * t53 / 0.2e1) * t53 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,2) * t27 / 0.2e1) * t27 + (t40 * mrSges(3,1) - t41 * mrSges(3,2) + Ifges(3,3) * t59 / 0.2e1) * t59 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t42 / 0.2e1) * t42 + (V_base(2) * mrSges(1,1) + t60 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t61 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t63 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t62 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t24 * mrSges(4,1) - t25 * mrSges(4,2) + Ifges(4,3) * t47 / 0.2e1) * t47 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(3) * (t40 ^ 2 + t41 ^ 2 + t45 ^ 2) / 0.2e1 + m(2) * (t60 ^ 2 + t61 ^ 2 + t69 ^ 2) / 0.2e1;
T  = t8;
