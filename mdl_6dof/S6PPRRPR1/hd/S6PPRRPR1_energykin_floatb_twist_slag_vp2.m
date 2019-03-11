% Calculate kinetic energy for
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:21
% EndTime: 2019-03-08 18:45:23
% DurationCPUTime: 1.60s
% Computational Cost: add. (8555->166), mult. (16221->247), div. (0->0), fcn. (14174->16), ass. (0->71)
t67 = qJ(1) * V_base(5) + V_base(1);
t68 = -qJ(1) * V_base(4) + V_base(2);
t72 = sin(pkin(11));
t77 = cos(pkin(11));
t60 = -t67 * t72 + t77 * t68;
t79 = cos(pkin(6));
t63 = t72 * V_base(5) + t77 * V_base(4);
t88 = qJ(2) * t63;
t54 = V_base(6) * pkin(1) - t79 * t88 + t60;
t62 = -t72 * V_base(4) + t77 * V_base(5);
t69 = V_base(3) + qJD(1);
t74 = sin(pkin(6));
t58 = -pkin(1) * t62 - t74 * t88 + t69;
t97 = t54 * t79 + t58 * t74;
t61 = t77 * t67 + t72 * t68;
t85 = t62 * t79 + t74 * V_base(6);
t51 = t85 * qJ(2) + t61;
t71 = sin(pkin(12));
t76 = cos(pkin(12));
t40 = -t51 * t71 + t97 * t76;
t59 = -t62 * t74 + t79 * V_base(6);
t78 = cos(pkin(7));
t53 = t63 * t76 + t85 * t71;
t95 = pkin(8) * t53;
t33 = pkin(2) * t59 - t78 * t95 + t40;
t45 = -t54 * t74 + t79 * t58 + qJD(2);
t52 = -t63 * t71 + t85 * t76;
t73 = sin(pkin(7));
t39 = -pkin(2) * t52 - t73 * t95 + t45;
t96 = t33 * t78 + t39 * t73;
t41 = t76 * t51 + t97 * t71;
t86 = t52 * t78 + t59 * t73;
t32 = t86 * pkin(8) + t41;
t82 = sin(qJ(3));
t84 = cos(qJ(3));
t24 = -t82 * t32 + t96 * t84;
t43 = -t82 * t53 + t86 * t84;
t25 = t84 * t32 + t96 * t82;
t47 = -t52 * t73 + t59 * t78 + qJD(3);
t19 = pkin(9) * t47 + t25;
t26 = -t33 * t73 + t78 * t39;
t44 = t53 * t84 + t86 * t82;
t23 = -pkin(3) * t43 - pkin(9) * t44 + t26;
t81 = sin(qJ(4));
t94 = cos(qJ(4));
t12 = t94 * t19 + t81 * t23;
t42 = qJD(4) - t43;
t10 = qJ(5) * t42 + t12;
t18 = -t47 * pkin(3) - t24;
t35 = t44 * t81 - t94 * t47;
t36 = t94 * t44 + t81 * t47;
t15 = t35 * pkin(4) - t36 * qJ(5) + t18;
t70 = sin(pkin(13));
t75 = cos(pkin(13));
t6 = t75 * t10 + t70 * t15;
t5 = -t10 * t70 + t75 * t15;
t11 = -t81 * t19 + t94 * t23;
t9 = -t42 * pkin(4) + qJD(5) - t11;
t83 = cos(qJ(6));
t80 = sin(qJ(6));
t34 = qJD(6) + t35;
t28 = t36 * t75 + t42 * t70;
t27 = -t36 * t70 + t42 * t75;
t21 = t27 * t80 + t28 * t83;
t20 = t27 * t83 - t28 * t80;
t7 = -t27 * pkin(5) + t9;
t4 = pkin(10) * t27 + t6;
t3 = pkin(5) * t35 - pkin(10) * t28 + t5;
t2 = t3 * t80 + t4 * t83;
t1 = t3 * t83 - t4 * t80;
t8 = (V_base(2) * mrSges(1,1) + t60 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t61 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t63 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t62 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + m(3) * (t40 ^ 2 + t41 ^ 2 + t45 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t18 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (-t26 * mrSges(4,1) + t25 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,6) * t47 + Ifges(4,2) * t43 / 0.2e1) * t43 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t18 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t36 + Ifges(6,5) * t28 - Ifges(5,6) * t42 + Ifges(6,6) * t27 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t69 * mrSges(2,1) + t61 * mrSges(2,3) + Ifges(2,4) * t63 + Ifges(2,2) * t62 / 0.2e1) * t62 + (t40 * mrSges(3,1) - t41 * mrSges(3,2) + Ifges(3,3) * t59 / 0.2e1) * t59 + (t69 * mrSges(2,2) - t60 * mrSges(2,3) + Ifges(2,1) * t63 / 0.2e1) * t63 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t21 / 0.2e1) * t21 + (t24 * mrSges(4,1) - t25 * mrSges(4,2) + Ifges(4,3) * t47 / 0.2e1) * t47 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t34 + Ifges(7,2) * t20 / 0.2e1) * t20 + (t26 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,5) * t47 + Ifges(4,1) * t44 / 0.2e1) * t44 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t45 * mrSges(3,1) + t41 * mrSges(3,3) + Ifges(3,4) * t53 + Ifges(3,6) * t59 + Ifges(3,2) * t52 / 0.2e1) * t52 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t28 / 0.2e1) * t28 + (t45 * mrSges(3,2) - t40 * mrSges(3,3) + Ifges(3,5) * t59 + Ifges(3,1) * t53 / 0.2e1) * t53 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t60 ^ 2 + t61 ^ 2 + t69 ^ 2) / 0.2e1 + (t18 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t42 + Ifges(5,1) * t36 / 0.2e1) * t36 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t42 / 0.2e1) * t42;
T  = t8;
