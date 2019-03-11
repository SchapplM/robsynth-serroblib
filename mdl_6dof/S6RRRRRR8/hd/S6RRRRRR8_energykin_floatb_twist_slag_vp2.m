% Calculate kinetic energy for
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:11
% EndTime: 2019-03-10 04:48:13
% DurationCPUTime: 1.91s
% Computational Cost: add. (10465->166), mult. (16419->250), div. (0->0), fcn. (14354->16), ass. (0->73)
t67 = pkin(8) * V_base(5) + V_base(1);
t68 = -pkin(8) * V_base(4) + V_base(2);
t79 = sin(qJ(1));
t85 = cos(qJ(1));
t60 = -t67 * t79 + t68 * t85;
t69 = V_base(6) + qJD(1);
t73 = cos(pkin(6));
t63 = t79 * V_base(5) + t85 * V_base(4);
t95 = pkin(9) * t63;
t54 = pkin(1) * t69 - t73 * t95 + t60;
t62 = -t79 * V_base(4) + t85 * V_base(5);
t71 = sin(pkin(6));
t58 = -pkin(1) * t62 - t71 * t95 + V_base(3);
t97 = t54 * t73 + t58 * t71;
t61 = t67 * t85 + t68 * t79;
t87 = t62 * t73 + t69 * t71;
t51 = pkin(9) * t87 + t61;
t78 = sin(qJ(2));
t84 = cos(qJ(2));
t39 = -t51 * t78 + t84 * t97;
t59 = -t62 * t71 + t69 * t73 + qJD(2);
t72 = cos(pkin(7));
t53 = t63 * t84 + t78 * t87;
t94 = pkin(10) * t53;
t35 = pkin(2) * t59 - t72 * t94 + t39;
t45 = -t54 * t71 + t58 * t73;
t52 = -t63 * t78 + t84 * t87;
t70 = sin(pkin(7));
t38 = -pkin(2) * t52 - t70 * t94 + t45;
t96 = t35 * t72 + t38 * t70;
t40 = t84 * t51 + t78 * t97;
t88 = t52 * t72 + t59 * t70;
t34 = pkin(10) * t88 + t40;
t77 = sin(qJ(3));
t83 = cos(qJ(3));
t21 = -t34 * t77 + t83 * t96;
t43 = -t53 * t77 + t83 * t88;
t22 = t83 * t34 + t77 * t96;
t46 = -t52 * t70 + t59 * t72 + qJD(3);
t17 = pkin(11) * t46 + t22;
t28 = -t35 * t70 + t38 * t72;
t44 = t53 * t83 + t77 * t88;
t20 = -pkin(3) * t43 - pkin(11) * t44 + t28;
t76 = sin(qJ(4));
t82 = cos(qJ(4));
t13 = t17 * t82 + t20 * t76;
t32 = -t44 * t76 + t46 * t82;
t11 = pkin(12) * t32 + t13;
t75 = sin(qJ(5));
t81 = cos(qJ(5));
t12 = -t17 * t76 + t20 * t82;
t33 = t44 * t82 + t46 * t76;
t42 = qJD(4) - t43;
t9 = pkin(4) * t42 - pkin(12) * t33 + t12;
t6 = t11 * t81 + t75 * t9;
t5 = -t11 * t75 + t81 * t9;
t26 = t32 * t81 - t33 * t75;
t16 = -pkin(3) * t46 - t21;
t14 = -pkin(4) * t32 + t16;
t86 = V_base(3) ^ 2;
t80 = cos(qJ(6));
t74 = sin(qJ(6));
t41 = qJD(5) + t42;
t27 = t32 * t75 + t33 * t81;
t25 = qJD(6) - t26;
t24 = t27 * t80 + t41 * t74;
t23 = -t27 * t74 + t41 * t80;
t7 = -pkin(5) * t26 - pkin(13) * t27 + t14;
t4 = pkin(13) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t80 + t7 * t74;
t1 = -t4 * t74 + t7 * t80;
t8 = (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t41 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t24 / 0.2e1) * t24 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t24 + Ifges(7,6) * t25 + Ifges(7,2) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t28 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,6) * t46 + Ifges(4,2) * t43 / 0.2e1) * t43 + (t45 * mrSges(3,2) - t39 * mrSges(3,3) + Ifges(3,5) * t59 + Ifges(3,1) * t53 / 0.2e1) * t53 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t16 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t42 + Ifges(5,1) * t33 / 0.2e1) * t33 + (t28 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t44 / 0.2e1) * t44 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t41 / 0.2e1) * t41 + (-V_base(3) * mrSges(2,1) + t61 * mrSges(2,3) + Ifges(2,4) * t63 + Ifges(2,6) * t69 + Ifges(2,2) * t62 / 0.2e1) * t62 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t41 + Ifges(6,1) * t27 / 0.2e1) * t27 + (-t16 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,6) * t42 + Ifges(5,2) * t32 / 0.2e1) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (t60 * mrSges(2,1) - t61 * mrSges(2,2) + Ifges(2,3) * t69 / 0.2e1) * t69 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(3) * mrSges(2,2) - t60 * mrSges(2,3) + Ifges(2,5) * t69 + Ifges(2,1) * t63 / 0.2e1) * t63 + (-t45 * mrSges(3,1) + t40 * mrSges(3,3) + Ifges(3,4) * t53 + Ifges(3,6) * t59 + Ifges(3,2) * t52 / 0.2e1) * t52 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t28 ^ 2) / 0.2e1 + m(3) * (t39 ^ 2 + t40 ^ 2 + t45 ^ 2) / 0.2e1 + (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t42 / 0.2e1) * t42 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t86) / 0.2e1 + m(2) * (t60 ^ 2 + t61 ^ 2 + t86) / 0.2e1 + (t39 * mrSges(3,1) - t40 * mrSges(3,2) + Ifges(3,3) * t59 / 0.2e1) * t59;
T  = t8;
