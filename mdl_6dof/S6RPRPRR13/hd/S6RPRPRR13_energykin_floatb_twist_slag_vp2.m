% Calculate kinetic energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR13_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:30
% EndTime: 2019-03-09 04:21:32
% DurationCPUTime: 1.44s
% Computational Cost: add. (6251->165), mult. (10197->237), div. (0->0), fcn. (8792->14), ass. (0->70)
t61 = V_base(5) * pkin(8) + V_base(1);
t62 = -V_base(4) * pkin(8) + V_base(2);
t73 = sin(qJ(1));
t76 = cos(qJ(1));
t54 = -t61 * t73 + t76 * t62;
t63 = V_base(6) + qJD(1);
t69 = cos(pkin(6));
t57 = t73 * V_base(5) + t76 * V_base(4);
t84 = qJ(2) * t57;
t47 = pkin(1) * t63 - t69 * t84 + t54;
t56 = -t73 * V_base(4) + t76 * V_base(5);
t66 = sin(pkin(6));
t51 = -pkin(1) * t56 - t66 * t84 + V_base(3);
t90 = t47 * t69 + t51 * t66;
t89 = pkin(3) + pkin(10);
t64 = sin(pkin(12));
t67 = cos(pkin(12));
t80 = t56 * t69 + t63 * t66;
t46 = t57 * t67 + t64 * t80;
t88 = pkin(9) * t46;
t45 = -t57 * t64 + t67 * t80;
t53 = -t56 * t66 + t63 * t69;
t72 = sin(qJ(3));
t68 = cos(pkin(7));
t87 = cos(qJ(3));
t82 = t68 * t87;
t65 = sin(pkin(7));
t83 = t65 * t87;
t35 = -t45 * t82 + t46 * t72 - t53 * t83;
t55 = t76 * t61 + t73 * t62;
t44 = qJ(2) * t80 + t55;
t32 = -t44 * t64 + t90 * t67;
t28 = pkin(2) * t53 - t68 * t88 + t32;
t37 = -t47 * t66 + t69 * t51 + qJD(2);
t31 = -pkin(2) * t45 - t65 * t88 + t37;
t18 = -t28 * t65 + t68 * t31;
t81 = t45 * t68 + t53 * t65;
t36 = t46 * t87 + t72 * t81;
t79 = -qJ(4) * t36 + t18;
t12 = t35 * t89 + t79;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t39 = -t45 * t65 + t53 * t68 + qJD(3);
t33 = t67 * t44 + t90 * t64;
t25 = pkin(9) * t81 + t33;
t16 = -t72 * t25 + t28 * t82 + t31 * t83;
t78 = qJD(4) - t16;
t9 = t36 * pkin(4) - t39 * t89 + t78;
t6 = t75 * t12 + t71 * t9;
t17 = t87 * t25 + (t28 * t68 + t31 * t65) * t72;
t14 = -t39 * qJ(4) - t17;
t5 = -t12 * t71 + t75 * t9;
t26 = t35 * t75 - t39 * t71;
t10 = -pkin(4) * t35 - t14;
t77 = V_base(3) ^ 2;
t74 = cos(qJ(6));
t70 = sin(qJ(6));
t34 = qJD(5) + t36;
t27 = t35 * t71 + t39 * t75;
t24 = qJD(6) - t26;
t20 = t27 * t74 + t34 * t70;
t19 = -t27 * t70 + t34 * t74;
t15 = pkin(3) * t35 + t79;
t13 = -t39 * pkin(3) + t78;
t7 = -pkin(5) * t26 - pkin(11) * t27 + t10;
t4 = pkin(11) * t34 + t6;
t3 = -pkin(5) * t34 - t5;
t2 = t4 * t74 + t7 * t70;
t1 = -t4 * t70 + t7 * t74;
t8 = m(3) * (t32 ^ 2 + t33 ^ 2 + t37 ^ 2) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (t37 * mrSges(3,2) - t32 * mrSges(3,3) + Ifges(3,5) * t53 + Ifges(3,1) * t46 / 0.2e1) * t46 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (-t10 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t34 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t18 * mrSges(4,1) + t14 * mrSges(5,1) - t15 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t35 + (Ifges(5,5) - Ifges(4,6)) * t39 + (-Ifges(4,4) - Ifges(5,6)) * t36) * t35 + (-t37 * mrSges(3,1) + t33 * mrSges(3,3) + Ifges(3,4) * t46 + Ifges(3,6) * t53 + Ifges(3,2) * t45 / 0.2e1) * t45 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t24 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t13 * mrSges(5,1) + t18 * mrSges(4,2) - t16 * mrSges(4,3) - t15 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t36 + (-Ifges(5,4) + Ifges(4,5)) * t39) * t36 + (t54 * mrSges(2,1) - t55 * mrSges(2,2) + Ifges(2,3) * t63 / 0.2e1) * t63 + (t10 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(2,2) - t54 * mrSges(2,3) + Ifges(2,5) * t63 + Ifges(2,1) * t57 / 0.2e1) * t57 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + t13 * mrSges(5,2) - t14 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t39) * t39 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t77) / 0.2e1 + m(2) * (t54 ^ 2 + t55 ^ 2 + t77) / 0.2e1 + (t32 * mrSges(3,1) - t33 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + (-V_base(3) * mrSges(2,1) + t55 * mrSges(2,3) + Ifges(2,4) * t57 + Ifges(2,6) * t63 + Ifges(2,2) * t56 / 0.2e1) * t56 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t20 / 0.2e1) * t20;
T  = t8;
