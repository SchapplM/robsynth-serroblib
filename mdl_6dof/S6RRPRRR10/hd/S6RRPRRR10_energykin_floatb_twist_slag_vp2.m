% Calculate kinetic energy for
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:03
% EndTime: 2019-03-09 14:18:05
% DurationCPUTime: 1.53s
% Computational Cost: add. (6459->161), mult. (9743->238), div. (0->0), fcn. (8340->14), ass. (0->66)
t76 = sin(qJ(1));
t81 = cos(qJ(1));
t59 = -t76 * V_base(4) + t81 * V_base(5);
t67 = V_base(6) + qJD(1);
t69 = sin(pkin(6));
t71 = cos(pkin(6));
t83 = t59 * t71 + t67 * t69;
t65 = V_base(5) * pkin(7) + V_base(1);
t66 = -V_base(4) * pkin(7) + V_base(2);
t56 = -t65 * t76 + t81 * t66;
t60 = t76 * V_base(5) + t81 * V_base(4);
t89 = pkin(8) * t60;
t51 = pkin(1) * t67 - t71 * t89 + t56;
t54 = -pkin(1) * t59 - t69 * t89 + V_base(3);
t90 = t51 * t71 + t54 * t69;
t57 = t81 * t65 + t76 * t66;
t48 = t83 * pkin(8) + t57;
t75 = sin(qJ(2));
t80 = cos(qJ(2));
t38 = -t75 * t48 + t90 * t80;
t40 = -t51 * t69 + t71 * t54;
t49 = t60 * t75 - t83 * t80;
t50 = t60 * t80 + t83 * t75;
t30 = pkin(2) * t49 - qJ(3) * t50 + t40;
t39 = t80 * t48 + t90 * t75;
t55 = -t59 * t69 + t67 * t71 + qJD(2);
t35 = qJ(3) * t55 + t39;
t68 = sin(pkin(12));
t70 = cos(pkin(12));
t23 = t70 * t30 - t35 * t68;
t43 = t50 * t70 + t55 * t68;
t17 = pkin(3) * t49 - pkin(9) * t43 + t23;
t24 = t68 * t30 + t70 * t35;
t42 = -t50 * t68 + t55 * t70;
t22 = pkin(9) * t42 + t24;
t74 = sin(qJ(4));
t79 = cos(qJ(4));
t12 = t74 * t17 + t79 * t22;
t47 = qJD(4) + t49;
t10 = pkin(10) * t47 + t12;
t33 = -pkin(2) * t55 + qJD(3) - t38;
t25 = -pkin(3) * t42 + t33;
t36 = t42 * t79 - t74 * t43;
t37 = t42 * t74 + t43 * t79;
t15 = -pkin(4) * t36 - pkin(10) * t37 + t25;
t73 = sin(qJ(5));
t78 = cos(qJ(5));
t6 = t78 * t10 + t73 * t15;
t5 = -t10 * t73 + t78 * t15;
t11 = t17 * t79 - t74 * t22;
t34 = qJD(5) - t36;
t9 = -pkin(4) * t47 - t11;
t82 = V_base(3) ^ 2;
t77 = cos(qJ(6));
t72 = sin(qJ(6));
t32 = qJD(6) + t34;
t27 = t37 * t78 + t47 * t73;
t26 = -t37 * t73 + t47 * t78;
t19 = t26 * t72 + t27 * t77;
t18 = t26 * t77 - t27 * t72;
t7 = -pkin(5) * t26 + t9;
t4 = pkin(11) * t26 + t6;
t3 = pkin(5) * t34 - pkin(11) * t27 + t5;
t2 = t3 * t72 + t4 * t77;
t1 = t3 * t77 - t4 * t72;
t8 = (-V_base(3) * mrSges(2,1) + t57 * mrSges(2,3) + Ifges(2,4) * t60 + Ifges(2,6) * t67 + Ifges(2,2) * t59 / 0.2e1) * t59 + (-t25 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,6) * t47 + Ifges(5,2) * t36 / 0.2e1) * t36 + (t40 * mrSges(3,1) + t23 * mrSges(4,1) - t24 * mrSges(4,2) - t39 * mrSges(3,3) - Ifges(3,4) * t50 + Ifges(4,5) * t43 - Ifges(3,6) * t55 + Ifges(4,6) * t42 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t49) * t49 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (V_base(3) * mrSges(2,2) - t56 * mrSges(2,3) + Ifges(2,5) * t67 + Ifges(2,1) * t60 / 0.2e1) * t60 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t34 + Ifges(6,2) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t19 / 0.2e1) * t19 + (t56 * mrSges(2,1) - t57 * mrSges(2,2) + Ifges(2,3) * t67 / 0.2e1) * t67 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t40 * mrSges(3,2) - t38 * mrSges(3,3) + Ifges(3,5) * t55 + Ifges(3,1) * t50 / 0.2e1) * t50 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t32 + Ifges(7,2) * t18 / 0.2e1) * t18 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t25 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t47 + Ifges(5,1) * t37 / 0.2e1) * t37 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (-t33 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,2) * t42 / 0.2e1) * t42 + (t33 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,1) * t43 / 0.2e1) * t43 + (t38 * mrSges(3,1) - t39 * mrSges(3,2) + Ifges(3,3) * t55 / 0.2e1) * t55 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t82) / 0.2e1 + m(2) * (t56 ^ 2 + t57 ^ 2 + t82) / 0.2e1;
T  = t8;
