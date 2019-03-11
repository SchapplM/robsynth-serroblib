% Calculate kinetic energy for
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:19
% EndTime: 2019-03-09 00:42:20
% DurationCPUTime: 1.54s
% Computational Cost: add. (5685->161), mult. (9743->239), div. (0->0), fcn. (8340->14), ass. (0->66)
t65 = V_base(5) * qJ(1) + V_base(1);
t66 = -V_base(4) * qJ(1) + V_base(2);
t68 = sin(pkin(12));
t70 = cos(pkin(12));
t57 = -t65 * t68 + t70 * t66;
t71 = cos(pkin(6));
t61 = t68 * V_base(5) + t70 * V_base(4);
t86 = pkin(7) * t61;
t52 = V_base(6) * pkin(1) - t71 * t86 + t57;
t60 = -t68 * V_base(4) + t70 * V_base(5);
t67 = V_base(3) + qJD(1);
t69 = sin(pkin(6));
t55 = -pkin(1) * t60 - t69 * t86 + t67;
t87 = t52 * t71 + t55 * t69;
t58 = t70 * t65 + t68 * t66;
t82 = t60 * t71 + t69 * V_base(6);
t48 = t82 * pkin(7) + t58;
t76 = sin(qJ(2));
t81 = cos(qJ(2));
t38 = -t76 * t48 + t87 * t81;
t50 = -t76 * t61 + t82 * t81;
t40 = -t52 * t69 + t71 * t55;
t51 = t61 * t81 + t82 * t76;
t30 = -pkin(2) * t50 - pkin(8) * t51 + t40;
t39 = t81 * t48 + t87 * t76;
t56 = -t60 * t69 + t71 * V_base(6) + qJD(2);
t33 = pkin(8) * t56 + t39;
t75 = sin(qJ(3));
t80 = cos(qJ(3));
t23 = t80 * t30 - t33 * t75;
t43 = t51 * t80 + t56 * t75;
t49 = qJD(3) - t50;
t17 = pkin(3) * t49 - pkin(9) * t43 + t23;
t24 = t75 * t30 + t80 * t33;
t42 = -t51 * t75 + t56 * t80;
t22 = pkin(9) * t42 + t24;
t74 = sin(qJ(4));
t79 = cos(qJ(4));
t12 = t74 * t17 + t79 * t22;
t47 = qJD(4) + t49;
t10 = pkin(10) * t47 + t12;
t32 = -t56 * pkin(2) - t38;
t25 = -t42 * pkin(3) + t32;
t36 = t42 * t79 - t74 * t43;
t37 = t42 * t74 + t43 * t79;
t15 = -t36 * pkin(4) - t37 * pkin(10) + t25;
t73 = sin(qJ(5));
t78 = cos(qJ(5));
t6 = t78 * t10 + t73 * t15;
t5 = -t10 * t73 + t78 * t15;
t11 = t17 * t79 - t74 * t22;
t35 = qJD(5) - t36;
t9 = -pkin(4) * t47 - t11;
t77 = cos(qJ(6));
t72 = sin(qJ(6));
t34 = qJD(6) + t35;
t27 = t37 * t78 + t47 * t73;
t26 = -t37 * t73 + t47 * t78;
t19 = t26 * t72 + t27 * t77;
t18 = t26 * t77 - t27 * t72;
t7 = -pkin(5) * t26 + t9;
t4 = pkin(11) * t26 + t6;
t3 = pkin(5) * t35 - pkin(11) * t27 + t5;
t2 = t3 * t72 + t4 * t77;
t1 = t3 * t77 - t4 * t72;
t8 = (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (-t40 * mrSges(3,1) + t39 * mrSges(3,3) + Ifges(3,4) * t51 + Ifges(3,6) * t56 + Ifges(3,2) * t50 / 0.2e1) * t50 + (t38 * mrSges(3,1) - t39 * mrSges(3,2) + Ifges(3,3) * t56 / 0.2e1) * t56 + (-t32 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t49 + Ifges(4,2) * t42 / 0.2e1) * t42 + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + Ifges(4,3) * t49 / 0.2e1) * t49 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t25 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,6) * t47 + Ifges(5,2) * t36 / 0.2e1) * t36 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t35 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t67 * mrSges(2,1) + t58 * mrSges(2,3) + Ifges(2,4) * t61 + Ifges(2,2) * t60 / 0.2e1) * t60 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t19 / 0.2e1) * t19 + (t32 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,5) * t49 + Ifges(4,1) * t43 / 0.2e1) * t43 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t35 + Ifges(6,1) * t27 / 0.2e1) * t27 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t34 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t67 * mrSges(2,2) - t57 * mrSges(2,3) + Ifges(2,1) * t61 / 0.2e1) * t61 + (t25 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t47 + Ifges(5,1) * t37 / 0.2e1) * t37 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t35 / 0.2e1) * t35 + (t40 * mrSges(3,2) - t38 * mrSges(3,3) + Ifges(3,5) * t56 + Ifges(3,1) * t51 / 0.2e1) * t51 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) + t57 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t58 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t61 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t60 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + m(5) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(2) * (t57 ^ 2 + t58 ^ 2 + t67 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T  = t8;
