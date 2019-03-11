% Calculate kinetic energy for
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:25
% EndTime: 2019-03-08 19:29:27
% DurationCPUTime: 1.38s
% Computational Cost: add. (5501->161), mult. (9957->235), div. (0->0), fcn. (8538->14), ass. (0->64)
t63 = V_base(5) * qJ(1) + V_base(1);
t64 = -V_base(4) * qJ(1) + V_base(2);
t68 = sin(pkin(10));
t72 = cos(pkin(10));
t56 = -t63 * t68 + t72 * t64;
t73 = cos(pkin(6));
t59 = t68 * V_base(5) + t72 * V_base(4);
t83 = pkin(7) * t59;
t49 = V_base(6) * pkin(1) - t73 * t83 + t56;
t58 = -t68 * V_base(4) + t72 * V_base(5);
t65 = V_base(3) + qJD(1);
t69 = sin(pkin(6));
t53 = -pkin(1) * t58 - t69 * t83 + t65;
t84 = t49 * t73 + t53 * t69;
t57 = t72 * t63 + t68 * t64;
t79 = t58 * t73 + t69 * V_base(6);
t46 = t79 * pkin(7) + t57;
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t33 = -t46 * t76 + t84 * t78;
t48 = t59 * t78 + t79 * t76;
t55 = -t58 * t69 + t73 * V_base(6) + qJD(2);
t29 = pkin(2) * t55 - qJ(3) * t48 + t33;
t34 = t78 * t46 + t84 * t76;
t47 = -t76 * t59 + t79 * t78;
t32 = qJ(3) * t47 + t34;
t67 = sin(pkin(11));
t71 = cos(pkin(11));
t23 = t67 * t29 + t71 * t32;
t21 = pkin(8) * t55 + t23;
t42 = -t49 * t69 + t73 * t53;
t35 = -pkin(2) * t47 + qJD(3) + t42;
t40 = t47 * t71 - t48 * t67;
t41 = t47 * t67 + t48 * t71;
t25 = -pkin(3) * t40 - pkin(8) * t41 + t35;
t75 = sin(qJ(4));
t82 = cos(qJ(4));
t12 = t82 * t21 + t75 * t25;
t39 = qJD(4) - t40;
t10 = qJ(5) * t39 + t12;
t22 = t29 * t71 - t67 * t32;
t20 = -pkin(3) * t55 - t22;
t37 = t41 * t75 - t82 * t55;
t38 = t82 * t41 + t75 * t55;
t15 = pkin(4) * t37 - qJ(5) * t38 + t20;
t66 = sin(pkin(12));
t70 = cos(pkin(12));
t6 = t70 * t10 + t66 * t15;
t5 = -t10 * t66 + t70 * t15;
t11 = -t75 * t21 + t82 * t25;
t9 = -t39 * pkin(4) + qJD(5) - t11;
t77 = cos(qJ(6));
t74 = sin(qJ(6));
t36 = qJD(6) + t37;
t27 = t38 * t70 + t39 * t66;
t26 = -t38 * t66 + t39 * t70;
t17 = t26 * t74 + t27 * t77;
t16 = t26 * t77 - t27 * t74;
t7 = -t26 * pkin(5) + t9;
t4 = pkin(9) * t26 + t6;
t3 = pkin(5) * t37 - pkin(9) * t27 + t5;
t2 = t3 * t74 + t4 * t77;
t1 = t3 * t77 - t4 * t74;
t8 = (-t42 * mrSges(3,1) + t34 * mrSges(3,3) + Ifges(3,4) * t48 + Ifges(3,2) * t47 / 0.2e1) * t47 + (t33 * mrSges(3,1) + t22 * mrSges(4,1) - t34 * mrSges(3,2) - t23 * mrSges(4,2) + Ifges(3,5) * t48 + Ifges(4,5) * t41 + Ifges(3,6) * t47 + Ifges(4,6) * t40 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t55) * t55 + (-t65 * mrSges(2,1) + t57 * mrSges(2,3) + Ifges(2,4) * t59 + Ifges(2,2) * t58 / 0.2e1) * t58 + (t20 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t38 + Ifges(6,5) * t27 - Ifges(5,6) * t39 + Ifges(6,6) * t26 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t37) * t37 + (-t35 * mrSges(4,1) + t23 * mrSges(4,3) + Ifges(4,4) * t41 + Ifges(4,2) * t40 / 0.2e1) * t40 + (t20 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t39 + Ifges(5,1) * t38 / 0.2e1) * t38 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t65 * mrSges(2,2) - t56 * mrSges(2,3) + Ifges(2,1) * t59 / 0.2e1) * t59 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t36 + Ifges(7,1) * t17 / 0.2e1) * t17 + (t42 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,1) * t48 / 0.2e1) * t48 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t36 + Ifges(7,2) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (V_base(2) * mrSges(1,1) + t56 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t57 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t59 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t58 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t36 / 0.2e1) * t36 + (t35 * mrSges(4,2) - t22 * mrSges(4,3) + Ifges(4,1) * t41 / 0.2e1) * t41 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t39 / 0.2e1) * t39 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t20 ^ 2) / 0.2e1 + m(4) * (t22 ^ 2 + t23 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t42 ^ 2) / 0.2e1 + m(2) * (t56 ^ 2 + t57 ^ 2 + t65 ^ 2) / 0.2e1;
T  = t8;
