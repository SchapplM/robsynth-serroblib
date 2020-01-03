% Calculate joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:32
% EndTime: 2020-01-03 12:03:33
% DurationCPUTime: 0.41s
% Computational Cost: add. (870->135), mult. (1590->182), div. (0->0), fcn. (1627->8), ass. (0->59)
t58 = cos(pkin(9));
t86 = t58 ^ 2;
t57 = sin(pkin(9));
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t43 = -t60 * t57 + t63 * t58;
t44 = t63 * t57 + t60 * t58;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t22 = t62 * t43 - t59 * t44;
t85 = t59 * pkin(4) * t22 * mrSges(6,3) + Ifges(5,5) * t44 + Ifges(5,6) * t43;
t23 = t59 * t43 + t62 * t44;
t9 = -t22 * mrSges(6,1) + t23 * mrSges(6,2);
t84 = 0.2e1 * t9;
t24 = -t43 * mrSges(5,1) + t44 * mrSges(5,2);
t83 = 0.2e1 * t24;
t82 = t44 * pkin(8);
t64 = cos(qJ(2));
t81 = t64 * pkin(1);
t80 = Ifges(6,5) * t23 + Ifges(6,6) * t22;
t61 = sin(qJ(2));
t50 = t61 * pkin(1) + qJ(3);
t32 = (-pkin(7) - t50) * t57;
t54 = t58 * pkin(7);
t33 = t58 * t50 + t54;
t18 = t60 * t32 + t63 * t33;
t46 = (-pkin(7) - qJ(3)) * t57;
t48 = t58 * qJ(3) + t54;
t26 = t60 * t46 + t63 * t48;
t79 = t57 ^ 2 + t86;
t78 = 2 * mrSges(5,3);
t77 = 0.2e1 * mrSges(6,3);
t76 = t62 * t23 * mrSges(6,3);
t51 = -t58 * pkin(3) - pkin(2);
t47 = -t58 * mrSges(4,1) + t57 * mrSges(4,2);
t17 = t63 * t32 - t60 * t33;
t25 = t63 * t46 - t60 * t48;
t75 = t79 * qJ(3);
t10 = t17 - t82;
t37 = t43 * pkin(8);
t11 = t37 + t18;
t2 = t62 * t10 - t59 * t11;
t3 = t59 * t10 + t62 * t11;
t74 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t80;
t14 = t25 - t82;
t15 = t37 + t26;
t7 = t62 * t14 - t59 * t15;
t8 = t59 * t14 + t62 * t15;
t73 = t7 * mrSges(6,1) - t8 * mrSges(6,2) + t80;
t72 = 0.2e1 * t79 * mrSges(4,3);
t28 = -t43 * pkin(4) + t51;
t71 = Ifges(5,1) * t44 ^ 2 + Ifges(6,1) * t23 ^ 2 + Ifges(4,2) * t86 + Ifges(3,3) + (Ifges(4,1) * t57 + 0.2e1 * Ifges(4,4) * t58) * t57 + (0.2e1 * Ifges(5,4) * t44 + Ifges(5,2) * t43) * t43 + (0.2e1 * Ifges(6,4) * t23 + Ifges(6,2) * t22) * t22;
t70 = (t64 * mrSges(3,1) - t61 * mrSges(3,2)) * pkin(1);
t69 = (t62 * mrSges(6,1) - t59 * mrSges(6,2)) * pkin(4);
t68 = t24 + t47 + t9;
t52 = -pkin(2) - t81;
t45 = t51 - t81;
t27 = t28 - t81;
t1 = [t71 + m(4) * (t79 * t50 ^ 2 + t52 ^ 2) + (-t17 * t44 + t18 * t43) * t78 + (-t2 * t23 + t3 * t22) * t77 + m(3) * (t61 ^ 2 + t64 ^ 2) * pkin(1) ^ 2 + t50 * t72 + m(6) * (t2 ^ 2 + t27 ^ 2 + t3 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t45 ^ 2) + t45 * t83 + 0.2e1 * t52 * t47 + t27 * t84 + Ifges(2,3) + 0.2e1 * t70; t70 + t71 + m(4) * (-pkin(2) * t52 + t50 * t75) + (t79 * t50 + t75) * mrSges(4,3) + ((-t17 - t25) * t44 + (t18 + t26) * t43) * mrSges(5,3) + ((-t2 - t7) * t23 + (t3 + t8) * t22) * mrSges(6,3) + m(6) * (t7 * t2 + t28 * t27 + t8 * t3) + m(5) * (t25 * t17 + t26 * t18 + t51 * t45) + (t27 + t28) * t9 + (-pkin(2) + t52) * t47 + (t45 + t51) * t24; -0.2e1 * pkin(2) * t47 + t51 * t83 + t28 * t84 + (t8 * t22 - t7 * t23) * t77 + (-t25 * t44 + t26 * t43) * t78 + qJ(3) * t72 + m(6) * (t28 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2 + t51 ^ 2) + m(4) * (t79 * qJ(3) ^ 2 + pkin(2) ^ 2) + t71; m(4) * t52 + m(5) * t45 + m(6) * t27 + t68; -m(4) * pkin(2) + m(5) * t51 + m(6) * t28 + t68; m(4) + m(5) + m(6); t17 * mrSges(5,1) - t18 * mrSges(5,2) + (-t76 + m(6) * (t2 * t62 + t3 * t59)) * pkin(4) + t74 + t85; t25 * mrSges(5,1) - t26 * mrSges(5,2) + (-t76 + m(6) * (t59 * t8 + t62 * t7)) * pkin(4) + t73 + t85; 0; Ifges(5,3) + Ifges(6,3) + m(6) * (t59 ^ 2 + t62 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t69; t74; t73; 0; Ifges(6,3) + t69; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
