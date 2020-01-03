% Calculate joint inertia matrix for
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:05
% DurationCPUTime: 0.68s
% Computational Cost: add. (1085->189), mult. (2078->290), div. (0->0), fcn. (2199->8), ass. (0->71)
t66 = sin(pkin(9));
t67 = cos(pkin(9));
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t45 = t66 * t70 - t67 * t73;
t46 = t66 * t73 + t67 * t70;
t72 = cos(qJ(4));
t88 = t46 * t72;
t97 = Ifges(5,5) * t88 + Ifges(5,3) * t45;
t87 = -qJ(3) - pkin(6);
t52 = t87 * t73;
t81 = t87 * t70;
t33 = -t66 * t52 - t67 * t81;
t96 = t33 ^ 2;
t95 = 0.2e1 * t33;
t59 = -t73 * pkin(2) - pkin(1);
t94 = 0.2e1 * t59;
t57 = t66 * pkin(2) + pkin(7);
t92 = pkin(8) + t57;
t69 = sin(qJ(4));
t91 = Ifges(5,4) * t69;
t90 = Ifges(5,4) * t72;
t89 = t46 * t69;
t26 = t45 * pkin(3) - t46 * pkin(7) + t59;
t35 = -t67 * t52 + t66 * t81;
t10 = t69 * t26 + t72 * t35;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t48 = -t68 * t69 + t71 * t72;
t49 = t68 * t72 + t71 * t69;
t86 = Ifges(6,5) * t49 + Ifges(6,6) * t48;
t85 = Ifges(5,5) * t69 + Ifges(5,6) * t72;
t84 = t69 ^ 2 + t72 ^ 2;
t83 = t70 ^ 2 + t73 ^ 2;
t18 = t49 * t46;
t19 = t48 * t46;
t82 = Ifges(6,5) * t19 - Ifges(6,6) * t18 + Ifges(6,3) * t45;
t58 = -t67 * pkin(2) - pkin(3);
t30 = -t48 * mrSges(6,1) + t49 * mrSges(6,2);
t9 = t72 * t26 - t69 * t35;
t51 = -t72 * mrSges(5,1) + t69 * mrSges(5,2);
t80 = t69 * mrSges(5,1) + t72 * mrSges(5,2);
t40 = t92 * t69;
t41 = t92 * t72;
t24 = -t71 * t40 - t68 * t41;
t25 = -t68 * t40 + t71 * t41;
t79 = t24 * mrSges(6,1) - t25 * mrSges(6,2) + t86;
t4 = t45 * pkin(4) - pkin(8) * t88 + t9;
t7 = -pkin(8) * t89 + t10;
t2 = t71 * t4 - t68 * t7;
t3 = t68 * t4 + t71 * t7;
t78 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t82;
t77 = (t71 * mrSges(6,1) - t68 * mrSges(6,2)) * pkin(4);
t54 = Ifges(5,1) * t69 + t90;
t53 = Ifges(5,2) * t72 + t91;
t50 = -t72 * pkin(4) + t58;
t37 = t46 * mrSges(4,2);
t32 = Ifges(6,1) * t49 + Ifges(6,4) * t48;
t31 = Ifges(6,4) * t49 + Ifges(6,2) * t48;
t28 = t45 * mrSges(5,1) - mrSges(5,3) * t88;
t27 = -t45 * mrSges(5,2) - mrSges(5,3) * t89;
t23 = t80 * t46;
t15 = pkin(4) * t89 + t33;
t14 = Ifges(5,5) * t45 + (Ifges(5,1) * t72 - t91) * t46;
t13 = Ifges(5,6) * t45 + (-Ifges(5,2) * t69 + t90) * t46;
t12 = t45 * mrSges(6,1) - t19 * mrSges(6,3);
t11 = -t45 * mrSges(6,2) - t18 * mrSges(6,3);
t8 = t18 * mrSges(6,1) + t19 * mrSges(6,2);
t6 = Ifges(6,1) * t19 - Ifges(6,4) * t18 + Ifges(6,5) * t45;
t5 = Ifges(6,4) * t19 - Ifges(6,2) * t18 + Ifges(6,6) * t45;
t1 = [t70 * (Ifges(3,1) * t70 + Ifges(3,4) * t73) + t73 * (Ifges(3,4) * t70 + Ifges(3,2) * t73) - 0.2e1 * pkin(1) * (-t73 * mrSges(3,1) + t70 * mrSges(3,2)) + t37 * t94 + t23 * t95 - t18 * t5 + t19 * t6 + 0.2e1 * t10 * t27 + 0.2e1 * t9 * t28 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t15 * t8 + Ifges(2,3) + 0.2e1 * t83 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t94 - 0.2e1 * t35 * mrSges(4,3) + Ifges(4,2) * t45 + t82 + t97) * t45 + (mrSges(4,3) * t95 + Ifges(4,1) * t46 - t69 * t13 + t72 * t14 + (-Ifges(5,6) * t69 - (2 * Ifges(4,4))) * t45) * t46 + m(3) * (t83 * pkin(6) ^ 2 + pkin(1) ^ 2) + m(4) * (t35 ^ 2 + t59 ^ 2 + t96) + m(6) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t10 ^ 2 + t9 ^ 2 + t96); Ifges(3,5) * t70 + Ifges(3,6) * t73 + t49 * t6 / 0.2e1 + t50 * t8 + t58 * t23 - Ifges(4,6) * t45 + t48 * t5 / 0.2e1 - t18 * t31 / 0.2e1 + t19 * t32 / 0.2e1 - t35 * mrSges(4,2) + t24 * t12 + t25 * t11 + t15 * t30 + (t51 - mrSges(4,1)) * t33 + (-t70 * mrSges(3,1) - t73 * mrSges(3,2)) * pkin(6) + (-t2 * t49 + t3 * t48) * mrSges(6,3) + (t13 / 0.2e1 + t57 * t27 + t10 * mrSges(5,3)) * t72 + (t14 / 0.2e1 - t57 * t28 - t9 * mrSges(5,3)) * t69 + m(5) * (t58 * t33 + (t10 * t72 - t9 * t69) * t57) + m(6) * (t50 * t15 + t24 * t2 + t25 * t3) + (t72 * t54 / 0.2e1 - t69 * t53 / 0.2e1 + Ifges(4,5)) * t46 + (m(4) * (-t33 * t67 + t35 * t66) + (-t66 * t45 - t67 * t46) * mrSges(4,3)) * pkin(2) + (t85 + t86) * t45 / 0.2e1; 0.2e1 * t50 * t30 + t48 * t31 + t49 * t32 + 0.2e1 * t58 * t51 + t72 * t53 + t69 * t54 + Ifges(3,3) + Ifges(4,3) + m(6) * (t24 ^ 2 + t25 ^ 2 + t50 ^ 2) + m(5) * (t84 * t57 ^ 2 + t58 ^ 2) + m(4) * (t66 ^ 2 + t67 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t67 * mrSges(4,1) - t66 * mrSges(4,2)) * pkin(2) + 0.2e1 * (-t24 * t49 + t25 * t48) * mrSges(6,3) + 0.2e1 * t84 * t57 * mrSges(5,3); t45 * mrSges(4,1) + t49 * t11 + t48 * t12 + t69 * t27 + t72 * t28 + t37 + m(6) * (t48 * t2 + t49 * t3) + m(5) * (t69 * t10 + t72 * t9) + m(4) * t59; m(6) * (t48 * t24 + t49 * t25); m(4) + m(5) * t84 + m(6) * (t48 ^ 2 + t49 ^ 2); -Ifges(5,6) * t89 + t9 * mrSges(5,1) - t10 * mrSges(5,2) + (m(6) * (t2 * t71 + t3 * t68) + t68 * t11 + t71 * t12) * pkin(4) + t78 + t97; -t80 * t57 + (m(6) * (t24 * t71 + t25 * t68) + (t68 * t48 - t71 * t49) * mrSges(6,3)) * pkin(4) + t79 + t85; m(6) * (t48 * t71 + t49 * t68) * pkin(4) - t51 - t30; Ifges(5,3) + Ifges(6,3) + m(6) * (t68 ^ 2 + t71 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t77; t78; t79; -t30; Ifges(6,3) + t77; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
