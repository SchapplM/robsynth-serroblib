% Calculate joint inertia matrix for
% S5RRPRR6
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:42
% EndTime: 2019-12-05 18:35:44
% DurationCPUTime: 0.58s
% Computational Cost: add. (828->162), mult. (1596->219), div. (0->0), fcn. (1436->8), ass. (0->81)
t64 = sin(pkin(9));
t70 = cos(qJ(4));
t92 = t64 * t70;
t108 = Ifges(5,5) * t92;
t66 = sin(qJ(5));
t67 = sin(qJ(4));
t69 = cos(qJ(5));
t48 = t66 * t70 + t67 * t69;
t30 = t48 * t64;
t47 = -t66 * t67 + t69 * t70;
t31 = t47 * t64;
t89 = Ifges(6,5) * t31 - Ifges(6,6) * t30;
t93 = t64 * t67;
t109 = -Ifges(5,6) * t93 + t108 + t89;
t65 = cos(pkin(9));
t51 = -pkin(3) * t65 - pkin(7) * t64 - pkin(2);
t71 = cos(qJ(2));
t97 = pkin(1) * t71;
t38 = t51 - t97;
t33 = t70 * t38;
t84 = pkin(8) * t92;
t68 = sin(qJ(2));
t58 = pkin(1) * t68 + qJ(3);
t94 = t58 * t67;
t10 = -t84 + t33 + (-pkin(4) - t94) * t65;
t91 = t65 * t70;
t17 = t67 * t38 + t58 * t91;
t85 = pkin(8) * t93;
t14 = -t85 + t17;
t2 = t10 * t69 - t14 * t66;
t3 = t10 * t66 + t14 * t69;
t107 = t2 * mrSges(6,1) - t3 * mrSges(6,2);
t41 = t70 * t51;
t86 = qJ(3) * t67;
t15 = -t84 + t41 + (-pkin(4) - t86) * t65;
t25 = qJ(3) * t91 + t67 * t51;
t20 = -t85 + t25;
t5 = t15 * t69 - t20 * t66;
t6 = t15 * t66 + t20 * t69;
t106 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t62 = t64 ^ 2;
t63 = t65 ^ 2;
t105 = (t62 + t63) * mrSges(4,3);
t104 = qJ(3) + t58;
t11 = mrSges(6,1) * t30 + mrSges(6,2) * t31;
t103 = 0.2e1 * t11;
t22 = mrSges(6,2) * t65 - mrSges(6,3) * t30;
t102 = 0.2e1 * t22;
t23 = -mrSges(6,1) * t65 - mrSges(6,3) * t31;
t101 = 0.2e1 * t23;
t44 = mrSges(5,2) * t65 - mrSges(5,3) * t93;
t100 = 0.2e1 * t44;
t45 = -mrSges(5,1) * t65 - mrSges(5,3) * t92;
t99 = 0.2e1 * t45;
t98 = m(6) * pkin(4);
t90 = -Ifges(6,3) - Ifges(5,3);
t87 = qJ(3) * t58;
t52 = -t65 * mrSges(4,1) + t64 * mrSges(4,2);
t83 = t47 * mrSges(6,1) - t48 * mrSges(6,2);
t9 = -Ifges(6,3) * t65 + t89;
t81 = (mrSges(3,1) * t71 - mrSges(3,2) * t68) * pkin(1);
t80 = (mrSges(6,1) * t69 - mrSges(6,2) * t66) * pkin(4);
t79 = t48 * t22 + t47 * t23 + t67 * t44 + t70 * t45 + t52;
t37 = (mrSges(5,1) * t67 + mrSges(5,2) * t70) * t64;
t78 = 0.2e1 * t64 * t37 + 0.2e1 * t105;
t77 = Ifges(6,1) * t31 ^ 2 + Ifges(3,3) + (-0.2e1 * Ifges(6,4) * t31 + Ifges(6,2) * t30) * t30 + ((Ifges(5,1) * t70 - Ifges(5,4) * t67) * t92 + Ifges(4,1) * t64) * t64 + (0.2e1 * Ifges(4,4) * t64 - t108 - t9 + (Ifges(4,2) + Ifges(5,3)) * t65 - t109) * t65;
t76 = t90 * t65 + (t22 * t66 + t23 * t69) * pkin(4) + t109;
t29 = -Ifges(5,6) * t65 + (Ifges(5,4) * t70 - Ifges(5,2) * t67) * t64;
t75 = -t29 * t93 + t77;
t72 = qJ(3) ^ 2;
t60 = t62 * t72;
t59 = -pkin(2) - t97;
t57 = t58 ^ 2;
t56 = pkin(4) * t93;
t53 = t62 * t57;
t50 = t62 * t87;
t46 = qJ(3) * t64 + t56;
t34 = t58 * t64 + t56;
t24 = -t65 * t86 + t41;
t16 = -t65 * t94 + t33;
t1 = [t75 + m(6) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t53) + m(4) * (t57 * t63 + t59 ^ 2 + t53) + t17 * t100 + t16 * t99 + 0.2e1 * t59 * t52 + t34 * t103 + t3 * t102 + t2 * t101 + m(3) * (t68 ^ 2 + t71 ^ 2) * pkin(1) ^ 2 + 0.2e1 * t81 + Ifges(2,3) + t78 * t58; t77 + m(6) * (t2 * t5 + t3 * t6 + t34 * t46) + m(5) * (t16 * t24 + t17 * t25 + t50) + m(4) * (-pkin(2) * t59 + t63 * t87 + t50) + (-pkin(2) + t59) * t52 + (t24 + t16) * t45 + (t25 + t17) * t44 + (t5 + t2) * t23 + (t6 + t3) * t22 + (t46 + t34) * t11 + t81 + (t104 * t37 - t29 * t67) * t64 + t104 * t105; t75 + m(6) * (t46 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t60) + m(4) * (pkin(2) ^ 2 + t63 * t72 + t60) + t25 * t100 + t24 * t99 + t46 * t103 - 0.2e1 * pkin(2) * t52 + t6 * t102 + t5 * t101 + t78 * qJ(3); m(6) * (t2 * t47 + t3 * t48) + m(5) * (t16 * t70 + t17 * t67) + m(4) * t59 + t79; -m(4) * pkin(2) + m(6) * (t47 * t5 + t48 * t6) + m(5) * (t24 * t70 + t25 * t67) + t79; m(4) + m(5) * (t67 ^ 2 + t70 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2); t16 * mrSges(5,1) - t17 * mrSges(5,2) + (t2 * t69 + t3 * t66) * t98 + t76 + t107; t24 * mrSges(5,1) - t25 * mrSges(5,2) + (t5 * t69 + t6 * t66) * t98 + t76 + t106; t70 * mrSges(5,1) - t67 * mrSges(5,2) + (t47 * t69 + t48 * t66) * t98 + t83; m(6) * (t66 ^ 2 + t69 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t80 - t90; t9 + t107; t9 + t106; t83; Ifges(6,3) + t80; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
