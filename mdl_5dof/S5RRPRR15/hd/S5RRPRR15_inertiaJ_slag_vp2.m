% Calculate joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR15_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:15
% DurationCPUTime: 0.70s
% Computational Cost: add. (689->184), mult. (1247->259), div. (0->0), fcn. (1089->6), ass. (0->76)
t100 = pkin(3) + pkin(6);
t65 = sin(qJ(2));
t68 = cos(qJ(2));
t99 = t65 ^ 2 + t68 ^ 2;
t98 = -m(4) * pkin(2) + mrSges(4,2);
t64 = sin(qJ(4));
t97 = -t64 / 0.2e1;
t69 = -pkin(2) - pkin(7);
t96 = -pkin(8) + t69;
t95 = Ifges(5,4) * t64;
t67 = cos(qJ(4));
t94 = Ifges(5,4) * t67;
t93 = t67 * mrSges(5,1);
t92 = t67 * t68;
t84 = -t65 * qJ(3) - pkin(1);
t28 = t69 * t68 + t84;
t46 = t100 * t65;
t10 = t67 * t28 + t64 * t46;
t91 = t99 * pkin(6) ^ 2;
t47 = t100 * t68;
t90 = t64 ^ 2 + t67 ^ 2;
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t33 = -t63 * t67 - t66 * t64;
t76 = t63 * t64 - t66 * t67;
t89 = t33 ^ 2 + t76 ^ 2;
t23 = t76 * t68;
t24 = t33 * t68;
t88 = Ifges(6,5) * t24 + Ifges(6,6) * t23 + Ifges(6,3) * t65;
t87 = m(5) * t90;
t86 = t90 * mrSges(5,3);
t85 = -mrSges(6,1) * t76 + t33 * mrSges(6,2);
t37 = t67 * t46;
t6 = t65 * pkin(4) + t37 + (pkin(8) * t68 - t28) * t64;
t8 = -pkin(8) * t92 + t10;
t2 = t66 * t6 - t63 * t8;
t3 = t63 * t6 + t66 * t8;
t82 = -t2 * t76 - t33 * t3;
t9 = -t64 * t28 + t37;
t81 = t64 * t10 + t67 * t9;
t80 = -t64 * mrSges(5,2) + t93;
t79 = -Ifges(5,5) * t64 - Ifges(5,6) * t67;
t40 = t96 * t64;
t41 = t96 * t67;
t15 = -t63 * t40 + t66 * t41;
t16 = t66 * t40 + t63 * t41;
t78 = -t15 * t76 - t33 * t16;
t77 = t33 * t63 + t66 * t76;
t30 = Ifges(6,6) * t33;
t31 = Ifges(6,5) * t76;
t75 = t15 * mrSges(6,1) - t16 * mrSges(6,2) + t30 - t31;
t74 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t88;
t73 = (t66 * mrSges(6,1) - t63 * mrSges(6,2)) * pkin(4);
t70 = qJ(3) ^ 2;
t53 = Ifges(5,5) * t67;
t52 = Ifges(5,3) * t65;
t48 = t64 * pkin(4) + qJ(3);
t45 = Ifges(5,1) * t67 - t95;
t44 = -Ifges(5,2) * t64 + t94;
t43 = t64 * mrSges(5,1) + t67 * mrSges(5,2);
t42 = -t68 * pkin(2) + t84;
t39 = -t65 * mrSges(5,2) - mrSges(5,3) * t92;
t38 = t64 * t68 * mrSges(5,3) + t65 * mrSges(5,1);
t27 = t80 * t68;
t26 = pkin(4) * t92 + t47;
t22 = Ifges(5,5) * t65 + (-Ifges(5,1) * t64 - t94) * t68;
t21 = Ifges(5,6) * t65 + (-Ifges(5,2) * t67 - t95) * t68;
t18 = t65 * mrSges(6,1) - t24 * mrSges(6,3);
t17 = -t65 * mrSges(6,2) + t23 * mrSges(6,3);
t14 = -Ifges(6,1) * t76 + Ifges(6,4) * t33;
t13 = -Ifges(6,4) * t76 + Ifges(6,2) * t33;
t12 = -t33 * mrSges(6,1) - mrSges(6,2) * t76;
t7 = -t23 * mrSges(6,1) + t24 * mrSges(6,2);
t5 = Ifges(6,1) * t24 + Ifges(6,4) * t23 + Ifges(6,5) * t65;
t4 = Ifges(6,4) * t24 + Ifges(6,2) * t23 + Ifges(6,6) * t65;
t1 = [0.2e1 * t10 * t39 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + t23 * t4 + t24 * t5 + 0.2e1 * t26 * t7 + 0.2e1 * t47 * t27 + 0.2e1 * t9 * t38 + Ifges(2,3) + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t42 * mrSges(4,3) + t52 + (Ifges(4,2) + Ifges(3,1)) * t65 + t88) * t65 + (0.2e1 * t42 * mrSges(4,2) + 0.2e1 * pkin(1) * mrSges(3,1) - t64 * t22 - t67 * t21 + (Ifges(3,2) + Ifges(4,3)) * t68 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t79) * t65) * t68 + m(6) * (t2 ^ 2 + t26 ^ 2 + t3 ^ 2) + m(5) * (t10 ^ 2 + t47 ^ 2 + t9 ^ 2) + m(3) * (pkin(1) ^ 2 + t91) + m(4) * (t42 ^ 2 + t91) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(6) * t99; -t76 * t5 / 0.2e1 + t47 * t43 + t48 * t7 + t26 * t12 + qJ(3) * t27 + t33 * t4 / 0.2e1 + t16 * t17 + t15 * t18 + t23 * t13 / 0.2e1 + t24 * t14 / 0.2e1 - t82 * mrSges(6,3) + (t22 / 0.2e1 - t9 * mrSges(5,3) + t69 * t38) * t67 + (-t21 / 0.2e1 - t10 * mrSges(5,3) + t69 * t39) * t64 + m(6) * (t15 * t2 + t16 * t3 + t48 * t26) + m(5) * (qJ(3) * t47 + t69 * t81) + (Ifges(5,6) * t97 + t53 / 0.2e1 - t31 / 0.2e1 + t30 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1)) * t65 + (-Ifges(4,5) + Ifges(3,6) + t45 * t97 - t67 * t44 / 0.2e1 + qJ(3) * mrSges(4,1)) * t68 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t68 + (-mrSges(3,1) + t98) * t65) * pkin(6); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t48 * t12 + t33 * t13 - t76 * t14 - t64 * t44 + t67 * t45 + Ifges(4,1) + Ifges(3,3) + m(6) * (t15 ^ 2 + t16 ^ 2 + t48 ^ 2) + m(5) * (t90 * t69 ^ 2 + t70) + m(4) * (pkin(2) ^ 2 + t70) + 0.2e1 * (t43 + mrSges(4,3)) * qJ(3) - 0.2e1 * t78 * mrSges(6,3) - 0.2e1 * t69 * t86; -t33 * t17 - t76 * t18 + t67 * t38 + t64 * t39 + (m(4) * pkin(6) + mrSges(4,1)) * t65 + m(6) * t82 + m(5) * t81; m(6) * t78 - t89 * mrSges(6,3) + t69 * t87 - t86 + t98; m(6) * t89 + m(4) + t87; t9 * mrSges(5,1) - t10 * mrSges(5,2) + t52 + t79 * t68 + (m(6) * (t2 * t66 + t3 * t63) + t66 * t18 + t63 * t17) * pkin(4) + t74; t69 * t93 + t53 + (-t69 * mrSges(5,2) - Ifges(5,6)) * t64 + (m(6) * (t15 * t66 + t16 * t63) + t77 * mrSges(6,3)) * pkin(4) + t75; -m(6) * pkin(4) * t77 + t80 + t85; Ifges(5,3) + Ifges(6,3) + m(6) * (t63 ^ 2 + t66 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t73; t74; t75; t85; Ifges(6,3) + t73; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
