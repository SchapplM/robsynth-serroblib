% Calculate joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:04
% EndTime: 2020-01-03 12:13:05
% DurationCPUTime: 0.40s
% Computational Cost: add. (771->146), mult. (1396->206), div. (0->0), fcn. (1182->8), ass. (0->82)
t68 = cos(qJ(4));
t62 = t68 ^ 2;
t63 = sin(qJ(5));
t64 = sin(qJ(4));
t67 = cos(qJ(5));
t39 = -t63 * t64 + t67 * t68;
t92 = mrSges(6,3) * t39;
t102 = t63 * pkin(4) * t92 + Ifges(5,5) * t64 + Ifges(5,6) * t68;
t70 = cos(qJ(2));
t53 = t70 * pkin(1) + pkin(2);
t65 = sin(qJ(3));
t69 = cos(qJ(3));
t66 = sin(qJ(2));
t99 = pkin(1) * t66;
t31 = t65 * t53 + t69 * t99;
t29 = pkin(8) + t31;
t18 = (-pkin(9) - t29) * t64;
t60 = t68 * pkin(9);
t19 = t68 * t29 + t60;
t4 = t63 * t18 + t67 * t19;
t1 = t4 * t92;
t30 = t69 * t53 - t65 * t99;
t28 = -pkin(3) - t30;
t44 = -t68 * mrSges(5,1) + t64 * mrSges(5,2);
t16 = t28 * t44;
t61 = t64 ^ 2;
t94 = mrSges(5,3) * t61;
t22 = t29 * t94;
t93 = mrSges(5,3) * t62;
t23 = t29 * t93;
t25 = t30 * mrSges(4,1);
t40 = t63 * t68 + t67 * t64;
t15 = -t39 * mrSges(6,1) + t40 * mrSges(6,2);
t96 = t68 * pkin(4);
t24 = t28 - t96;
t5 = t24 * t15;
t101 = t1 + t16 + t22 + t23 + t25 + t5;
t95 = t69 * pkin(2);
t52 = -pkin(3) - t95;
t32 = t52 * t44;
t97 = t65 * pkin(2);
t51 = pkin(8) + t97;
t41 = t51 * t94;
t42 = t51 * t93;
t57 = mrSges(4,1) * t95;
t35 = (-pkin(9) - t51) * t64;
t36 = t68 * t51 + t60;
t14 = t63 * t35 + t67 * t36;
t6 = t14 * t92;
t54 = -pkin(3) - t96;
t43 = t54 - t95;
t9 = t43 * t15;
t100 = t32 + t41 + t42 + t57 + t6 + t9;
t98 = pkin(3) * t44;
t91 = t31 * mrSges(4,2);
t90 = t40 * mrSges(6,3);
t89 = Ifges(6,5) * t40 + Ifges(6,6) * t39;
t88 = t61 + t62;
t87 = -0.2e1 * t90;
t86 = mrSges(4,2) * t97;
t85 = t67 * t90;
t84 = t88 * t51;
t83 = Ifges(6,1) * t40 ^ 2 + Ifges(5,2) * t62 + Ifges(4,3) + (Ifges(5,1) * t64 + 0.2e1 * Ifges(5,4) * t68) * t64 + (0.2e1 * Ifges(6,4) * t40 + Ifges(6,2) * t39) * t39;
t82 = -t64 * mrSges(5,1) - t68 * mrSges(5,2);
t3 = t67 * t18 - t63 * t19;
t81 = t3 * mrSges(6,1) - t4 * mrSges(6,2) + t89;
t80 = Ifges(3,3) + t83;
t13 = t67 * t35 - t63 * t36;
t79 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t89;
t45 = (-pkin(9) - pkin(8)) * t64;
t46 = t68 * pkin(8) + t60;
t20 = t67 * t45 - t63 * t46;
t21 = t63 * t45 + t67 * t46;
t78 = t20 * mrSges(6,1) - t21 * mrSges(6,2) + t89;
t77 = (t70 * mrSges(3,1) - t66 * mrSges(3,2)) * pkin(1);
t76 = (t67 * mrSges(6,1) - t63 * mrSges(6,2)) * pkin(4);
t10 = t21 * t92;
t12 = t54 * t15;
t55 = pkin(8) * t94;
t56 = pkin(8) * t93;
t75 = t10 + t12 + t55 + t56 + t83 - t98;
t2 = [t80 + m(5) * (t88 * t29 ^ 2 + t28 ^ 2) + m(3) * (t66 ^ 2 + t70 ^ 2) * pkin(1) ^ 2 + m(6) * (t24 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + 0.2e1 * t25 + 0.2e1 * t22 + 0.2e1 * t23 + 0.2e1 * t16 + 0.2e1 * t5 - 0.2e1 * t91 + 0.2e1 * t1 + t3 * t87 + Ifges(2,3) + 0.2e1 * t77; t100 + (-t31 - t97) * mrSges(4,2) + m(5) * (t52 * t28 + t29 * t84) + t80 + t77 + m(6) * (t13 * t3 + t14 * t4 + t43 * t24) + m(4) * (t30 * t69 + t31 * t65) * pkin(2) + (-t13 - t3) * t90 + t101; -0.2e1 * t86 + t13 * t87 + 0.2e1 * t32 + 0.2e1 * t41 + 0.2e1 * t42 + 0.2e1 * t57 + 0.2e1 * t6 + 0.2e1 * t9 + m(6) * (t13 ^ 2 + t14 ^ 2 + t43 ^ 2) + m(5) * (t88 * t51 ^ 2 + t52 ^ 2) + m(4) * (t65 ^ 2 + t69 ^ 2) * pkin(2) ^ 2 + t80; t75 + m(5) * (t88 * t29 * pkin(8) - pkin(3) * t28) + m(6) * (t20 * t3 + t21 * t4 + t54 * t24) - t91 + (-t20 - t3) * t90 + t101; t75 + m(5) * (-pkin(3) * t52 + pkin(8) * t84) + m(6) * (t20 * t13 + t21 * t14 + t54 * t43) - t86 + (-t13 - t20) * t90 + t100; t20 * t87 - 0.2e1 * t98 + 0.2e1 * t10 + 0.2e1 * t12 + 0.2e1 * t55 + 0.2e1 * t56 + m(6) * (t20 ^ 2 + t21 ^ 2 + t54 ^ 2) + m(5) * (t88 * pkin(8) ^ 2 + pkin(3) ^ 2) + t83; t82 * t29 + (-t85 + m(6) * (t3 * t67 + t4 * t63)) * pkin(4) + t81 + t102; t82 * t51 + (-t85 + m(6) * (t13 * t67 + t14 * t63)) * pkin(4) + t79 + t102; t82 * pkin(8) + (-t85 + m(6) * (t20 * t67 + t21 * t63)) * pkin(4) + t78 + t102; Ifges(5,3) + Ifges(6,3) + m(6) * (t63 ^ 2 + t67 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t76; t81; t79; t78; Ifges(6,3) + t76; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
