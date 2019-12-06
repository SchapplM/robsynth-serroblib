% Calculate joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaJ_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:57
% EndTime: 2019-12-05 18:52:58
% DurationCPUTime: 0.38s
% Computational Cost: add. (468->109), mult. (1197->159), div. (0->0), fcn. (1143->8), ass. (0->64)
t61 = cos(qJ(3));
t54 = t61 ^ 2;
t55 = sin(qJ(5));
t59 = cos(qJ(5));
t36 = -t59 * mrSges(6,1) + t55 * mrSges(6,2);
t101 = mrSges(5,1) - t36;
t56 = sin(qJ(4));
t60 = cos(qJ(4));
t79 = t55 ^ 2 + t59 ^ 2;
t100 = (t101 * t60 + (t79 * mrSges(6,3) - mrSges(5,2)) * t56) * pkin(2);
t57 = sin(qJ(3));
t58 = sin(qJ(2));
t62 = cos(qJ(2));
t78 = t57 ^ 2 + t54;
t99 = ((t61 * mrSges(4,1) - t57 * mrSges(4,2) + mrSges(3,1)) * t62 + (t78 * mrSges(4,3) - mrSges(3,2)) * t58) * pkin(1);
t34 = t56 * t57 - t60 * t61;
t35 = t56 * t61 + t60 * t57;
t82 = t55 * mrSges(6,3);
t14 = -t34 * mrSges(6,2) - t35 * t82;
t80 = t56 * t59;
t95 = pkin(2) * t56;
t98 = -t34 * mrSges(5,3) * t95 + pkin(2) * t14 * t80 + Ifges(4,5) * t57 + Ifges(4,6) * t61;
t96 = pkin(1) * t58;
t22 = t35 * t96;
t97 = t22 ^ 2;
t94 = t61 * pkin(2);
t93 = Ifges(6,4) * t55;
t92 = Ifges(6,4) * t59;
t24 = t34 * t96;
t41 = -t62 * pkin(1) - t94;
t11 = t55 * t24 + t59 * t41;
t86 = t35 * t59;
t15 = t34 * mrSges(6,1) - mrSges(6,3) * t86;
t91 = t11 * t15;
t12 = -t59 * t24 + t55 * t41;
t90 = t12 * t14;
t72 = mrSges(6,1) * t55 + mrSges(6,2) * t59;
t13 = t72 * t35;
t89 = t22 * t13;
t88 = t22 * t60;
t87 = t35 * t55;
t17 = t34 * mrSges(5,1) + t35 * mrSges(5,2);
t85 = t41 * t17;
t63 = pkin(2) ^ 2;
t84 = t56 ^ 2 * t63;
t64 = pkin(1) ^ 2;
t83 = t58 ^ 2 * t64;
t81 = t55 * t56;
t38 = Ifges(6,5) * t55 + Ifges(6,6) * t59;
t39 = Ifges(6,2) * t59 + t93;
t40 = Ifges(6,1) * t55 + t92;
t77 = t59 * t39 + t55 * t40 + Ifges(5,3);
t71 = Ifges(6,5) * t86 - Ifges(6,6) * t87 + Ifges(6,3) * t34;
t70 = -t14 * t55 - t15 * t59 - t17;
t5 = Ifges(6,6) * t34 + (-Ifges(6,2) * t55 + t92) * t35;
t6 = Ifges(6,5) * t34 + (Ifges(6,1) * t59 - t93) * t35;
t69 = -t39 * t87 / 0.2e1 + t40 * t86 / 0.2e1 + t55 * t6 / 0.2e1 + t59 * t5 / 0.2e1 + Ifges(5,5) * t35 + (t38 / 0.2e1 - Ifges(5,6)) * t34;
t68 = (t22 * t35 + t24 * t34) * mrSges(5,3);
t67 = -t15 * t81 + (-t35 * mrSges(5,3) - t13) * t60;
t66 = Ifges(5,1) * t35 ^ 2 + Ifges(4,2) * t54 - t5 * t87 + t6 * t86 + Ifges(3,3) + (Ifges(4,1) * t57 + 0.2e1 * Ifges(4,4) * t61) * t57 + (-0.2e1 * Ifges(5,4) * t35 + Ifges(5,2) * t34 + t71) * t34;
t65 = t12 * t59 * mrSges(6,3) + t24 * mrSges(5,2) - t101 * t22 - t11 * t82 + t69;
t48 = t62 ^ 2 * t64;
t47 = t60 ^ 2 * t63;
t1 = [t66 + 0.2e1 * t68 + 0.2e1 * t99 + m(4) * (t78 * t83 + t48) + m(6) * (t11 ^ 2 + t12 ^ 2 + t97) + m(5) * (t24 ^ 2 + t41 ^ 2 + t97) + m(3) * (t48 + t83) + 0.2e1 * t85 + 0.2e1 * t89 + 0.2e1 * t90 + 0.2e1 * t91 + Ifges(2,3); t91 + t90 + t89 + t85 + t68 + t99 + (-m(5) * t41 + m(6) * (-t11 * t59 - t12 * t55) + t70) * t94 + t66; (m(6) * t79 + m(5)) * t63 * t54 + 0.2e1 * t70 * t94 + t66; t65 + (m(6) * (-t11 * t81 + t12 * t80 - t88) + m(5) * (-t24 * t56 - t88) + t67) * pkin(2) + (-mrSges(4,1) * t57 - mrSges(4,2) * t61) * t96 + t98; t67 * pkin(2) + t69 + t98; Ifges(4,3) + m(6) * (t79 * t84 + t47) + m(5) * (t47 + t84) + 0.2e1 * t100 + t77; t65; t69; t100 + t77; t77; t11 * mrSges(6,1) - t12 * mrSges(6,2) + t71; t36 * t94 + t71; -t72 * t95 + t38; t38; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
