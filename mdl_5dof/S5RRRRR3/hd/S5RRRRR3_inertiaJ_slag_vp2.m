% Calculate joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:07
% EndTime: 2019-12-05 18:55:10
% DurationCPUTime: 0.57s
% Computational Cost: add. (873->174), mult. (1859->265), div. (0->0), fcn. (1977->8), ass. (0->77)
t71 = sin(qJ(5));
t72 = sin(qJ(4));
t103 = t71 * t72;
t75 = cos(qJ(5));
t76 = cos(qJ(4));
t54 = t75 * t76 - t103;
t107 = t54 * mrSges(6,3);
t99 = Ifges(5,5) * t72 + Ifges(5,6) * t76;
t115 = t71 * pkin(3) * t107 + t99;
t73 = sin(qJ(3));
t74 = sin(qJ(2));
t77 = cos(qJ(3));
t78 = cos(qJ(2));
t57 = t73 * t78 + t77 * t74;
t104 = t57 * t76;
t55 = t73 * t74 - t77 * t78;
t114 = Ifges(5,5) * t104 + Ifges(5,3) * t55;
t102 = t75 * t72;
t56 = t71 * t76 + t102;
t32 = -t54 * mrSges(6,1) + t56 * mrSges(6,2);
t113 = 0.2e1 * t32;
t111 = m(4) * pkin(1) ^ 2;
t110 = t77 * pkin(1);
t109 = Ifges(5,4) * t72;
t108 = Ifges(5,4) * t76;
t106 = t56 * mrSges(6,3);
t105 = t57 * t72;
t29 = -t78 * pkin(1) + t55 * pkin(2) - t57 * pkin(5);
t101 = t76 * t29;
t100 = Ifges(6,5) * t56 + Ifges(6,6) * t54;
t69 = t72 ^ 2;
t98 = t76 ^ 2 + t69;
t97 = 0.2e1 * mrSges(6,3);
t96 = pkin(3) * t105;
t95 = t75 * t106;
t25 = t56 * t57;
t26 = t54 * t57;
t94 = Ifges(6,5) * t26 - Ifges(6,6) * t25 + Ifges(6,3) * t55;
t66 = -t76 * pkin(3) - pkin(2);
t64 = t73 * pkin(1) + pkin(5);
t93 = t98 * t64;
t33 = Ifges(6,4) * t56 + Ifges(6,2) * t54;
t34 = Ifges(6,1) * t56 + Ifges(6,4) * t54;
t60 = Ifges(5,2) * t76 + t109;
t61 = Ifges(5,1) * t72 + t108;
t92 = t54 * t33 + t56 * t34 + t76 * t60 + t72 * t61 + Ifges(4,3);
t91 = t76 * mrSges(5,1) - t72 * mrSges(5,2);
t90 = t72 * mrSges(5,1) + t76 * mrSges(5,2);
t30 = -t55 * mrSges(5,2) - mrSges(5,3) * t105;
t31 = t55 * mrSges(5,1) - mrSges(5,3) * t104;
t89 = t76 * t30 - t72 * t31;
t88 = 0.2e1 * t98 * mrSges(5,3);
t38 = t56 * t64;
t39 = t54 * t64;
t87 = -t38 * mrSges(6,1) - t39 * mrSges(6,2) + t100;
t44 = t56 * pkin(5);
t45 = t54 * pkin(5);
t86 = -t44 * mrSges(6,1) - t45 * mrSges(6,2) + t100;
t16 = t55 * pkin(3) + t101;
t5 = -t29 * t103 + t75 * t16;
t6 = t29 * t102 + t71 * t16;
t85 = t5 * mrSges(6,1) - t6 * mrSges(6,2) + t94;
t84 = (t77 * mrSges(4,1) - t73 * mrSges(4,2)) * pkin(1);
t83 = (t75 * mrSges(6,1) - t71 * mrSges(6,2)) * pkin(3);
t17 = Ifges(5,6) * t55 + (-Ifges(5,2) * t72 + t108) * t57;
t18 = Ifges(5,5) * t55 + (Ifges(5,1) * t76 - t109) * t57;
t7 = Ifges(6,4) * t26 - Ifges(6,2) * t25 + Ifges(6,6) * t55;
t8 = Ifges(6,1) * t26 - Ifges(6,4) * t25 + Ifges(6,5) * t55;
t82 = -t5 * t106 + t6 * t107 + t54 * t7 / 0.2e1 - t25 * t33 / 0.2e1 + t26 * t34 / 0.2e1 + t72 * t18 / 0.2e1 + t76 * t17 / 0.2e1 + t32 * t96 + t56 * t8 / 0.2e1 - t60 * t105 / 0.2e1 + t61 * t104 / 0.2e1 - Ifges(4,6) * t55 + Ifges(4,5) * t57 + (t100 + t99) * t55 / 0.2e1;
t80 = pkin(3) ^ 2;
t65 = -pkin(2) - t110;
t58 = t66 - t110;
t27 = t90 * t57;
t13 = t55 * mrSges(6,1) - t26 * mrSges(6,3);
t12 = -t55 * mrSges(6,2) - t25 * mrSges(6,3);
t9 = t25 * mrSges(6,1) + t26 * mrSges(6,2);
t1 = [Ifges(3,1) * t74 ^ 2 + 0.2e1 * t31 * t101 + 0.2e1 * t6 * t12 + 0.2e1 * t5 * t13 - t25 * t7 + t26 * t8 + Ifges(2,3) + (Ifges(4,1) * t57 + t76 * t18) * t57 + m(6) * (t69 * t57 ^ 2 * t80 + t5 ^ 2 + t6 ^ 2) + m(5) * t98 * t29 ^ 2 + (0.2e1 * t29 * t30 + (0.2e1 * pkin(3) * t9 - t17) * t57) * t72 + (Ifges(4,2) * t55 + (-Ifges(5,6) * t72 - (2 * Ifges(4,4))) * t57 + t94 + t114) * t55 + (-0.2e1 * pkin(1) * (t55 * mrSges(4,1) + t57 * mrSges(4,2)) + 0.2e1 * Ifges(3,4) * t74 + (Ifges(3,2) + t111) * t78) * t78; t82 + m(6) * (-t38 * t5 + t39 * t6 + t58 * t96) + t89 * t64 + (-t55 * t73 - t57 * t77) * pkin(1) * mrSges(4,3) + Ifges(3,5) * t74 + Ifges(3,6) * t78 + t65 * t27 + t58 * t9 + t39 * t12 - t38 * t13; t58 * t113 - 0.2e1 * t65 * t91 + Ifges(3,3) + 0.2e1 * t84 + (t38 * t56 + t39 * t54) * t97 + t64 * t88 + m(6) * (t38 ^ 2 + t39 ^ 2 + t58 ^ 2) + m(5) * (t98 * t64 ^ 2 + t65 ^ 2) + (t73 ^ 2 + t77 ^ 2) * t111 + t92; t82 + m(6) * (-t44 * t5 + t45 * t6 + t66 * t96) + t89 * pkin(5) + t66 * t9 - t44 * t13 + t45 * t12 - pkin(2) * t27; -(-pkin(2) + t65) * t91 + (t58 + t66) * t32 + t84 + m(5) * (-pkin(2) * t65 + pkin(5) * t93) + m(6) * (t44 * t38 + t45 * t39 + t66 * t58) + ((t38 + t44) * t56 + (t39 + t45) * t54) * mrSges(6,3) + (t98 * pkin(5) + t93) * mrSges(5,3) + t92; 0.2e1 * pkin(2) * t91 + t66 * t113 + (t44 * t56 + t45 * t54) * t97 + pkin(5) * t88 + m(6) * (t44 ^ 2 + t45 ^ 2 + t66 ^ 2) + m(5) * (t98 * pkin(5) ^ 2 + pkin(2) ^ 2) + t92; -Ifges(5,6) * t105 + t91 * t29 + (t75 * t13 + m(6) * (t5 * t75 + t6 * t71) + t71 * t12) * pkin(3) + t85 + t114; -t90 * t64 + (-t95 + m(6) * (-t38 * t75 + t39 * t71)) * pkin(3) + t87 + t115; -t90 * pkin(5) + (-t95 + m(6) * (-t44 * t75 + t45 * t71)) * pkin(3) + t86 + t115; Ifges(5,3) + Ifges(6,3) + m(6) * (t71 ^ 2 + t75 ^ 2) * t80 + 0.2e1 * t83; t85; t87; t86; Ifges(6,3) + t83; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
