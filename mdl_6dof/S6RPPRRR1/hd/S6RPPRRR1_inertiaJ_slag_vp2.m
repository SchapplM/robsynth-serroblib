% Calculate joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:47:51
% EndTime: 2018-11-23 15:47:52
% DurationCPUTime: 0.70s
% Computational Cost: add. (1468->169), mult. (2682->249), div. (0->0), fcn. (2964->10), ass. (0->76)
t69 = sin(qJ(6));
t72 = cos(qJ(6));
t92 = t69 ^ 2 + t72 ^ 2;
t117 = mrSges(7,3) * t92;
t65 = sin(pkin(11));
t67 = cos(pkin(11));
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t44 = -t71 * t65 + t74 * t67;
t45 = t74 * t65 + t71 * t67;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t36 = -t73 * t44 + t70 * t45;
t38 = t70 * t44 + t73 * t45;
t99 = t69 * mrSges(7,3);
t15 = -t36 * mrSges(7,2) - t38 * t99;
t100 = t38 * t72;
t16 = t36 * mrSges(7,1) - mrSges(7,3) * t100;
t116 = t72 * t15 - t69 * t16;
t49 = -t72 * mrSges(7,1) + t69 * mrSges(7,2);
t115 = t38 * t117 + t36 * t49;
t101 = t38 * t69;
t14 = -mrSges(7,1) * t101 - mrSges(7,2) * t100;
t66 = sin(pkin(10));
t53 = t66 * pkin(1) + qJ(3);
t105 = pkin(7) + t53;
t89 = t105 * t67;
t90 = t105 * t65;
t29 = -t71 * t90 + t74 * t89;
t18 = t44 * pkin(8) + t29;
t28 = -t71 * t89 - t74 * t90;
t79 = -t45 * pkin(8) + t28;
t8 = t70 * t18 - t73 * t79;
t114 = m(7) * t8 - t14;
t10 = t73 * t18 + t70 * t79;
t108 = pkin(5) * t36;
t68 = cos(pkin(10));
t55 = -t68 * pkin(1) - pkin(2);
t48 = -t67 * pkin(3) + t55;
t39 = -t44 * pkin(4) + t48;
t13 = -t38 * pkin(9) + t108 + t39;
t3 = t72 * t10 + t69 * t13;
t107 = t3 * t72;
t2 = -t69 * t10 + t72 * t13;
t84 = -t2 * t69 + t107;
t113 = m(7) * t84 + t116;
t112 = t8 ^ 2;
t111 = t36 ^ 2;
t110 = t44 ^ 2;
t62 = t67 ^ 2;
t109 = 0.2e1 * t44;
t106 = t8 * t36;
t103 = Ifges(7,4) * t69;
t102 = Ifges(7,4) * t72;
t96 = Ifges(7,5) * t100 + Ifges(7,3) * t36;
t95 = t36 * mrSges(6,1) + t38 * mrSges(6,2);
t94 = Ifges(7,5) * t69 + Ifges(7,6) * t72;
t93 = t65 ^ 2 + t62;
t50 = Ifges(7,2) * t72 + t103;
t51 = Ifges(7,1) * t69 + t102;
t91 = t72 * t50 + t69 * t51 + Ifges(6,3);
t88 = t92 * pkin(9);
t56 = t70 * pkin(4) + pkin(9);
t87 = t92 * t56;
t86 = -t67 * mrSges(4,1) + t65 * mrSges(4,2);
t85 = -t44 * mrSges(5,1) + t45 * mrSges(5,2);
t83 = -mrSges(7,1) * t69 - mrSges(7,2) * t72;
t82 = 0.2e1 * t117;
t81 = -t85 - t95;
t80 = (t73 * mrSges(6,1) - t70 * mrSges(6,2)) * pkin(4);
t11 = Ifges(7,6) * t36 + (-Ifges(7,2) * t69 + t102) * t38;
t12 = Ifges(7,5) * t36 + (Ifges(7,1) * t72 - t103) * t38;
t78 = -t10 * mrSges(6,2) + mrSges(7,3) * t107 - t2 * t99 - t50 * t101 / 0.2e1 + t51 * t100 / 0.2e1 + Ifges(6,5) * t38 + t69 * t12 / 0.2e1 + t72 * t11 / 0.2e1 + (t49 - mrSges(6,1)) * t8 + (t94 / 0.2e1 - Ifges(6,6)) * t36;
t57 = -t73 * pkin(4) - pkin(5);
t35 = t38 ^ 2;
t1 = [Ifges(4,2) * t62 + 0.2e1 * t55 * t86 + Ifges(5,2) * t110 + 0.2e1 * t48 * t85 + 0.2e1 * t39 * t95 - 0.2e1 * t8 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + t29 * mrSges(5,3) * t109 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t65 + 0.2e1 * Ifges(4,4) * t67) * t65 + (-0.2e1 * t28 * mrSges(5,3) + Ifges(5,1) * t45 + Ifges(5,4) * t109) * t45 + (-0.2e1 * t10 * mrSges(6,3) + Ifges(6,2) * t36 + t96) * t36 + (0.2e1 * t8 * mrSges(6,3) + Ifges(6,1) * t38 - t69 * t11 + t72 * t12 + (-Ifges(7,6) * t69 - (2 * Ifges(6,4))) * t36) * t38 + m(4) * (t93 * t53 ^ 2 + t55 ^ 2) + m(6) * (t10 ^ 2 + t39 ^ 2 + t112) + m(5) * (t28 ^ 2 + t29 ^ 2 + t48 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t112) + m(3) * (t66 ^ 2 + t68 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t68 * mrSges(3,1) - t66 * mrSges(3,2)) * pkin(1) + 0.2e1 * t93 * t53 * mrSges(4,3); -t36 * t14 + t116 * t38 + m(6) * (t10 * t38 + t106) + m(7) * (t84 * t38 + t106) + m(5) * (t28 * t44 + t29 * t45); m(3) + m(6) * (t35 + t111) + m(7) * (t92 * t35 + t111) + m(5) * (t45 ^ 2 + t110) + m(4) * t93; t69 * t15 + t72 * t16 + m(7) * (t72 * t2 + t69 * t3) + m(6) * t39 + m(5) * t48 + m(4) * t55 - t81 + t86; 0; m(7) * t92 + m(4) + m(5) + m(6); (m(6) * (t10 * t70 - t73 * t8) + (-t70 * t36 - t73 * t38) * mrSges(6,3)) * pkin(4) + t78 + Ifges(5,6) * t44 + Ifges(5,5) * t45 + t28 * mrSges(5,1) - t29 * mrSges(5,2) + t114 * t57 + t113 * t56; m(7) * (t57 * t36 + t38 * t87) + m(6) * (-t36 * t73 + t38 * t70) * pkin(4) + t81 + t115; 0; 0.2e1 * t57 * t49 + Ifges(5,3) + 0.2e1 * t80 + t56 * t82 + m(7) * (t92 * t56 ^ 2 + t57 ^ 2) + m(6) * (t70 ^ 2 + t73 ^ 2) * pkin(4) ^ 2 + t91; -t114 * pkin(5) + t113 * pkin(9) + t78; m(7) * (t38 * t88 - t108) - t95 + t115; 0; m(7) * (-pkin(5) * t57 + pkin(9) * t87) + (t57 - pkin(5)) * t49 + t80 + (t87 + t88) * mrSges(7,3) + t91; -0.2e1 * pkin(5) * t49 + m(7) * (t92 * pkin(9) ^ 2 + pkin(5) ^ 2) + pkin(9) * t82 + t91; t2 * mrSges(7,1) - t3 * mrSges(7,2) - Ifges(7,6) * t101 + t96; t14; -t49; t83 * t56 + t94; t83 * pkin(9) + t94; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
