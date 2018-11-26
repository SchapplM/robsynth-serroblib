% Calculate joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2018-11-23 15:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:43:53
% EndTime: 2018-11-23 15:43:54
% DurationCPUTime: 0.80s
% Computational Cost: add. (928->192), mult. (1734->261), div. (0->0), fcn. (1716->8), ass. (0->84)
t114 = Ifges(6,5) + Ifges(7,5);
t113 = Ifges(6,6) + Ifges(7,6);
t112 = Ifges(6,3) + Ifges(7,3);
t111 = -2 * mrSges(7,3);
t110 = m(6) + m(7);
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t109 = t113 * t68 + t114 * t66;
t63 = sin(pkin(9));
t48 = t63 * pkin(1) + qJ(3);
t100 = pkin(7) + t48;
t64 = cos(pkin(10));
t29 = t100 * t64;
t67 = sin(qJ(4));
t62 = sin(pkin(10));
t99 = cos(qJ(4));
t75 = t99 * t62;
t15 = t100 * t75 + t67 * t29;
t108 = t15 ^ 2;
t88 = t67 * t62;
t34 = -t99 * t64 + t88;
t107 = t34 ^ 2;
t59 = t64 ^ 2;
t106 = 0.2e1 * t15;
t65 = cos(pkin(9));
t50 = -t65 * pkin(1) - pkin(2);
t37 = -t64 * pkin(3) + t50;
t105 = 0.2e1 * t37;
t104 = m(7) * pkin(5);
t102 = pkin(4) * t34;
t36 = t67 * t64 + t75;
t101 = t36 * pkin(8);
t98 = mrSges(6,2) * t68;
t97 = Ifges(6,4) * t66;
t96 = Ifges(6,4) * t68;
t95 = Ifges(7,4) * t66;
t94 = Ifges(7,4) * t68;
t93 = t15 * t34;
t92 = t36 * t66;
t91 = t36 * t68;
t90 = t66 * mrSges(6,2);
t89 = t66 * mrSges(7,3);
t86 = -qJ(6) - pkin(8);
t14 = -t101 + t37 + t102;
t17 = -t100 * t88 + t99 * t29;
t4 = t66 * t14 + t68 * t17;
t20 = -t34 * mrSges(7,2) - t36 * t89;
t21 = -t34 * mrSges(6,2) - mrSges(6,3) * t92;
t85 = t20 + t21;
t22 = t34 * mrSges(7,1) - mrSges(7,3) * t91;
t23 = t34 * mrSges(6,1) - mrSges(6,3) * t91;
t84 = t22 + t23;
t18 = mrSges(7,1) * t92 + mrSges(7,2) * t91;
t40 = -t68 * mrSges(6,1) + t90;
t83 = t40 - mrSges(5,1);
t80 = t62 ^ 2 + t59;
t79 = t66 ^ 2 + t68 ^ 2;
t78 = qJ(6) * t36;
t76 = -mrSges(6,1) - t104;
t74 = t79 * mrSges(6,3);
t73 = -t64 * mrSges(4,1) + t62 * mrSges(4,2);
t53 = t66 * mrSges(7,2);
t39 = -t68 * mrSges(7,1) + t53;
t3 = t68 * t14 - t66 * t17;
t72 = t112 * t34 + t114 * t91;
t71 = mrSges(6,1) * t66 + t98;
t51 = -t68 * pkin(5) - pkin(4);
t45 = Ifges(6,1) * t66 + t96;
t44 = Ifges(7,1) * t66 + t94;
t43 = Ifges(6,2) * t68 + t97;
t42 = Ifges(7,2) * t68 + t95;
t41 = t86 * t68;
t38 = t86 * t66;
t33 = t36 ^ 2;
t30 = t36 * mrSges(5,2);
t19 = t71 * t36;
t9 = Ifges(6,5) * t34 + (Ifges(6,1) * t68 - t97) * t36;
t8 = Ifges(7,5) * t34 + (Ifges(7,1) * t68 - t95) * t36;
t7 = Ifges(6,6) * t34 + (-Ifges(6,2) * t66 + t96) * t36;
t6 = Ifges(7,6) * t34 + (-Ifges(7,2) * t66 + t94) * t36;
t5 = pkin(5) * t92 + t15;
t2 = -t66 * t78 + t4;
t1 = t34 * pkin(5) - t68 * t78 + t3;
t10 = [Ifges(4,2) * t59 + 0.2e1 * t50 * t73 + t30 * t105 + 0.2e1 * t5 * t18 + t19 * t106 + 0.2e1 * t2 * t20 + 0.2e1 * t4 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t3 * t23 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t62 + 0.2e1 * Ifges(4,4) * t64) * t62 + (mrSges(5,1) * t105 - 0.2e1 * t17 * mrSges(5,3) + Ifges(5,2) * t34 + t72) * t34 + (mrSges(5,3) * t106 + Ifges(5,1) * t36 - 0.2e1 * Ifges(5,4) * t34 + (t8 + t9) * t68 + (-t113 * t34 - t6 - t7) * t66) * t36 + m(4) * (t80 * t48 ^ 2 + t50 ^ 2) + m(5) * (t17 ^ 2 + t37 ^ 2 + t108) + m(6) * (t3 ^ 2 + t4 ^ 2 + t108) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (t63 ^ 2 + t65 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t65 * mrSges(3,1) - t63 * mrSges(3,2)) * pkin(1) + 0.2e1 * t80 * t48 * mrSges(4,3); (t18 + t19) * t34 + (-t84 * t66 + t85 * t68) * t36 + m(6) * (t93 + (-t3 * t66 + t4 * t68) * t36) + m(7) * (t5 * t34 + (-t1 * t66 + t2 * t68) * t36) + m(5) * (t17 * t36 + t93); m(3) + m(5) * (t33 + t107) + m(4) * t80 + (t79 * t33 + t107) * t110; t34 * mrSges(5,1) + t30 + t84 * t68 + t85 * t66 + m(7) * (t68 * t1 + t66 * t2) + m(6) * (t68 * t3 + t66 * t4) + m(5) * t37 + m(4) * t50 + t73; 0; t79 * t110 + m(4) + m(5); m(7) * (t38 * t1 - t41 * t2 + t51 * t5) + t51 * t18 + t38 * t22 + t5 * t39 - t41 * t20 - t17 * mrSges(5,2) - pkin(4) * t19 + (-m(6) * pkin(4) + t83) * t15 + (t4 * mrSges(6,3) + t2 * mrSges(7,3) + t6 / 0.2e1 + t7 / 0.2e1 + (m(6) * t4 + t21) * pkin(8)) * t68 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + t8 / 0.2e1 + t9 / 0.2e1 + (-m(6) * t3 - t23) * pkin(8)) * t66 + (Ifges(5,5) + (t44 / 0.2e1 + t45 / 0.2e1) * t68 + (-t42 / 0.2e1 - t43 / 0.2e1) * t66) * t36 + (-Ifges(5,6) + t109 / 0.2e1) * t34; -t30 + (t39 + t83) * t34 + (t79 * mrSges(7,3) + t74) * t36 + m(6) * (t79 * t101 - t102) + m(7) * (t51 * t34 + (-t38 * t66 - t41 * t68) * t36); m(7) * (t38 * t68 - t41 * t66); -0.2e1 * pkin(4) * t40 + 0.2e1 * t51 * t39 + Ifges(5,3) + 0.2e1 * pkin(8) * t74 + m(6) * (t79 * pkin(8) ^ 2 + pkin(4) ^ 2) + m(7) * (t38 ^ 2 + t41 ^ 2 + t51 ^ 2) + (t111 * t41 + t42 + t43) * t68 + (t111 * t38 + t44 + t45) * t66; t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) - t113 * t92 + (m(7) * t1 + t22) * pkin(5) + t72; (t76 * t66 - t98) * t36 - t18; -t90 - t53 + (mrSges(7,1) - t76) * t68; t38 * mrSges(7,1) + t41 * mrSges(7,2) - t71 * pkin(8) + (m(7) * t38 - t89) * pkin(5) + t109; (0.2e1 * mrSges(7,1) + t104) * pkin(5) + t112; m(7) * t5 + t18; m(7) * t34; 0; m(7) * t51 + t39; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
