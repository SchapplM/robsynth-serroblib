% Calculate joint inertia matrix for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:18:12
% EndTime: 2018-11-23 16:18:13
% DurationCPUTime: 0.66s
% Computational Cost: add. (1595->189), mult. (2938->252), div. (0->0), fcn. (3306->8), ass. (0->81)
t115 = mrSges(6,2) - mrSges(5,1);
t69 = sin(qJ(6));
t72 = cos(qJ(6));
t99 = t69 ^ 2 + t72 ^ 2;
t94 = t99 * mrSges(7,3);
t51 = mrSges(7,1) * t69 + mrSges(7,2) * t72;
t114 = mrSges(6,3) + t51;
t70 = sin(qJ(4));
t56 = pkin(3) * t70 + qJ(5);
t113 = t56 ^ 2;
t68 = cos(pkin(10));
t63 = t68 ^ 2;
t67 = sin(pkin(10));
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t45 = -t67 * t71 + t68 * t74;
t46 = t67 * t74 + t68 * t71;
t73 = cos(qJ(4));
t35 = -t73 * t45 + t46 * t70;
t36 = t45 * t70 + t46 * t73;
t58 = -pkin(2) * t68 - pkin(1);
t39 = -pkin(3) * t45 + t58;
t80 = -qJ(5) * t36 + t39;
t13 = pkin(4) * t35 + t80;
t112 = -0.2e1 * t13;
t111 = 0.2e1 * t39;
t110 = 0.2e1 * t45;
t109 = pkin(4) + pkin(9);
t48 = m(7) * t99;
t108 = m(6) + t48;
t107 = mrSges(7,1) * t72;
t106 = Ifges(7,4) * t69;
t105 = Ifges(7,4) * t72;
t104 = t35 * t69;
t103 = t35 * t72;
t102 = mrSges(6,1) + mrSges(5,3);
t101 = pkin(7) + qJ(2);
t49 = t101 * t67;
t50 = t101 * t68;
t38 = -t71 * t49 + t74 * t50;
t100 = t67 ^ 2 + t63;
t98 = qJ(5) * t56;
t96 = Ifges(7,5) * t104 + Ifges(7,6) * t103 + Ifges(7,3) * t36;
t25 = pkin(8) * t45 + t38;
t37 = -t74 * t49 - t50 * t71;
t82 = -pkin(8) * t46 + t37;
t14 = t25 * t70 - t73 * t82;
t16 = t73 * t25 + t70 * t82;
t95 = t14 ^ 2 + t16 ^ 2;
t59 = -pkin(3) * t73 - pkin(4);
t93 = t99 * t109;
t92 = -t68 * mrSges(3,1) + t67 * mrSges(3,2);
t91 = -t45 * mrSges(4,1) + t46 * mrSges(4,2);
t61 = Ifges(7,5) * t72;
t90 = -Ifges(7,6) * t69 + t61;
t89 = 0.2e1 * t114;
t4 = t109 * t35 + t80;
t5 = pkin(5) * t36 + t14;
t1 = -t4 * t69 + t5 * t72;
t2 = t4 * t72 + t5 * t69;
t88 = t1 * t72 + t2 * t69;
t87 = mrSges(6,2) - t94;
t86 = -t69 * mrSges(7,2) + t107;
t19 = mrSges(7,1) * t36 - mrSges(7,3) * t104;
t20 = -mrSges(7,2) * t36 + mrSges(7,3) * t103;
t85 = t72 * t19 + t69 * t20;
t84 = -0.2e1 * t94;
t52 = -Ifges(7,2) * t69 + t105;
t53 = Ifges(7,1) * t72 - t106;
t83 = -t69 * t52 + t72 * t53 + Ifges(6,1) + Ifges(5,3);
t81 = (mrSges(5,1) * t73 - mrSges(5,2) * t70) * pkin(3);
t10 = Ifges(7,5) * t36 + (Ifges(7,1) * t69 + t105) * t35;
t6 = -t35 * pkin(5) + t16;
t9 = Ifges(7,6) * t36 + (Ifges(7,2) * t72 + t106) * t35;
t79 = t53 * t104 / 0.2e1 + t52 * t103 / 0.2e1 + t6 * t51 - t69 * t9 / 0.2e1 + t72 * t10 / 0.2e1 + (-mrSges(5,2) + mrSges(6,3)) * t16 - t88 * mrSges(7,3) + (Ifges(6,5) - Ifges(5,6)) * t35 + t115 * t14 + (t90 / 0.2e1 + Ifges(5,5) - Ifges(6,4)) * t36;
t76 = qJ(5) ^ 2;
t55 = -pkin(9) + t59;
t29 = t36 * mrSges(6,3);
t28 = t36 * mrSges(5,2);
t18 = t86 * t35;
t3 = [Ifges(3,2) * t63 - 0.2e1 * pkin(1) * t92 + 0.2e1 * t58 * t91 + Ifges(4,2) * t45 ^ 2 + t28 * t111 + t29 * t112 - 0.2e1 * t6 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + t38 * mrSges(4,3) * t110 + Ifges(2,3) + (Ifges(3,1) * t67 + 0.2e1 * Ifges(3,4) * t68) * t67 + 0.2e1 * t100 * qJ(2) * mrSges(3,3) + (-0.2e1 * t37 * mrSges(4,3) + Ifges(4,1) * t46 + Ifges(4,4) * t110) * t46 + (mrSges(5,1) * t111 + mrSges(6,2) * t112 + t69 * t10 + t72 * t9 + (Ifges(5,2) + Ifges(6,3)) * t35 - 0.2e1 * t102 * t16) * t35 + m(3) * (t100 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2 + t58 ^ 2) + m(5) * (t39 ^ 2 + t95) + m(6) * (t13 ^ 2 + t95) + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) + ((Ifges(5,1) + Ifges(6,2)) * t36 + t96 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t35 + 0.2e1 * t102 * t14) * t36; -m(3) * pkin(1) - t69 * t19 + t72 * t20 + t28 - t29 - t115 * t35 + m(7) * (-t1 * t69 + t2 * t72) + m(6) * t13 + m(5) * t39 + m(4) * t58 + t91 + t92; m(3) + m(4) + m(5) + t108; t85 * t55 + m(6) * (t14 * t59 + t16 * t56) + (-t35 * t56 + t36 * t59) * mrSges(6,1) + (m(5) * (-t14 * t73 + t16 * t70) + (-t35 * t70 - t36 * t73) * mrSges(5,3)) * pkin(3) + m(7) * (t88 * t55 + t56 * t6) - t56 * t18 + Ifges(4,6) * t45 + Ifges(4,5) * t46 + t37 * mrSges(4,1) - t38 * mrSges(4,2) + t79; 0; 0.2e1 * t59 * mrSges(6,2) + Ifges(4,3) + t56 * t89 + 0.2e1 * t81 + t55 * t84 + m(7) * (t99 * t55 ^ 2 + t113) + m(6) * (t59 ^ 2 + t113) + m(5) * (t70 ^ 2 + t73 ^ 2) * pkin(3) ^ 2 + t83; -t85 * t109 + m(6) * (-pkin(4) * t14 + qJ(5) * t16) + (-pkin(4) * t36 - qJ(5) * t35) * mrSges(6,1) + m(7) * (qJ(5) * t6 - t109 * t88) - qJ(5) * t18 + t79; 0; t81 + (-pkin(4) + t59) * mrSges(6,2) + m(7) * (-t55 * t93 + t98) + m(6) * (-pkin(4) * t59 + t98) + t83 + (-t55 + t109) * t94 + t114 * (qJ(5) + t56); -0.2e1 * pkin(4) * mrSges(6,2) + qJ(5) * t89 - t109 * t84 + m(6) * (pkin(4) ^ 2 + t76) + m(7) * (t109 ^ 2 * t99 + t76) + t83; m(6) * t14 + m(7) * t88 + t36 * mrSges(6,1) + t85; 0; m(6) * t59 + t55 * t48 + t87; -m(6) * pkin(4) - m(7) * t93 + t87; t108; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t96; -t51; t86 * t55 + t90; -t109 * t107 + t61 + (mrSges(7,2) * t109 - Ifges(7,6)) * t69; t86; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
