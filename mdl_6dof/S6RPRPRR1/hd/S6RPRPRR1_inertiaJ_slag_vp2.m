% Calculate joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:54
% EndTime: 2019-03-09 03:33:55
% DurationCPUTime: 0.84s
% Computational Cost: add. (1566->188), mult. (2850->273), div. (0->0), fcn. (3116->10), ass. (0->81)
t71 = sin(qJ(6));
t74 = cos(qJ(6));
t95 = t71 ^ 2 + t74 ^ 2;
t122 = mrSges(7,3) * t95;
t75 = cos(qJ(3));
t121 = t75 ^ 2;
t120 = m(5) * pkin(3);
t101 = t71 * mrSges(7,3);
t110 = cos(qJ(5));
t67 = sin(pkin(11));
t69 = cos(pkin(11));
t73 = sin(qJ(3));
t47 = -t67 * t73 + t69 * t75;
t48 = t67 * t75 + t69 * t73;
t72 = sin(qJ(5));
t34 = -t110 * t47 + t48 * t72;
t36 = t110 * t48 + t72 * t47;
t15 = -mrSges(7,2) * t34 - t36 * t101;
t104 = t36 * t74;
t16 = mrSges(7,1) * t34 - mrSges(7,3) * t104;
t83 = t74 * t15 - t71 * t16;
t52 = -mrSges(7,1) * t74 + mrSges(7,2) * t71;
t119 = t36 * t122 + t34 * t52;
t68 = sin(pkin(10));
t58 = pkin(1) * t68 + pkin(7);
t93 = qJ(4) + t58;
t87 = t93 * t75;
t88 = t93 * t73;
t27 = -t67 * t88 + t69 * t87;
t18 = pkin(8) * t47 + t27;
t26 = -t67 * t87 - t69 * t88;
t80 = -t48 * pkin(8) + t26;
t10 = -t110 * t80 + t18 * t72;
t118 = t10 ^ 2;
t117 = t34 ^ 2;
t116 = t47 ^ 2;
t115 = 0.2e1 * t10;
t114 = 0.2e1 * t47;
t113 = pkin(3) * t67;
t112 = pkin(5) * t34;
t12 = t110 * t18 + t72 * t80;
t70 = cos(pkin(10));
t60 = -pkin(1) * t70 - pkin(2);
t51 = -pkin(3) * t75 + t60;
t37 = -pkin(4) * t47 + t51;
t13 = -pkin(9) * t36 + t112 + t37;
t3 = t12 * t74 + t13 * t71;
t111 = t3 * t74;
t108 = Ifges(7,4) * t71;
t107 = Ifges(7,4) * t74;
t106 = t10 * t34;
t105 = t36 * t71;
t59 = pkin(3) * t69 + pkin(4);
t43 = t110 * t59 - t72 * t113;
t103 = t43 * mrSges(6,1);
t44 = t110 * t113 + t72 * t59;
t102 = t44 * mrSges(6,2);
t98 = Ifges(7,5) * t104 + Ifges(7,3) * t34;
t97 = t34 * mrSges(6,1) + t36 * mrSges(6,2);
t96 = Ifges(7,5) * t71 + Ifges(7,6) * t74;
t94 = t73 ^ 2 + t121;
t53 = Ifges(7,2) * t74 + t108;
t54 = Ifges(7,1) * t71 + t107;
t92 = t74 * t53 + t71 * t54 + Ifges(6,3);
t91 = t95 * pkin(9);
t42 = pkin(9) + t44;
t90 = t95 * t42;
t89 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t2 = -t12 * t71 + t13 * t74;
t86 = -t2 * t71 + t111;
t85 = -t75 * mrSges(4,1) + t73 * mrSges(4,2);
t84 = mrSges(7,1) * t71 + mrSges(7,2) * t74;
t82 = 0.2e1 * t122;
t81 = t89 + t97;
t8 = Ifges(7,6) * t34 + (-Ifges(7,2) * t71 + t107) * t36;
t9 = Ifges(7,5) * t34 + (Ifges(7,1) * t74 - t108) * t36;
t79 = -t12 * mrSges(6,2) + mrSges(7,3) * t111 - t2 * t101 - t53 * t105 / 0.2e1 + t54 * t104 / 0.2e1 + Ifges(6,5) * t36 + t71 * t9 / 0.2e1 + t74 * t8 / 0.2e1 + (t96 / 0.2e1 - Ifges(6,6)) * t34 + (t52 - mrSges(6,1)) * t10;
t41 = -pkin(5) - t43;
t33 = t36 ^ 2;
t14 = t84 * t36;
t1 = [0.2e1 * t60 * t85 + Ifges(4,2) * t121 + 0.2e1 * t51 * t89 + Ifges(5,2) * t116 + 0.2e1 * t37 * t97 + 0.2e1 * t2 * t16 + t27 * mrSges(5,3) * t114 + t14 * t115 + 0.2e1 * t3 * t15 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * mrSges(5,3) * t26 + Ifges(5,1) * t48 + Ifges(5,4) * t114) * t48 + (-0.2e1 * mrSges(6,3) * t12 + Ifges(6,2) * t34 + t98) * t34 + (mrSges(6,3) * t115 + Ifges(6,1) * t36 - t71 * t8 + t74 * t9 + (-Ifges(7,6) * t71 - (2 * Ifges(6,4))) * t34) * t36 + m(4) * (t94 * t58 ^ 2 + t60 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2 + t51 ^ 2) + m(6) * (t12 ^ 2 + t37 ^ 2 + t118) + m(7) * (t2 ^ 2 + t3 ^ 2 + t118) + m(3) * (t68 ^ 2 + t70 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t73 + 0.2e1 * Ifges(4,4) * t75) * t73 + 0.2e1 * (mrSges(3,1) * t70 - mrSges(3,2) * t68) * pkin(1) + 0.2e1 * t94 * t58 * mrSges(4,3); t34 * t14 + t83 * t36 + m(6) * (t12 * t36 + t106) + m(7) * (t86 * t36 + t106) + m(5) * (t26 * t47 + t27 * t48); m(3) + m(6) * (t33 + t117) + m(7) * (t95 * t33 + t117) + m(5) * (t48 ^ 2 + t116) + m(4) * t94; t83 * t42 + m(6) * (-t10 * t43 + t12 * t44) + (-mrSges(4,1) * t73 - mrSges(4,2) * t75) * t58 + (-t34 * t44 - t36 * t43) * mrSges(6,3) + Ifges(4,5) * t73 + Ifges(4,6) * t75 + t41 * t14 + Ifges(5,6) * t47 + Ifges(5,5) * t48 + t26 * mrSges(5,1) - t27 * mrSges(5,2) + t79 + m(7) * (t10 * t41 + t86 * t42) + (m(5) * (t26 * t69 + t27 * t67) + (t47 * t67 - t48 * t69) * mrSges(5,3)) * pkin(3); m(6) * (-t34 * t43 + t36 * t44) + m(7) * (t34 * t41 + t36 * t90) + (t47 * t69 + t48 * t67) * t120 - t81 - t85 + t119; 0.2e1 * t103 - 0.2e1 * t102 + 0.2e1 * t41 * t52 + Ifges(4,3) + Ifges(5,3) + t42 * t82 + m(7) * (t95 * t42 ^ 2 + t41 ^ 2) + m(6) * (t43 ^ 2 + t44 ^ 2) + t92 + (0.2e1 * mrSges(5,1) * t69 - 0.2e1 * mrSges(5,2) * t67 + (t67 ^ 2 + t69 ^ 2) * t120) * pkin(3); t71 * t15 + t74 * t16 + m(7) * (t2 * t74 + t3 * t71) + m(6) * t37 + m(5) * t51 + t81; 0; 0; m(7) * t95 + m(5) + m(6); t79 + (-m(7) * t10 - t14) * pkin(5) + (m(7) * t86 + t83) * pkin(9); m(7) * (t36 * t91 - t112) - t97 + t119; m(7) * (-pkin(5) * t41 + pkin(9) * t90) - t102 + t103 + (t41 - pkin(5)) * t52 + (t90 + t91) * mrSges(7,3) + t92; 0; -0.2e1 * pkin(5) * t52 + m(7) * (t95 * pkin(9) ^ 2 + pkin(5) ^ 2) + pkin(9) * t82 + t92; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t105 + t98; -t14; -t84 * t42 + t96; -t52; -t84 * pkin(9) + t96; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
