% Calculate joint inertia matrix for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:06:41
% EndTime: 2018-11-23 16:06:41
% DurationCPUTime: 0.76s
% Computational Cost: add. (1567->195), mult. (2693->273), div. (0->0), fcn. (2968->8), ass. (0->83)
t75 = sin(qJ(6));
t78 = cos(qJ(6));
t104 = t75 ^ 2 + t78 ^ 2;
t133 = mrSges(7,3) * t104;
t56 = -mrSges(7,1) * t78 + mrSges(7,2) * t75;
t132 = t56 - mrSges(6,1);
t119 = cos(qJ(5));
t73 = sin(pkin(10));
t74 = cos(pkin(10));
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t51 = -t73 * t79 - t74 * t77;
t53 = -t73 * t77 + t74 * t79;
t76 = sin(qJ(5));
t36 = t119 * t53 + t51 * t76;
t121 = pkin(3) * t73;
t63 = pkin(3) * t74 + pkin(4);
t43 = t119 * t63 - t76 * t121;
t44 = t119 * t121 + t76 * t63;
t86 = -t119 * t51 + t76 * t53;
t131 = t43 * t36 + t44 * t86;
t130 = t79 ^ 2;
t129 = m(5) * pkin(3);
t109 = t75 * mrSges(7,3);
t15 = -mrSges(7,2) * t86 - t36 * t109;
t113 = t36 * t78;
t16 = mrSges(7,1) * t86 - mrSges(7,3) * t113;
t90 = t78 * t15 - t75 * t16;
t80 = -pkin(1) - pkin(7);
t102 = -qJ(4) + t80;
t93 = t102 * t79;
t94 = t102 * t77;
t38 = t73 * t93 + t74 * t94;
t21 = pkin(8) * t51 + t38;
t37 = -t73 * t94 + t74 * t93;
t85 = -t53 * pkin(8) + t37;
t10 = -t119 * t85 + t21 * t76;
t127 = t10 ^ 2;
t126 = t36 ^ 2;
t125 = t51 ^ 2;
t124 = 0.2e1 * t10;
t64 = t77 * pkin(3) + qJ(2);
t39 = -pkin(4) * t51 + t64;
t123 = 0.2e1 * t39;
t12 = t119 * t21 + t76 * t85;
t13 = pkin(5) * t86 - pkin(9) * t36 + t39;
t3 = t12 * t78 + t13 * t75;
t120 = t3 * t78;
t117 = Ifges(7,4) * t75;
t116 = Ifges(7,4) * t78;
t115 = t10 * t36;
t114 = t36 * t75;
t112 = t38 * t51;
t111 = t43 * mrSges(6,1);
t110 = t44 * mrSges(6,2);
t106 = Ifges(7,5) * t113 + Ifges(7,3) * t86;
t105 = Ifges(7,5) * t75 + Ifges(7,6) * t78;
t103 = t77 ^ 2 + t130;
t101 = t53 ^ 2 + t125;
t57 = Ifges(7,2) * t78 + t117;
t58 = Ifges(7,1) * t75 + t116;
t100 = t78 * t57 + t75 * t58 + Ifges(6,3);
t99 = m(4) * t103;
t98 = t104 * pkin(9);
t42 = pkin(9) + t44;
t97 = t104 * t42;
t96 = t103 * mrSges(4,3);
t95 = -t51 * mrSges(5,1) + t53 * mrSges(5,2);
t2 = -t12 * t75 + t13 * t78;
t92 = -t2 * t75 + t120;
t91 = mrSges(7,1) * t75 + mrSges(7,2) * t78;
t89 = t73 * t51 - t74 * t53;
t88 = 0.2e1 * t133;
t87 = -t132 * t36 + (-mrSges(6,2) + t133) * t86;
t6 = Ifges(7,6) * t86 + (-Ifges(7,2) * t75 + t116) * t36;
t7 = Ifges(7,5) * t86 + (Ifges(7,1) * t78 - t117) * t36;
t84 = -t12 * mrSges(6,2) + mrSges(7,3) * t120 - t2 * t109 - t57 * t114 / 0.2e1 + t58 * t113 / 0.2e1 + Ifges(6,5) * t36 + t75 * t7 / 0.2e1 + t78 * t6 / 0.2e1 + (t105 / 0.2e1 - Ifges(6,6)) * t86 + t132 * t10;
t81 = qJ(2) ^ 2;
t41 = -pkin(5) - t43;
t31 = t86 ^ 2;
t27 = t36 * mrSges(6,2);
t14 = t91 * t36;
t1 = [0.2e1 * mrSges(5,3) * t112 + Ifges(4,1) * t130 + 0.2e1 * t64 * t95 + Ifges(5,2) * t125 + t27 * t123 + t14 * t124 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 - (2 * pkin(1) * mrSges(3,2)) + Ifges(3,1) + Ifges(2,3) + (mrSges(6,1) * t123 - 0.2e1 * mrSges(6,3) * t12 + Ifges(6,2) * t86 + t106) * t86 - 0.2e1 * t80 * t96 + (-0.2e1 * mrSges(5,3) * t37 + Ifges(5,1) * t53 + 0.2e1 * Ifges(5,4) * t51) * t53 + (mrSges(6,3) * t124 + Ifges(6,1) * t36 - t75 * t6 + t78 * t7 + (-Ifges(7,6) * t75 - (2 * Ifges(6,4))) * t86) * t36 + m(4) * (t103 * t80 ^ 2 + t81) + m(3) * ((pkin(1) ^ 2) + t81) + m(5) * (t37 ^ 2 + t38 ^ 2 + t64 ^ 2) + m(6) * (t12 ^ 2 + t39 ^ 2 + t127) + m(7) * (t2 ^ 2 + t3 ^ 2 + t127) + (-0.2e1 * Ifges(4,4) * t79 + Ifges(4,2) * t77) * t77 + 0.2e1 * (mrSges(4,1) * t77 + mrSges(4,2) * t79 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) - (t36 * mrSges(6,3) + t14) * t36 - t101 * mrSges(5,3) - t96 + (-mrSges(6,3) * t86 + t90) * t86 + m(7) * (t86 * t92 - t115) + m(6) * (t12 * t86 - t115) + m(5) * (t37 * t53 - t112) + t80 * t99; m(3) + m(7) * (t104 * t31 + t126) + m(6) * (t31 + t126) + m(5) * t101 + t99; m(6) * (-t10 * t43 + t12 * t44) + (t80 * mrSges(4,1) + Ifges(4,5)) * t79 + (-t80 * mrSges(4,2) - Ifges(4,6)) * t77 - t131 * mrSges(6,3) + Ifges(5,6) * t51 + Ifges(5,5) * t53 + t37 * mrSges(5,1) - t38 * mrSges(5,2) + t41 * t14 + t84 + t90 * t42 + m(7) * (t10 * t41 + t92 * t42) + (m(5) * (t37 * t74 + t38 * t73) + t89 * mrSges(5,3)) * pkin(3); t79 * mrSges(4,1) + t53 * mrSges(5,1) - t77 * mrSges(4,2) + t51 * mrSges(5,2) + m(7) * (-t36 * t41 + t86 * t97) + m(6) * t131 - t89 * t129 + t87; 0.2e1 * t111 - 0.2e1 * t110 + 0.2e1 * t41 * t56 + Ifges(4,3) + Ifges(5,3) + t42 * t88 + m(7) * (t104 * t42 ^ 2 + t41 ^ 2) + m(6) * (t43 ^ 2 + t44 ^ 2) + t100 + (0.2e1 * mrSges(5,1) * t74 - 0.2e1 * mrSges(5,2) * t73 + (t73 ^ 2 + t74 ^ 2) * t129) * pkin(3); t86 * mrSges(6,1) + t75 * t15 + t78 * t16 + t27 + m(7) * (t2 * t78 + t3 * t75) + m(6) * t39 + m(5) * t64 + t95; 0; 0; m(7) * t104 + m(5) + m(6); t84 + (-m(7) * t10 - t14) * pkin(5) + (m(7) * t92 + t90) * pkin(9); m(7) * (pkin(5) * t36 + t86 * t98) + t87; m(7) * (-pkin(5) * t41 + pkin(9) * t97) - t110 + t111 + (t41 - pkin(5)) * t56 + (t97 + t98) * mrSges(7,3) + t100; 0; -0.2e1 * pkin(5) * t56 + m(7) * (t104 * pkin(9) ^ 2 + pkin(5) ^ 2) + pkin(9) * t88 + t100; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t114 + t106; -t91 * t86; -t91 * t42 + t105; -t56; -t91 * pkin(9) + t105; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
