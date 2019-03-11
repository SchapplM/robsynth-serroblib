% Calculate joint inertia matrix for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:09
% EndTime: 2019-03-09 02:33:10
% DurationCPUTime: 0.73s
% Computational Cost: add. (1468->176), mult. (2527->248), div. (0->0), fcn. (2821->8), ass. (0->79)
t73 = sin(qJ(6));
t75 = cos(qJ(6));
t97 = t73 ^ 2 + t75 ^ 2;
t128 = mrSges(7,3) * t97;
t51 = -t75 * mrSges(7,1) + t73 * mrSges(7,2);
t127 = t51 - mrSges(6,1);
t112 = sin(qJ(4));
t113 = cos(qJ(4));
t70 = sin(pkin(10));
t71 = cos(pkin(10));
t46 = -t112 * t71 - t113 * t70;
t48 = -t112 * t70 + t113 * t71;
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t36 = t46 * t74 + t76 * t48;
t85 = -t76 * t46 + t48 * t74;
t126 = t36 * t76 + t74 * t85;
t104 = t73 * mrSges(7,3);
t15 = -mrSges(7,2) * t85 - t36 * t104;
t106 = t36 * t75;
t16 = mrSges(7,1) * t85 - mrSges(7,3) * t106;
t124 = t75 * t15 - t73 * t16;
t72 = -pkin(1) - qJ(3);
t114 = -pkin(7) + t72;
t88 = t114 * t112;
t89 = t114 * t113;
t38 = t70 * t89 + t71 * t88;
t20 = pkin(8) * t46 + t38;
t37 = -t70 * t88 + t71 * t89;
t80 = -t48 * pkin(8) + t37;
t10 = t74 * t20 - t76 * t80;
t86 = mrSges(7,1) * t73 + mrSges(7,2) * t75;
t14 = t86 * t36;
t123 = m(7) * t10 + t14;
t12 = t76 * t20 + t74 * t80;
t54 = t70 * pkin(3) + qJ(2);
t39 = -t46 * pkin(4) + t54;
t13 = pkin(5) * t85 - pkin(9) * t36 + t39;
t3 = t12 * t75 + t13 * t73;
t115 = t3 * t75;
t2 = -t12 * t73 + t13 * t75;
t87 = -t2 * t73 + t115;
t122 = m(7) * t87 + t124;
t121 = t10 ^ 2;
t120 = t36 ^ 2;
t119 = t46 ^ 2;
t66 = t71 ^ 2;
t118 = 0.2e1 * t10;
t117 = 0.2e1 * t39;
t110 = Ifges(7,4) * t73;
t109 = Ifges(7,4) * t75;
t108 = t36 * t10;
t107 = t36 * t73;
t105 = t38 * t46;
t101 = Ifges(7,5) * t106 + Ifges(7,3) * t85;
t100 = t70 * mrSges(4,1) + t71 * mrSges(4,2);
t99 = Ifges(7,5) * t73 + Ifges(7,6) * t75;
t98 = t70 ^ 2 + t66;
t96 = t48 ^ 2 + t119;
t52 = Ifges(7,2) * t75 + t110;
t53 = Ifges(7,1) * t73 + t109;
t95 = t75 * t52 + t73 * t53 + Ifges(6,3);
t94 = m(4) * t98;
t93 = t97 * pkin(9);
t92 = t98 * mrSges(4,3);
t58 = t74 * pkin(4) + pkin(9);
t91 = t97 * t58;
t90 = -t46 * mrSges(5,1) + t48 * mrSges(5,2);
t84 = 0.2e1 * t128;
t83 = -t127 * t36 + (-mrSges(6,2) + t128) * t85;
t82 = (mrSges(6,1) * t76 - mrSges(6,2) * t74) * pkin(4);
t6 = Ifges(7,6) * t85 + (-Ifges(7,2) * t73 + t109) * t36;
t7 = Ifges(7,5) * t85 + (Ifges(7,1) * t75 - t110) * t36;
t81 = -t12 * mrSges(6,2) + mrSges(7,3) * t115 - t2 * t104 - t52 * t107 / 0.2e1 + t53 * t106 / 0.2e1 + Ifges(6,5) * t36 + t73 * t7 / 0.2e1 + t75 * t6 / 0.2e1 + (t99 / 0.2e1 - Ifges(6,6)) * t85 + t127 * t10;
t77 = qJ(2) ^ 2;
t59 = -t76 * pkin(4) - pkin(5);
t31 = t85 ^ 2;
t27 = t36 * mrSges(6,2);
t1 = [0.2e1 * mrSges(5,3) * t105 + Ifges(4,1) * t66 + 0.2e1 * t54 * t90 + Ifges(5,2) * t119 + t27 * t117 + 0.2e1 * t2 * t16 + t14 * t118 + 0.2e1 * t3 * t15 - (2 * pkin(1) * mrSges(3,2)) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t71 + Ifges(4,2) * t70) * t70 + (mrSges(6,1) * t117 - 0.2e1 * mrSges(6,3) * t12 + Ifges(6,2) * t85 + t101) * t85 - 0.2e1 * t72 * t92 + (-0.2e1 * mrSges(5,3) * t37 + Ifges(5,1) * t48 + 0.2e1 * Ifges(5,4) * t46) * t48 + (mrSges(6,3) * t118 + Ifges(6,1) * t36 - t73 * t6 + t75 * t7 + (-Ifges(7,6) * t73 - (2 * Ifges(6,4))) * t85) * t36 + m(4) * (t98 * t72 ^ 2 + t77) + m(3) * ((pkin(1) ^ 2) + t77) + m(5) * (t37 ^ 2 + t38 ^ 2 + t54 ^ 2) + m(6) * (t12 ^ 2 + t39 ^ 2 + t121) + m(7) * (t2 ^ 2 + t3 ^ 2 + t121) + 0.2e1 * (t100 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) - (t36 * mrSges(6,3) + t14) * t36 - t96 * mrSges(5,3) - t92 + (-mrSges(6,3) * t85 + t124) * t85 + m(7) * (t85 * t87 - t108) + m(6) * (t12 * t85 - t108) + m(5) * (t37 * t48 - t105) + t72 * t94; m(3) + m(7) * (t97 * t31 + t120) + m(6) * (t31 + t120) + m(5) * t96 + t94; m(4) * qJ(2) + t85 * mrSges(6,1) + t73 * t15 + t75 * t16 + t27 + m(7) * (t2 * t75 + t3 * t73) + m(6) * t39 + m(5) * t54 + t90 + t100; 0; m(7) * t97 + m(4) + m(5) + m(6); (m(6) * (-t10 * t76 + t12 * t74) - t126 * mrSges(6,3)) * pkin(4) + t81 + Ifges(5,6) * t46 + Ifges(5,5) * t48 + t37 * mrSges(5,1) - t38 * mrSges(5,2) + t123 * t59 + t122 * t58; t48 * mrSges(5,1) + t46 * mrSges(5,2) + m(7) * (-t36 * t59 + t85 * t91) + m(6) * t126 * pkin(4) + t83; 0; 0.2e1 * t59 * t51 + Ifges(5,3) + 0.2e1 * t82 + t58 * t84 + m(7) * (t97 * t58 ^ 2 + t59 ^ 2) + m(6) * (t74 ^ 2 + t76 ^ 2) * pkin(4) ^ 2 + t95; -t123 * pkin(5) + t122 * pkin(9) + t81; m(7) * (pkin(5) * t36 + t85 * t93) + t83; 0; m(7) * (-pkin(5) * t59 + pkin(9) * t91) + (t59 - pkin(5)) * t51 + t82 + (t91 + t93) * mrSges(7,3) + t95; -0.2e1 * pkin(5) * t51 + m(7) * (t97 * pkin(9) ^ 2 + pkin(5) ^ 2) + pkin(9) * t84 + t95; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t107 + t101; -t86 * t85; -t51; -t86 * t58 + t99; -t86 * pkin(9) + t99; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
