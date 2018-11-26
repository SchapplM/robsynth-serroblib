% Calculate joint inertia matrix for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2018-11-23 16:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:15:09
% EndTime: 2018-11-23 16:15:10
% DurationCPUTime: 0.85s
% Computational Cost: add. (651->259), mult. (1227->320), div. (0->0), fcn. (867->4), ass. (0->95)
t120 = Ifges(7,4) + Ifges(6,5);
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t90 = t70 ^ 2 + t72 ^ 2;
t119 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t30 = -mrSges(5,1) * t72 + mrSges(5,2) * t70;
t118 = m(5) * pkin(3) + mrSges(4,1) - t30;
t74 = (-pkin(1) - pkin(7));
t117 = -2 * t74;
t116 = 2 * mrSges(7,1);
t115 = 2 * qJ(2);
t114 = m(6) + m(7);
t113 = pkin(5) + pkin(8);
t111 = Ifges(5,4) * t70;
t110 = Ifges(5,4) * t72;
t71 = sin(qJ(3));
t109 = Ifges(6,4) * t71;
t108 = Ifges(5,6) * t71;
t107 = Ifges(6,6) * t70;
t106 = Ifges(6,6) * t72;
t105 = Ifges(7,6) * t70;
t104 = Ifges(7,6) * t72;
t103 = t70 * t71;
t73 = cos(qJ(3));
t102 = t70 * t73;
t101 = t71 * t72;
t100 = t71 * t74;
t99 = t72 * t73;
t97 = mrSges(6,1) + mrSges(7,1);
t96 = mrSges(6,2) - mrSges(7,3);
t95 = mrSges(6,3) + mrSges(7,2);
t69 = -pkin(4) - qJ(6);
t22 = t71 * mrSges(5,1) - mrSges(5,3) * t99;
t26 = mrSges(6,1) * t99 + t71 * mrSges(6,2);
t94 = -t22 + t26;
t24 = mrSges(6,1) * t102 - mrSges(6,3) * t71;
t25 = -mrSges(7,1) * t102 + t71 * mrSges(7,2);
t93 = -t24 + t25;
t27 = pkin(3) * t71 - pkin(8) * t73 + qJ(2);
t7 = t72 * t100 + t70 * t27;
t91 = t90 * pkin(8) ^ 2;
t65 = t71 ^ 2;
t67 = t73 ^ 2;
t89 = t67 + t65;
t88 = qJ(5) * t72;
t86 = t89 * mrSges(4,3);
t23 = mrSges(7,1) * t99 - t71 * mrSges(7,3);
t85 = -t70 * qJ(5) - pkin(3);
t40 = t70 * t100;
t6 = t72 * t27 - t40;
t84 = -t74 - t88;
t4 = -qJ(5) * t71 - t7;
t5 = -t71 * pkin(4) - t6;
t82 = -t4 * t72 + t5 * t70;
t81 = -t6 * t70 + t7 * t72;
t80 = t70 * mrSges(5,1) + t72 * mrSges(5,2);
t79 = -t70 * mrSges(6,2) - t72 * mrSges(6,3);
t78 = (Ifges(5,5) + Ifges(7,5)) * t99 + t120 * t102 + t119 * t71;
t76 = qJ(2) ^ 2;
t75 = qJ(5) ^ 2;
t68 = t74 ^ 2;
t59 = Ifges(5,5) * t70;
t58 = Ifges(7,5) * t70;
t57 = Ifges(5,6) * t72;
t53 = t67 * t74;
t52 = t67 * t68;
t50 = pkin(4) * t102;
t42 = t71 * t88;
t39 = t113 * t72;
t38 = t113 * t70;
t37 = Ifges(5,1) * t70 + t110;
t36 = Ifges(5,2) * t72 + t111;
t35 = -Ifges(6,2) * t70 - t106;
t34 = -Ifges(7,2) * t72 + t105;
t33 = -Ifges(6,3) * t72 - t107;
t32 = Ifges(7,3) * t70 - t104;
t31 = -mrSges(7,2) * t70 - mrSges(7,3) * t72;
t29 = mrSges(6,2) * t72 - mrSges(6,3) * t70;
t28 = -pkin(4) * t72 + t85;
t21 = -mrSges(5,2) * t71 - mrSges(5,3) * t102;
t19 = t80 * t73;
t18 = t79 * t73;
t17 = (-t72 * mrSges(7,2) + t70 * mrSges(7,3)) * t73;
t16 = t69 * t72 + t85;
t14 = t109 + (-Ifges(6,2) * t72 + t107) * t73;
t13 = Ifges(7,4) * t71 + (Ifges(7,2) * t70 + t104) * t73;
t12 = Ifges(6,5) * t71 + (Ifges(6,3) * t70 - t106) * t73;
t11 = Ifges(7,5) * t71 + (Ifges(7,3) * t72 + t105) * t73;
t10 = Ifges(5,5) * t71 + (Ifges(5,1) * t72 - t111) * t73;
t9 = t108 + (-Ifges(5,2) * t70 + t110) * t73;
t8 = t84 * t73 + t50;
t3 = t50 + (qJ(6) * t70 + t84) * t73;
t2 = -pkin(5) * t102 - t4;
t1 = t40 + (pkin(5) * t73 - t27) * t72 + t69 * t71;
t15 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t115) + 0.2e1 * t1 * t23 + 0.2e1 * t3 * t17 + 0.2e1 * t8 * t18 + 0.2e1 * t2 * t25 + 0.2e1 * t7 * t21 + 0.2e1 * t6 * t22 + 0.2e1 * t4 * t24 + 0.2e1 * t5 * t26 + Ifges(3,1) + Ifges(2,3) + t86 * t117 + (mrSges(4,1) * t115 + Ifges(4,2) * t71 + t78) * t71 + m(4) * (t65 * t68 + t52 + t76) + (m(3) * (pkin(1) ^ 2 + t76)) + m(5) * (t6 ^ 2 + t7 ^ 2 + t52) + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + ((mrSges(4,2) * t115) + Ifges(4,1) * t73 - 0.2e1 * Ifges(4,4) * t71 + t19 * t117 + (t10 + t11 - t14 - t109) * t72 + (t12 + t13 - t9 - t108) * t70) * t73; -(m(3) * pkin(1)) + mrSges(3,2) - t86 + (-t17 - t18 - t19) * t73 + ((t21 + t93) * t72 + (t23 + t94) * t70) * t71 + m(6) * (t82 * t71 - t73 * t8) + m(7) * (-t73 * t3 + (t1 * t70 + t2 * t72) * t71) + m(5) * (t81 * t71 + t53) + m(4) * (t65 * t74 + t53); m(3) + m(4) * t89 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t90 * t65 + t67); -pkin(3) * t19 + t16 * t17 + t38 * t23 + t39 * t25 + t8 * t29 + t3 * t31 + m(7) * (t1 * t38 + t16 * t3 + t2 * t39) + (-(t74 * mrSges(4,2)) + t59 / 0.2e1 + t57 / 0.2e1 - Ifges(4,6) + t58 / 0.2e1) * t71 + (t7 * mrSges(5,3) + t2 * mrSges(7,1) - t4 * mrSges(6,1) + t9 / 0.2e1 - t12 / 0.2e1 - t13 / 0.2e1 + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1) * t71) * t72 + (-t6 * mrSges(5,3) + t1 * mrSges(7,1) + t5 * mrSges(6,1) - t109 / 0.2e1 + t10 / 0.2e1 + t11 / 0.2e1 - t14 / 0.2e1) * t70 + ((t21 - t24) * t72 + t94 * t70 + m(5) * t81 + m(6) * t82) * pkin(8) + (Ifges(4,5) + t118 * t74 + (t32 / 0.2e1 - t35 / 0.2e1 + t37 / 0.2e1) * t72 + (t33 / 0.2e1 + t34 / 0.2e1 - t36 / 0.2e1) * t70) * t73 + (m(6) * t8 + t18) * t28; (-m(6) * t28 - m(7) * t16 + t118 - t29 - t31) * t73 + (m(7) * (t38 * t70 + t39 * t72) - mrSges(4,2) + (mrSges(5,3) + t97 + (m(5) + m(6)) * pkin(8)) * t90) * t71; -0.2e1 * pkin(3) * t30 + 0.2e1 * t16 * t31 + 0.2e1 * t28 * t29 + Ifges(4,3) + m(6) * (t28 ^ 2 + t91) + m(7) * (t16 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (pkin(3) ^ 2 + t91) + (t39 * t116 - t33 - t34 + t36) * t72 + (t38 * t116 + t32 - t35 + t37) * t70 + 0.2e1 * (mrSges(6,1) + mrSges(5,3)) * pkin(8) * t90; t93 * qJ(5) + m(6) * (-pkin(4) * t5 - qJ(5) * t4) + m(7) * (qJ(5) * t2 + t1 * t69) + (-Ifges(6,4) * t72 - Ifges(5,6) * t70) * t73 + t78 + t69 * t23 - pkin(4) * t26 - t4 * mrSges(6,3) + t5 * mrSges(6,2) + t6 * mrSges(5,1) - t7 * mrSges(5,2) - t1 * mrSges(7,3) + t2 * mrSges(7,2); (-mrSges(5,1) + t96) * t103 + m(6) * (-pkin(4) * t103 + t42) + m(7) * (t69 * t103 + t42) + (-mrSges(5,2) + t95) * t101; m(7) * (qJ(5) * t39 + t38 * t69) + t58 - t38 * mrSges(7,3) + t59 + t57 + t39 * mrSges(7,2) + (-pkin(4) * mrSges(6,1) + t69 * mrSges(7,1) - Ifges(6,4)) * t70 + (t97 * qJ(5) - t120) * t72 + (m(6) * (-pkin(4) * t70 + t88) - t79 - t80) * pkin(8); -0.2e1 * pkin(4) * mrSges(6,2) - 0.2e1 * t69 * mrSges(7,3) + 0.2e1 * t95 * qJ(5) + m(6) * (pkin(4) ^ 2 + t75) + m(7) * (t69 ^ 2 + t75) + t119; m(6) * t5 + m(7) * t1 + t23 + t26; t114 * t103; m(7) * t38 + (m(6) * pkin(8) + t97) * t70; -m(6) * pkin(4) + m(7) * t69 + t96; t114; m(7) * t2 + t25; m(7) * t101; m(7) * t39 + t72 * mrSges(7,1); m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
