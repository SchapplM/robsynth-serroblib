% Calculate time derivative of joint inertia matrix for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:20
% EndTime: 2019-12-31 19:07:23
% DurationCPUTime: 1.42s
% Computational Cost: add. (3263->199), mult. (7145->309), div. (0->0), fcn. (7260->8), ass. (0->99)
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t92 = cos(qJ(4));
t147 = (t88 ^ 2 + t91 ^ 2) * t92;
t116 = qJD(5) * t91;
t87 = cos(pkin(9));
t93 = cos(qJ(3));
t124 = t87 * t93;
t86 = sin(pkin(9));
t90 = sin(qJ(3));
t68 = -t90 * t86 + t124;
t69 = t86 * t93 + t90 * t87;
t89 = sin(qJ(4));
t100 = t92 * t68 - t69 * t89;
t59 = t68 * qJD(3);
t60 = t69 * qJD(3);
t35 = t100 * qJD(4) + t59 * t92 - t60 * t89;
t52 = t68 * t89 + t69 * t92;
t99 = t52 * t116 + t35 * t88;
t120 = pkin(6) + qJ(2);
t73 = t120 * t86;
t66 = t93 * t73;
t74 = t120 * t87;
t53 = -t74 * t90 - t66;
t49 = -pkin(7) * t69 + t53;
t54 = -t90 * t73 + t93 * t74;
t50 = pkin(7) * t68 + t54;
t101 = t92 * t49 - t50 * t89;
t46 = -qJD(3) * t66 + qJD(2) * t124 + (-qJD(2) * t86 - qJD(3) * t74) * t90;
t40 = -pkin(7) * t60 + t46;
t47 = -t69 * qJD(2) - t54 * qJD(3);
t96 = -t59 * pkin(7) + t47;
t11 = t101 * qJD(4) + t92 * t40 + t89 * t96;
t24 = t49 * t89 + t50 * t92;
t111 = -pkin(2) * t87 - pkin(1);
t55 = -pkin(3) * t68 + t111;
t25 = -pkin(4) * t100 - pkin(8) * t52 + t55;
t15 = -t24 * t88 + t25 * t91;
t140 = pkin(3) * t60;
t36 = t52 * qJD(4) + t59 * t89 + t92 * t60;
t17 = pkin(4) * t36 - pkin(8) * t35 + t140;
t2 = qJD(5) * t15 + t11 * t91 + t17 * t88;
t16 = t24 * t91 + t25 * t88;
t3 = -qJD(5) * t16 - t11 * t88 + t17 * t91;
t149 = t2 * t91 - t3 * t88;
t146 = -t15 * t88 + t16 * t91;
t145 = 2 * m(6);
t144 = -2 * mrSges(5,3);
t12 = qJD(4) * t24 + t40 * t89 - t92 * t96;
t143 = 0.2e1 * t12;
t142 = -0.2e1 * t101;
t137 = Ifges(6,4) * t88;
t136 = Ifges(6,4) * t91;
t135 = Ifges(6,6) * t88;
t134 = t12 * t101;
t129 = t35 * t91;
t128 = t52 * t88;
t127 = t52 * t91;
t123 = t89 * mrSges(5,1);
t75 = -mrSges(6,1) * t91 + mrSges(6,2) * t88;
t122 = t89 * t75;
t121 = t92 * mrSges(5,2);
t119 = Ifges(6,5) * t129 + Ifges(6,3) * t36;
t118 = pkin(3) * qJD(4);
t117 = qJD(5) * t88;
t115 = 0.2e1 * t140;
t114 = t52 * t117;
t98 = t114 - t129;
t8 = t99 * mrSges(6,1) - t98 * mrSges(6,2);
t112 = m(6) * t12 + t8;
t109 = t36 * mrSges(5,1) + t35 * mrSges(5,2);
t108 = -(2 * Ifges(5,4)) - t135;
t106 = mrSges(6,3) * t147;
t105 = mrSges(6,1) * t88 + mrSges(6,2) * t91;
t104 = Ifges(6,1) * t91 - t137;
t103 = -Ifges(6,2) * t88 + t136;
t102 = Ifges(6,5) * t88 + Ifges(6,6) * t91;
t71 = t103 * qJD(5);
t72 = t104 * qJD(5);
t76 = Ifges(6,2) * t91 + t137;
t77 = Ifges(6,1) * t88 + t136;
t97 = t77 * t116 - t76 * t117 + t91 * t71 + t88 * t72;
t13 = mrSges(6,1) * t36 + t98 * mrSges(6,3);
t14 = -mrSges(6,2) * t36 - t99 * mrSges(6,3);
t37 = mrSges(6,2) * t100 - mrSges(6,3) * t128;
t38 = -mrSges(6,1) * t100 - mrSges(6,3) * t127;
t95 = -t88 * t13 + m(6) * (-t15 * t116 - t16 * t117 + t149) + t91 * t14 - t38 * t116 - t37 * t117;
t20 = -Ifges(6,6) * t100 + t103 * t52;
t21 = -Ifges(6,5) * t100 + t104 * t52;
t6 = -t98 * Ifges(6,4) - t99 * Ifges(6,2) + Ifges(6,6) * t36;
t7 = -t98 * Ifges(6,1) - t99 * Ifges(6,4) + Ifges(6,5) * t36;
t70 = t105 * qJD(5);
t81 = Ifges(6,5) * t116;
t94 = t21 * t116 / 0.2e1 - t101 * t70 + t77 * t129 / 0.2e1 + Ifges(5,5) * t35 + t88 * t7 / 0.2e1 - t71 * t128 / 0.2e1 + t72 * t127 / 0.2e1 - t100 * (-Ifges(6,6) * t117 + t81) / 0.2e1 + t91 * t6 / 0.2e1 - t11 * mrSges(5,2) + (t102 / 0.2e1 - Ifges(5,6)) * t36 - t99 * t76 / 0.2e1 + (t75 - mrSges(5,1)) * t12 - (t52 * t77 + t20) * t117 / 0.2e1 + ((-t15 * t91 - t16 * t88) * qJD(5) + t149) * mrSges(6,3);
t80 = -pkin(3) * t92 - pkin(4);
t79 = pkin(3) * t89 + pkin(8);
t57 = t59 * mrSges(4,2);
t30 = t105 * t52;
t1 = [t24 * t36 * t144 + 0.2e1 * t111 * (t60 * mrSges(4,1) + t57) + 0.2e1 * t15 * t13 + 0.2e1 * t16 * t14 + t8 * t142 + t30 * t143 + 0.2e1 * t2 * t37 + 0.2e1 * t3 * t38 + 0.2e1 * t55 * t109 - 0.2e1 * t68 * Ifges(4,2) * t60 + 0.2e1 * t69 * t59 * Ifges(4,1) + (mrSges(5,3) * t142 - t20 * t88 + t21 * t91) * t35 + 0.2e1 * m(4) * (t46 * t54 + t47 * t53) + 0.2e1 * m(5) * (t11 * t24 + t55 * t140 - t134) + (t15 * t3 + t16 * t2 - t134) * t145 - (mrSges(5,1) * t115 + t11 * t144 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t36 + t108 * t35 + t119) * t100 + (mrSges(5,2) * t115 + mrSges(5,3) * t143 + 0.2e1 * Ifges(5,1) * t35 - t88 * t6 + t91 * t7 + (Ifges(6,5) * t91 + t108) * t36 + (t100 * t102 - t91 * t20 - t88 * t21) * qJD(5)) * t52 + 0.2e1 * (t59 * t68 - t60 * t69) * Ifges(4,4) + 0.2e1 * (t46 * t68 - t47 * t69 - t53 * t59 - t54 * t60) * mrSges(4,3) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t86 ^ 2 + t87 ^ 2) * qJD(2); m(6) * (t146 * qJD(5) + t2 * t88 + t3 * t91) + t37 * t116 + t88 * t14 - t38 * t117 + t91 * t13 + t57 - (-m(5) * pkin(3) - mrSges(4,1)) * t60 + t109; 0; t94 + t112 * t80 + t95 * t79 + (m(5) * (t11 * t89 - t12 * t92) + (-t92 * t35 - t89 * t36) * mrSges(5,3) + ((m(5) * t24 + m(6) * t146 + mrSges(5,3) * t100 + t91 * t37 - t88 * t38) * t92 + (t52 * mrSges(5,3) + t30 - (m(5) + m(6)) * t101) * t89) * qJD(4)) * pkin(3) - t46 * mrSges(4,2) + t47 * mrSges(4,1) + Ifges(4,5) * t59 - Ifges(4,6) * t60; 0; 0.2e1 * t80 * t70 + (-0.2e1 * t121 - 0.2e1 * t123 + (t147 * t79 + t80 * t89) * t145 + 0.2e1 * t122 + 0.2e1 * t106) * t118 + t97; -t112 * pkin(4) + t95 * pkin(8) + t94; 0; (t80 - pkin(4)) * t70 + (-t121 - t123 + m(6) * (-pkin(4) * t89 + t147 * pkin(8)) + t122 + t106) * t118 + t97; -0.2e1 * pkin(4) * t70 + t97; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t114 - t99 * Ifges(6,6) + t119; -t70; t81 - t105 * t92 * t118 + (t75 * t79 - t135) * qJD(5); t81 + (t75 * pkin(8) - t135) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
