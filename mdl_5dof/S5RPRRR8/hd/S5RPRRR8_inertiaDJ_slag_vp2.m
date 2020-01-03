% Calculate time derivative of joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:43
% EndTime: 2019-12-31 19:05:46
% DurationCPUTime: 1.10s
% Computational Cost: add. (1415->180), mult. (2868->269), div. (0->0), fcn. (2199->6), ass. (0->91)
t74 = cos(qJ(5));
t75 = cos(qJ(4));
t102 = t74 * t75;
t71 = sin(qJ(5));
t72 = sin(qJ(4));
t104 = t71 * t72;
t43 = -t102 + t104;
t131 = 0.2e1 * t43;
t44 = t71 * t75 + t72 * t74;
t130 = -0.2e1 * t44;
t85 = mrSges(5,1) * t72 + mrSges(5,2) * t75;
t45 = t85 * qJD(4);
t115 = qJD(4) + qJD(5);
t95 = qJD(4) * t75;
t24 = -qJD(5) * t102 + t104 * t115 - t74 * t95;
t25 = t115 * t44;
t5 = mrSges(6,1) * t25 - mrSges(6,2) * t24;
t129 = t5 + t45;
t73 = sin(qJ(3));
t128 = qJD(3) * t73;
t96 = qJD(4) * t72;
t127 = t72 * (Ifges(5,1) * t95 - Ifges(5,4) * t96) + t24 * Ifges(6,1) * t130 + t25 * Ifges(6,2) * t131 + (t130 * t25 + t131 * t24) * Ifges(6,4);
t109 = t75 * pkin(4);
t76 = cos(qJ(3));
t77 = -pkin(1) - pkin(2);
t52 = -qJ(2) * t73 + t76 * t77;
t48 = pkin(3) - t52;
t41 = t48 + t109;
t126 = t41 * t5;
t62 = -pkin(3) - t109;
t125 = t62 * t5;
t26 = mrSges(6,1) * t43 + mrSges(6,2) * t44;
t53 = qJ(2) * t76 + t73 * t77;
t38 = t73 * qJD(2) + qJD(3) * t53;
t93 = pkin(4) * t96;
t34 = t38 - t93;
t124 = t34 * t26;
t98 = t72 ^ 2 + t75 ^ 2;
t91 = t98 * mrSges(5,3);
t123 = mrSges(4,2) - t91;
t54 = -mrSges(5,1) * t75 + mrSges(5,2) * t72;
t99 = t54 - mrSges(4,1);
t122 = -t26 - t99;
t121 = (-t74 * t24 + t71 * t25 + (t43 * t74 - t44 * t71) * qJD(5)) * mrSges(6,3);
t46 = Ifges(5,4) * t95 - Ifges(5,2) * t96;
t119 = -t75 * t46 - t127;
t40 = t43 * t73;
t39 = t44 * t73;
t97 = qJD(3) * t76;
t11 = -t115 * t39 - t43 * t97;
t12 = t115 * t40 - t44 * t97;
t118 = (t11 * t43 + t12 * t44 + t24 * t39 - t25 * t40) * mrSges(6,3);
t117 = t99 * t38;
t116 = -Ifges(6,5) * t24 - Ifges(6,6) * t25;
t113 = 2 * m(6);
t112 = -0.2e1 * t45;
t111 = -pkin(8) - pkin(7);
t49 = -pkin(7) + t53;
t108 = pkin(8) - t49;
t105 = t38 * t76;
t56 = Ifges(5,1) * t72 + Ifges(5,4) * t75;
t101 = t75 * t56;
t94 = 0.2e1 * mrSges(6,3);
t92 = qJD(4) * t111;
t37 = t76 * qJD(2) + qJD(3) * t52;
t90 = t98 * t37;
t89 = t98 * t49;
t88 = mrSges(6,1) * t12 - t11 * mrSges(6,2);
t87 = qJD(4) * t108;
t84 = -Ifges(5,5) * t75 + Ifges(5,6) * t72;
t35 = t108 * t72;
t36 = t108 * t75;
t16 = t35 * t74 + t36 * t71;
t17 = t35 * t71 - t36 * t74;
t57 = t111 * t72;
t58 = t111 * t75;
t32 = t57 * t74 + t58 * t71;
t33 = t57 * t71 - t58 * t74;
t18 = t37 * t75 + t72 * t87;
t19 = -t37 * t72 + t75 * t87;
t2 = qJD(5) * t16 + t18 * t74 + t19 * t71;
t3 = -qJD(5) * t17 - t18 * t71 + t19 * t74;
t80 = mrSges(6,1) * t3 - t2 * mrSges(6,2) - t116;
t50 = t72 * t92;
t51 = t75 * t92;
t14 = qJD(5) * t32 + t50 * t74 + t51 * t71;
t15 = -qJD(5) * t33 - t50 * t71 + t51 * t74;
t79 = mrSges(6,1) * t15 - t14 * mrSges(6,2) + t116;
t55 = Ifges(5,4) * t72 + Ifges(5,2) * t75;
t42 = (-mrSges(6,1) * t71 - mrSges(6,2) * t74) * qJD(5) * pkin(4);
t1 = [-0.2e1 * t124 - 0.2e1 * t126 + t48 * t112 + (-t72 * t55 + t101) * qJD(4) + 0.2e1 * t37 * mrSges(4,2) + 0.2e1 * m(5) * (t37 * t89 + t38 * t48) + 0.2e1 * m(4) * (t37 * t53 - t38 * t52) + (t16 * t3 + t17 * t2 + t34 * t41) * t113 + (-t16 * t24 + t17 * t25 + t2 * t43 + t3 * t44) * t94 - 0.2e1 * t117 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) - 0.2e1 * t91 * t37 - t119; t129 * t76 + (t122 * t73 + t123 * t76) * qJD(3) + m(6) * (t11 * t17 + t12 * t16 + t128 * t41 - t2 * t40 - t3 * t39 - t34 * t76) + m(5) * (-t105 + t73 * t90 + (t48 * t73 + t76 * t89) * qJD(3)) + m(4) * (t37 * t73 - t105 + (-t52 * t73 + t53 * t76) * qJD(3)) + t118; (-t11 * t40 - t12 * t39) * t113 + 0.4e1 * (m(5) * (-0.1e1 + t98) / 0.2e1 - m(6) / 0.2e1) * t73 * t97; t124 - t125 + t126 + (t48 + pkin(3)) * t45 + t117 - t123 * t37 + (-t101 + (-pkin(4) * t26 + t55) * t72) * qJD(4) + m(6) * (t14 * t17 + t15 * t16 + t2 * t33 + t3 * t32 + t34 * t62 + t41 * t93) + m(5) * (-pkin(3) * t38 + pkin(7) * t90) + ((t15 - t3) * t44 + (t14 - t2) * t43 + (-t17 + t33) * t25 + (t16 - t32) * t24) * mrSges(6,3) + t119; m(6) * (t11 * t33 + t12 * t32 - t14 * t40 - t15 * t39) - t118 + (-m(5) * pkin(3) + m(6) * t62 - t122) * t128 + (-m(6) * t93 + (m(5) * pkin(7) * t98 - t123) * qJD(3) - t129) * t76; -t55 * t96 + (t14 * t33 + t15 * t32 + t62 * t93) * t113 + 0.2e1 * t26 * t93 + 0.2e1 * t125 + pkin(3) * t112 + (qJD(4) * t56 + t46) * t75 + (-t14 * t43 - t15 * t44 + t32 * t24 - t33 * t25) * t94 + t127; -t85 * t37 + (t49 * t54 + t84) * qJD(4) + (m(6) * (t2 * t71 + t3 * t74 + (-t16 * t71 + t17 * t74) * qJD(5)) + t121) * pkin(4) + t80; (t73 * t96 - t75 * t97) * mrSges(5,2) + (-t72 * t97 - t73 * t95) * mrSges(5,1) + m(6) * (t11 * t71 + t12 * t74 + (t39 * t71 - t40 * t74) * qJD(5)) * pkin(4) + t88; (pkin(7) * t54 - t84) * qJD(4) + (m(6) * (t14 * t71 + t15 * t74 + (-t32 * t71 + t33 * t74) * qJD(5)) - t121) * pkin(4) + t79; 0.2e1 * t42; t80; t88; t79; t42; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
