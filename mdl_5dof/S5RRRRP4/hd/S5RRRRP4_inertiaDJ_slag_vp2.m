% Calculate time derivative of joint inertia matrix for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:41
% EndTime: 2019-12-31 21:50:44
% DurationCPUTime: 1.02s
% Computational Cost: add. (1610->185), mult. (3666->244), div. (0->0), fcn. (2793->6), ass. (0->95)
t149 = -mrSges(5,1) - mrSges(6,1);
t148 = -mrSges(5,2) + mrSges(6,3);
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t147 = t78 ^ 2 + t80 ^ 2;
t144 = (mrSges(5,3) + mrSges(6,2));
t146 = 2 * t144;
t109 = qJD(3) * t80;
t110 = qJD(3) * t78;
t143 = -mrSges(4,1) * t109 + mrSges(4,2) * t110;
t81 = cos(qJ(2));
t142 = t147 * t81;
t141 = qJD(3) + qJD(4);
t124 = cos(qJ(4));
t102 = t124 * t78;
t79 = sin(qJ(2));
t66 = pkin(1) * t79 + pkin(7);
t125 = -pkin(8) - t66;
t73 = t80 * pkin(8);
t52 = t66 * t80 + t73;
t77 = sin(qJ(4));
t27 = -t125 * t102 + t77 * t52;
t111 = pkin(1) * qJD(2);
t108 = t81 * t111;
t96 = t80 * t108;
t99 = qJD(3) * t125;
t43 = t78 * t99 + t96;
t97 = t78 * t108;
t89 = t80 * t99 - t97;
t5 = -t27 * qJD(4) + t124 * t43 + t77 * t89;
t116 = t77 * t78;
t28 = t125 * t116 + t124 * t52;
t6 = t28 * qJD(4) - t124 * t89 + t77 * t43;
t140 = t148 * t5 + t149 * t6;
t131 = -pkin(8) - pkin(7);
t100 = qJD(3) * t131;
t115 = t77 * t80;
t62 = pkin(7) * t80 + t73;
t92 = t131 * t124;
t41 = t77 * t62 - t78 * t92;
t57 = t78 * t100;
t23 = -t41 * qJD(4) + t100 * t115 + t124 * t57;
t42 = t131 * t116 + t124 * t62;
t24 = t42 * qJD(4) - t92 * t109 + t77 * t57;
t139 = t148 * t23 + t149 * t24;
t138 = (t147 * mrSges(4,3) - mrSges(3,2)) * t81 + (-mrSges(4,1) * t80 + mrSges(4,2) * t78 - mrSges(3,1)) * t79;
t137 = 2 * m(5);
t136 = 2 * m(6);
t91 = t124 * t80 - t116;
t36 = t141 * t91;
t54 = t102 + t115;
t37 = t141 * t54;
t15 = mrSges(6,1) * t37 - mrSges(6,3) * t36;
t135 = 0.2e1 * t15;
t16 = mrSges(5,1) * t37 + mrSges(5,2) * t36;
t134 = 0.2e1 * t16;
t38 = -mrSges(6,1) * t91 - mrSges(6,3) * t54;
t133 = 0.2e1 * t38;
t132 = m(5) * pkin(3);
t130 = pkin(1) * t81;
t39 = -mrSges(5,1) * t91 + mrSges(5,2) * t54;
t129 = pkin(3) * t39;
t128 = pkin(3) * t77;
t123 = Ifges(4,4) * t78;
t31 = t36 * mrSges(6,2);
t71 = pkin(3) * t110;
t72 = t79 * t111;
t58 = t72 + t71;
t119 = t58 * t39;
t107 = qJD(4) * t128;
t106 = t124 * pkin(3);
t69 = -t80 * pkin(3) - pkin(2);
t103 = t27 * t6 + t28 * t5;
t101 = t42 * t23 + t24 * t41;
t98 = (-Ifges(5,6) + Ifges(6,6)) * t37 + (Ifges(6,4) + Ifges(5,5)) * t36;
t95 = t54 * t107;
t94 = qJD(4) * t106;
t90 = t23 * t28 + t24 * t27 + t41 * t6 + t42 * t5;
t30 = -pkin(4) * t91 - t54 * qJ(5) + t69;
t10 = pkin(4) * t37 - qJ(5) * t36 - qJD(5) * t54 + t71;
t88 = (Ifges(4,1) * t80 - t123) * t110 + (0.2e1 * Ifges(4,4) * t80 + (Ifges(4,1) - Ifges(4,2)) * t78) * t109 + 0.2e1 * (-Ifges(5,2) - Ifges(6,3)) * t91 * t37 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t54 * t36 + 0.2e1 * (-Ifges(6,5) + Ifges(5,4)) * (t36 * t91 - t54 * t37);
t85 = (-pkin(4) * t36 - qJ(5) * t37 + qJD(5) * t91) * mrSges(6,2) + t98;
t63 = t94 + qJD(5);
t83 = -mrSges(5,2) * t94 + t63 * mrSges(6,3) + t149 * t107;
t64 = qJ(5) + t128;
t67 = -t106 - pkin(4);
t82 = Ifges(4,5) * t109 - Ifges(4,6) * t110 + t67 * t31 + t98 + (-t37 * t64 + t63 * t91 + t95) * mrSges(6,2) + (-t36 * t106 - t37 * t128 + t91 * t94 + t95) * mrSges(5,3);
t76 = qJD(5) * mrSges(6,3);
t68 = -pkin(2) - t130;
t61 = Ifges(4,2) * t80 + t123;
t59 = t69 - t130;
t56 = (mrSges(4,1) * t78 + mrSges(4,2) * t80) * qJD(3);
t25 = t30 - t130;
t8 = t10 + t72;
t1 = [0.2e1 * t119 + t59 * t134 + 0.2e1 * t68 * t56 + t8 * t133 + t25 * t135 - t61 * t110 + 0.2e1 * (m(4) * (t142 * t66 + t68 * t79) + t138) * t111 + t88 + (t58 * t59 + t103) * t137 + (t25 * t8 + t103) * t136 + (t27 * t36 - t28 * t37 + t5 * t91 + t6 * t54) * t146; t119 + m(5) * (t58 * t69 + t59 * t71 + t90) + m(6) * (t10 * t25 + t30 * t8 + t90) + (t68 - pkin(2)) * t56 + (t10 + t8) * t38 + (m(4) * (-pkin(2) * t79 + t142 * pkin(7)) + t138) * t111 + (t69 + t59) * t16 + (t25 + t30) * t15 + (-t61 + t129) * t110 + t88 + t144 * ((t24 + t6) * t54 - (-t23 - t5) * t91 + (-t28 - t42) * t37 + (t27 + t41) * t36); (-t61 + 0.2e1 * t129) * t110 + t69 * t134 - 0.2e1 * pkin(2) * t56 + t10 * t133 + t30 * t135 + t88 + (t69 * t71 + t101) * t137 + (t10 * t30 + t101) * t136 + (t23 * t91 + t24 * t54 + t36 * t41 - t37 * t42) * t146; m(6) * (t27 * t107 + t28 * t63 + t5 * t64 + t6 * t67) + t82 - mrSges(4,2) * t96 - mrSges(4,1) * t97 + (-t124 * t6 + t5 * t77 + (t124 * t28 + t27 * t77) * qJD(4)) * t132 + t143 * t66 + t140; m(6) * (t41 * t107 + t23 * t64 + t24 * t67 + t42 * t63) + t82 + (-t124 * t24 + t23 * t77 + (t124 * t42 + t41 * t77) * qJD(4)) * t132 + t143 * pkin(7) + t139; 0.2e1 * m(6) * (t67 * t107 + t63 * t64) + 0.2e1 * t83; m(6) * (-pkin(4) * t6 + qJ(5) * t5 + qJD(5) * t28) + t85 + t140; m(6) * (-pkin(4) * t24 + qJ(5) * t23 + qJD(5) * t42) + t85 + t139; t76 + m(6) * (-pkin(4) * t107 + qJ(5) * t63 + qJD(5) * t64) + t83; 0.2e1 * m(6) * qJ(5) * qJD(5) + 0.2e1 * t76; m(6) * t6 + t31; m(6) * t24 + t31; m(6) * t107; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
