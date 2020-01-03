% Calculate time derivative of joint inertia matrix for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:24
% EndTime: 2019-12-31 20:12:28
% DurationCPUTime: 1.50s
% Computational Cost: add. (892->242), mult. (2005->348), div. (0->0), fcn. (1289->4), ass. (0->110)
t122 = pkin(3) + pkin(6);
t129 = Ifges(6,4) + Ifges(5,5);
t128 = Ifges(6,2) + Ifges(5,3);
t72 = cos(qJ(2));
t103 = qJD(4) * t72;
t69 = sin(qJ(4));
t70 = sin(qJ(2));
t107 = qJD(2) * t70;
t71 = cos(qJ(4));
t99 = t71 * t107;
t77 = t69 * t103 + t99;
t127 = 2 * qJ(3);
t73 = -pkin(2) - pkin(7);
t94 = -qJ(3) * t70 - pkin(1);
t33 = t73 * t72 + t94;
t56 = t122 * t70;
t126 = t71 * t33 + t69 * t56;
t93 = pkin(2) * t107 - qJD(3) * t70;
t22 = (pkin(7) * t70 - qJ(3) * t72) * qJD(2) + t93;
t106 = qJD(2) * t72;
t47 = t122 * t106;
t4 = -qJD(4) * t126 - t22 * t69 + t47 * t71;
t125 = -0.2e1 * pkin(1);
t29 = -qJ(3) * t106 + t93;
t124 = 0.2e1 * t29;
t49 = -pkin(2) * t72 + t94;
t123 = -0.2e1 * t49;
t121 = Ifges(5,4) * t69;
t120 = Ifges(5,4) * t71;
t119 = Ifges(6,5) * t69;
t118 = Ifges(6,5) * t71;
t117 = Ifges(6,6) * t70;
t116 = t69 * t72;
t115 = t69 * t73;
t114 = t70 * Ifges(5,6);
t113 = t71 * t72;
t112 = t71 * t73;
t17 = -mrSges(5,2) * t106 + t77 * mrSges(5,3);
t20 = t77 * mrSges(6,2) + mrSges(6,3) * t106;
t110 = t17 + t20;
t100 = t69 * t107;
t76 = -t71 * t103 + t100;
t18 = mrSges(5,1) * t106 - t76 * mrSges(5,3);
t104 = qJD(4) * t71;
t19 = mrSges(6,2) * t100 + (-mrSges(6,1) * qJD(2) - mrSges(6,2) * t104) * t72;
t109 = t18 - t19;
t87 = Ifges(6,1) * t69 - t118;
t25 = t70 * Ifges(6,4) - t87 * t72;
t88 = Ifges(5,1) * t69 + t120;
t26 = t70 * Ifges(5,5) - t88 * t72;
t108 = t25 + t26;
t57 = t122 * t72;
t105 = qJD(4) * t69;
t102 = qJD(5) * t69;
t98 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t52 = Ifges(6,3) * t69 + t118;
t53 = -Ifges(5,2) * t69 + t120;
t97 = -t52 / 0.2e1 + t53 / 0.2e1;
t54 = Ifges(6,1) * t71 + t119;
t55 = Ifges(5,1) * t71 - t121;
t96 = t55 / 0.2e1 + t54 / 0.2e1;
t95 = m(4) * pkin(6) + mrSges(4,1);
t85 = -Ifges(6,3) * t71 + t119;
t23 = -t85 * t72 + t117;
t86 = Ifges(5,2) * t71 + t121;
t24 = -t86 * t72 + t114;
t92 = -t23 + t24 - t117;
t91 = (m(6) * t73 - mrSges(6,2)) * t69;
t90 = mrSges(5,1) * t71 - mrSges(5,2) * t69;
t51 = t69 * mrSges(5,1) + t71 * mrSges(5,2);
t89 = mrSges(6,1) * t71 + mrSges(6,3) * t69;
t50 = t69 * mrSges(6,1) - t71 * mrSges(6,3);
t84 = pkin(4) * t71 + qJ(5) * t69;
t83 = -pkin(4) * t69 + qJ(5) * t71;
t10 = qJ(5) * t70 + t126;
t14 = -t33 * t69 + t56 * t71;
t11 = -pkin(4) * t70 - t14;
t82 = t10 * t71 + t11 * t69;
t81 = t126 * t71 - t14 * t69;
t78 = t77 * Ifges(5,6) + t129 * t100 + t128 * t106;
t3 = t56 * t104 - t33 * t105 + t71 * t22 + t69 * t47;
t42 = mrSges(5,1) * t70 + mrSges(5,3) * t116;
t43 = -mrSges(6,1) * t70 - mrSges(6,2) * t116;
t44 = -mrSges(5,2) * t70 - mrSges(5,3) * t113;
t45 = -mrSges(6,2) * t113 + mrSges(6,3) * t70;
t75 = (t44 + t45) * t71 + (-t42 + t43) * t69;
t74 = m(6) * t83 - t50 - t51;
t64 = Ifges(6,6) * t104;
t48 = qJ(3) - t83;
t46 = t122 * t107;
t41 = t88 * qJD(4);
t40 = t87 * qJD(4);
t39 = t86 * qJD(4);
t38 = t85 * qJD(4);
t37 = t90 * qJD(4);
t36 = t89 * qJD(4);
t32 = t90 * t72;
t31 = t89 * t72;
t27 = t84 * qJD(4) - qJD(5) * t71 + qJD(3);
t21 = t84 * t72 + t57;
t13 = -t77 * mrSges(5,1) + t76 * mrSges(5,2);
t12 = -t77 * mrSges(6,1) - t76 * mrSges(6,3);
t9 = -t55 * t103 + (t72 * Ifges(5,5) + t88 * t70) * qJD(2);
t8 = -t54 * t103 + (t72 * Ifges(6,4) + t87 * t70) * qJD(2);
t7 = -t53 * t103 + (t72 * Ifges(5,6) + t86 * t70) * qJD(2);
t6 = -t52 * t103 + (t72 * Ifges(6,6) + t85 * t70) * qJD(2);
t5 = (t83 * qJD(4) + t102) * t72 + (-t84 - t122) * t107;
t2 = -pkin(4) * t106 - t4;
t1 = qJ(5) * t106 + qJD(5) * t70 + t3;
t15 = [m(4) * t49 * t124 + 0.2e1 * t1 * t45 + 0.2e1 * t10 * t20 + 0.2e1 * t11 * t19 + 0.2e1 * t21 * t12 + 0.2e1 * t57 * t13 + 0.2e1 * t14 * t18 + 0.2e1 * t126 * t17 + 0.2e1 * t2 * t43 + 0.2e1 * t3 * t44 + 0.2e1 * t5 * t31 - 0.2e1 * t46 * t32 + 0.2e1 * t4 * t42 + 0.2e1 * m(5) * (t126 * t3 + t14 * t4 - t46 * t57) + 0.2e1 * m(6) * (t1 * t10 + t11 * t2 + t21 * t5) + (-0.2e1 * t29 * mrSges(4,3) + (mrSges(3,1) * t125 + mrSges(4,2) * t123 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6)) * t70 + t92 * t71 + t108 * t69) * qJD(2) + t78) * t70 + (mrSges(4,2) * t124 + (t6 - t7) * t71 + (-t8 - t9) * t69 + (t92 * t69 + (-t129 * t70 - t108) * t71) * qJD(4) + (mrSges(3,2) * t125 + mrSges(4,3) * t123 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + (-Ifges(5,6) + Ifges(6,6)) * t71 - t129 * t69) * t72 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + t128) * t70) * qJD(2)) * t72; t70 * t64 / 0.2e1 + t57 * t37 + t21 * t36 + t48 * t12 + t5 * t50 - t46 * t51 + t27 * t31 + qJD(3) * t32 + qJ(3) * t13 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t4 * mrSges(5,3) + t109 * t73) * t71 + (-t1 * mrSges(6,2) - t3 * mrSges(5,3) + t6 / 0.2e1 - t7 / 0.2e1 + t110 * t73) * t69 + m(6) * (t1 * t115 - t2 * t112 + t27 * t21 + t48 * t5) + m(5) * (-qJ(3) * t46 + qJD(3) * t57 + t4 * t112 + t3 * t115) + ((-t38 / 0.2e1 + t39 / 0.2e1) * t71 + (t40 / 0.2e1 + t41 / 0.2e1) * t69 + t95 * qJD(3)) * t72 + ((-t10 * mrSges(6,2) - t126 * mrSges(5,3) - t114 / 0.2e1 + t23 / 0.2e1 - t24 / 0.2e1 - t96 * t72) * t71 + (-t25 / 0.2e1 - t26 / 0.2e1 - t11 * mrSges(6,2) + t14 * mrSges(5,3) + t97 * t72 - t98 * t70) * t69 + (m(5) * t81 + m(6) * t82 + t75) * t73) * qJD(4) + ((-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) + t98 * t71 + (-Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t69 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(6)) * t72 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t97 * t71 + t96 * t69 + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(6)) * t70) * qJD(2); t37 * t127 + 0.2e1 * t36 * t48 + (-t40 - t41) * t71 + (-t38 + t39) * t69 + ((t52 - t53) * t71 + (-t54 - t55) * t69) * qJD(4) + 0.2e1 * (m(6) * t48 + t50) * t27 + (0.2e1 * mrSges(4,3) + 0.2e1 * t51 + (m(4) + m(5)) * t127) * qJD(3); t109 * t71 + t110 * t69 + t95 * t106 + t75 * qJD(4) + m(6) * (t82 * qJD(4) + t1 * t69 - t2 * t71) + m(5) * (t81 * qJD(4) + t3 * t69 + t4 * t71); 0; 0; -Ifges(6,6) * t99 + qJD(5) * t45 + qJ(5) * t20 + m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t10) + t1 * mrSges(6,3) - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t19 + (-Ifges(6,6) * t69 - t129 * t71) * t103 + t78; t64 + qJD(5) * t91 + ((-qJ(5) * mrSges(6,2) - Ifges(5,6)) * t71 + (pkin(4) * mrSges(6,2) - t129) * t69 + t74 * t73) * qJD(4); m(6) * t102 + t74 * qJD(4); 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t2 + t19; qJD(4) * t91; m(6) * t105; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
