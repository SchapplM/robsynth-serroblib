% Calculate time derivative of joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:27
% EndTime: 2019-12-31 19:46:31
% DurationCPUTime: 1.24s
% Computational Cost: add. (1216->228), mult. (2706->349), div. (0->0), fcn. (2118->6), ass. (0->104)
t85 = sin(pkin(8));
t86 = cos(pkin(8));
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t95 = t90 * t85 + t88 * t86;
t53 = t95 * qJD(5);
t130 = -t85 * t88 + t86 * t90;
t84 = t86 ^ 2;
t102 = (t85 ^ 2 + t84) * qJD(4);
t129 = 2 * m(5);
t128 = 2 * m(6);
t127 = -2 * pkin(1);
t89 = sin(qJ(2));
t108 = t89 * qJ(3);
t91 = cos(qJ(2));
t98 = -pkin(2) * t91 - t108;
t70 = -pkin(1) + t98;
t125 = -0.2e1 * t70;
t124 = -t95 / 0.2e1;
t123 = t130 / 0.2e1;
t122 = t86 / 0.2e1;
t121 = pkin(3) + pkin(6);
t81 = t91 * pkin(6);
t87 = -pkin(2) - qJ(4);
t120 = -pkin(7) + t87;
t119 = mrSges(5,2) * t85;
t118 = Ifges(5,4) * t85;
t117 = Ifges(5,4) * t86;
t116 = t85 * Ifges(5,1);
t114 = t85 * t89;
t113 = t85 * t91;
t111 = t86 * t91;
t107 = qJD(2) * t89;
t101 = pkin(2) * t107 - t89 * qJD(3);
t37 = -qJD(4) * t91 + (-qJ(3) * t91 + qJ(4) * t89) * qJD(2) + t101;
t106 = qJD(2) * t91;
t67 = t121 * t106;
t15 = t86 * t37 + t85 * t67;
t59 = t87 * t91 - pkin(1) - t108;
t72 = t121 * t89;
t30 = t86 * t59 + t85 * t72;
t52 = t130 * qJD(5);
t110 = -Ifges(6,5) * t53 - Ifges(6,6) * t52;
t73 = t91 * pkin(3) + t81;
t105 = 2 * mrSges(6,3);
t43 = t130 * t91;
t22 = -qJD(5) * t43 + t95 * t107;
t23 = t107 * t130 + t91 * t53;
t104 = Ifges(6,5) * t22 + Ifges(6,6) * t23 + Ifges(6,3) * t106;
t103 = t86 * t107;
t7 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t26 = mrSges(6,1) * t52 - t53 * mrSges(6,2);
t14 = -t85 * t37 + t86 * t67;
t100 = t91 * mrSges(4,2) - t89 * mrSges(4,3);
t99 = -Ifges(5,5) * t85 - Ifges(5,6) * t86;
t97 = t14 * t86 + t15 * t85;
t63 = t86 * t72;
t21 = t89 * pkin(4) + t63 + (pkin(7) * t91 - t59) * t85;
t25 = -pkin(7) * t111 + t30;
t5 = t21 * t90 - t25 * t88;
t6 = t21 * t88 + t25 * t90;
t96 = t130 * t53 - t52 * t95;
t68 = t120 * t85;
t69 = t120 * t86;
t33 = t68 * t90 + t69 * t88;
t32 = -t68 * t88 + t69 * t90;
t46 = -mrSges(5,1) * t103 + t107 * t119;
t8 = (pkin(4) * t91 - pkin(7) * t114) * qJD(2) + t14;
t9 = pkin(7) * t103 + t15;
t1 = qJD(5) * t5 + t8 * t88 + t9 * t90;
t2 = -qJD(5) * t6 + t8 * t90 - t88 * t9;
t93 = -t1 * t95 - t130 * t2 + t5 * t53 - t6 * t52;
t12 = -qJD(4) * t95 + qJD(5) * t32;
t13 = -qJD(4) * t130 - qJD(5) * t33;
t92 = -t12 * t95 - t13 * t130 + t32 * t53 - t33 * t52;
t76 = pkin(4) * t85 + qJ(3);
t71 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t66 = t121 * t107;
t65 = -t89 * mrSges(5,2) - mrSges(5,3) * t111;
t64 = t89 * mrSges(5,1) + mrSges(5,3) * t113;
t58 = (mrSges(5,3) * t86 * t89 - mrSges(5,2) * t91) * qJD(2);
t57 = (mrSges(5,1) * t91 - mrSges(5,3) * t114) * qJD(2);
t56 = (mrSges(5,1) * t86 - t119) * t91;
t51 = pkin(4) * t111 + t73;
t44 = t95 * t91;
t42 = (-pkin(4) * t86 - t121) * t107;
t41 = (t91 * Ifges(5,5) + (t116 + t117) * t89) * qJD(2);
t40 = (t91 * Ifges(5,6) + (t86 * Ifges(5,2) + t118) * t89) * qJD(2);
t39 = mrSges(6,1) * t89 + mrSges(6,3) * t44;
t38 = -mrSges(6,2) * t89 - mrSges(6,3) * t43;
t36 = Ifges(6,1) * t130 - Ifges(6,4) * t95;
t35 = Ifges(6,4) * t130 - Ifges(6,2) * t95;
t34 = mrSges(6,1) * t95 + mrSges(6,2) * t130;
t29 = -t85 * t59 + t63;
t28 = -Ifges(6,1) * t53 - Ifges(6,4) * t52;
t27 = -Ifges(6,4) * t53 - Ifges(6,2) * t52;
t24 = mrSges(6,1) * t43 - mrSges(6,2) * t44;
t17 = -Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t89;
t16 = -Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t89;
t11 = -mrSges(6,2) * t106 + t23 * mrSges(6,3);
t10 = mrSges(6,1) * t106 - t22 * mrSges(6,3);
t4 = Ifges(6,1) * t22 + Ifges(6,4) * t23 + Ifges(6,5) * t106;
t3 = Ifges(6,4) * t22 + Ifges(6,2) * t23 + Ifges(6,6) * t106;
t18 = [-t40 * t111 - t41 * t113 + t89 * t104 + 0.2e1 * t14 * t64 + 0.2e1 * t15 * t65 - 0.2e1 * t66 * t56 + 0.2e1 * t73 * t46 + 0.2e1 * t29 * t57 + 0.2e1 * t30 * t58 - t43 * t3 - t44 * t4 + 0.2e1 * t51 * t7 + 0.2e1 * t1 * t38 + 0.2e1 * t2 * t39 + 0.2e1 * t42 * t24 + t22 * t17 + t23 * t16 + 0.2e1 * t5 * t10 + 0.2e1 * t6 * t11 + (t14 * t29 + t15 * t30 - t66 * t73) * t129 + (t1 * t6 + t2 * t5 + t42 * t51) * t128 + (((mrSges(3,2) * t127) + mrSges(4,3) * t125 - Ifges(6,5) * t44 - Ifges(6,6) * t43 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t99) * t91) * t91 + (mrSges(4,2) * t125 + (mrSges(3,1) * t127) + 0.2e1 * (-Ifges(3,4) - Ifges(4,6) - t99) * t89 + (-t84 * Ifges(5,2) + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + (2 * Ifges(5,3)) + Ifges(6,3) + (-t116 - 0.2e1 * t117) * t85) * t91) * t89) * qJD(2) + 0.2e1 * (m(4) * t70 + t100) * (-qJ(3) * t106 + t101); m(5) * (-qJ(3) * t66 + t97 * t87 + (-t29 * t86 - t30 * t85) * qJD(4)) + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t122 - Ifges(5,6) * t85 / 0.2e1 + Ifges(6,5) * t123 + Ifges(6,6) * t124 - Ifges(4,4) + Ifges(3,5)) * t91 + (t85 * (Ifges(5,1) * t86 - t118) / 0.2e1 + (-Ifges(5,2) * t85 + t117) * t122 - qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6)) * t89 + (m(4) * t98 - t91 * mrSges(3,1) + t89 * mrSges(3,2) + t100) * pkin(6)) * qJD(2) + t89 * t110 / 0.2e1 + t76 * t7 - t66 * t71 - t53 * t17 / 0.2e1 - t43 * t27 / 0.2e1 - t44 * t28 / 0.2e1 + qJ(3) * t46 + t51 * t26 - t52 * t16 / 0.2e1 + t32 * t10 + t33 * t11 + t23 * t35 / 0.2e1 + t22 * t36 / 0.2e1 + t12 * t38 + t13 * t39 + t42 * t34 + t93 * mrSges(6,3) + t4 * t123 + t3 * t124 + (t87 * t57 - qJD(4) * t64 - t14 * mrSges(5,3) + t41 / 0.2e1) * t86 + (t87 * t58 - qJD(4) * t65 - t15 * mrSges(5,3) - t40 / 0.2e1) * t85 + m(6) * (t1 * t33 + t12 * t6 + t13 * t5 + t2 * t32 + t42 * t76) + (m(4) * t81 + m(5) * t73 + m(6) * t51 + t91 * mrSges(4,1) + t24 + t56) * qJD(3); 0.2e1 * t76 * t26 - t95 * t27 + t130 * t28 - t52 * t35 - t53 * t36 + (qJD(3) * t76 + t12 * t33 + t13 * t32) * t128 + (qJ(3) * qJD(3) - t102 * t87) * t129 + t92 * t105 + 0.2e1 * mrSges(5,3) * t102 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3) + t34 + t71) * qJD(3); t130 * t10 + t95 * t11 + t52 * t38 - t53 * t39 + t86 * t57 + t85 * t58 + (m(4) * pkin(6) + mrSges(4,1)) * t106 - m(6) * t93 + m(5) * t97; -m(5) * t102 - m(6) * t92 + t96 * t105; -0.2e1 * m(6) * t96; -m(5) * t66 + m(6) * t42 + t46 + t7; (m(5) + m(6)) * qJD(3) + t26; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t104; mrSges(6,1) * t13 - mrSges(6,2) * t12 + t110; -mrSges(6,1) * t53 - mrSges(6,2) * t52; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;
