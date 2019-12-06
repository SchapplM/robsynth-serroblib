% Calculate time derivative of joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:52
% EndTime: 2019-12-05 17:02:55
% DurationCPUTime: 0.96s
% Computational Cost: add. (845->168), mult. (2453->279), div. (0->0), fcn. (2226->8), ass. (0->89)
t132 = 2 * pkin(2);
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t131 = -mrSges(6,1) * t67 + mrSges(6,2) * t63 - mrSges(5,1);
t124 = qJD(3) + qJD(4);
t64 = sin(qJ(4));
t65 = sin(qJ(3));
t108 = t64 * t65;
t68 = cos(qJ(4));
t69 = cos(qJ(3));
t46 = -t68 * t69 + t108;
t31 = t124 * t46;
t47 = t64 * t69 + t65 * t68;
t95 = qJD(5) * t67;
t76 = -t31 * t63 + t47 * t95;
t66 = sin(qJ(2));
t104 = qJD(2) * t66;
t70 = cos(qJ(2));
t103 = qJD(2) * t70;
t98 = qJD(4) * t47;
t32 = qJD(3) * t47 + t98;
t14 = -t46 * t103 - t32 * t66;
t42 = t46 * t66;
t38 = t42 * t63 - t67 * t70;
t7 = qJD(5) * t38 + t104 * t63 + t14 * t67;
t39 = -t42 * t67 - t63 * t70;
t8 = -qJD(5) * t39 + t104 * t67 - t14 * t63;
t130 = -t8 * t63 + t67 * t7;
t106 = t63 ^ 2 + t67 ^ 2;
t129 = t106 * mrSges(6,3) - mrSges(5,2);
t101 = qJD(3) * t66;
t128 = -t65 * t101 + t69 * t103;
t127 = t103 * t66;
t126 = -t38 * t63 + t39 * t67;
t123 = m(5) / 0.2e1;
t122 = m(6) / 0.2e1;
t118 = Ifges(6,4) * t63;
t117 = Ifges(6,4) * t67;
t92 = t65 * t103;
t15 = -qJD(4) * t66 * t108 + (t124 * t69 * t66 + t92) * t68 + t128 * t64;
t41 = t47 * t66;
t116 = t15 * t41;
t114 = t31 * t67;
t111 = t47 * t63;
t110 = t47 * t67;
t83 = mrSges(6,1) * t63 + mrSges(6,2) * t67;
t48 = t83 * qJD(5);
t109 = t48 * t68;
t107 = -Ifges(6,5) * t114 + Ifges(6,3) * t32;
t105 = t65 ^ 2 + t69 ^ 2;
t102 = qJD(3) * t65;
t100 = qJD(3) * t70;
t99 = qJD(4) * t41;
t97 = qJD(4) * t68;
t96 = qJD(5) * t63;
t94 = qJD(5) * t69;
t89 = t47 * t96;
t87 = t131 * t64;
t85 = -t96 / 0.2e1;
t84 = -Ifges(6,6) * t63 - (2 * Ifges(5,4));
t49 = Ifges(6,5) * t95 - Ifges(6,6) * t96;
t82 = Ifges(6,1) * t67 - t118;
t81 = -Ifges(6,2) * t63 + t117;
t80 = Ifges(6,5) * t63 + Ifges(6,6) * t67;
t25 = -mrSges(6,2) * t46 - mrSges(6,3) * t111;
t26 = mrSges(6,1) * t46 - mrSges(6,3) * t110;
t79 = -t25 * t67 + t26 * t63;
t78 = -t25 * t63 - t26 * t67;
t77 = t38 * t67 + t39 * t63;
t75 = t89 + t114;
t50 = t81 * qJD(5);
t51 = t82 * qJD(5);
t55 = Ifges(6,2) * t67 + t118;
t56 = Ifges(6,1) * t63 + t117;
t74 = t67 * t50 + t63 * t51 - t55 * t96 + t56 * t95;
t18 = Ifges(6,6) * t46 + t47 * t81;
t19 = Ifges(6,5) * t46 + t47 * t82;
t3 = -Ifges(6,4) * t75 - Ifges(6,2) * t76 + Ifges(6,6) * t32;
t4 = -Ifges(6,1) * t75 - Ifges(6,4) * t76 + Ifges(6,5) * t32;
t73 = t51 * t110 / 0.2e1 - t50 * t111 / 0.2e1 + t19 * t95 / 0.2e1 + t63 * t4 / 0.2e1 + t18 * t85 + t67 * t3 / 0.2e1 - Ifges(5,5) * t31 + t46 * t49 / 0.2e1 + (t47 * t85 - t114 / 0.2e1) * t56 + (-Ifges(5,6) + t80 / 0.2e1) * t32 - t76 * t55 / 0.2e1;
t72 = -t14 * mrSges(5,2) + t41 * t48 + t131 * t15 + (-qJD(5) * t77 + t130) * mrSges(6,3);
t71 = pkin(2) ^ 2;
t35 = mrSges(5,1) * t46 + mrSges(5,2) * t47;
t23 = t83 * t47;
t12 = mrSges(5,1) * t32 - mrSges(5,2) * t31;
t10 = -mrSges(6,2) * t32 - mrSges(6,3) * t76;
t9 = mrSges(6,1) * t32 + mrSges(6,3) * t75;
t5 = mrSges(6,1) * t76 - mrSges(6,2) * t75;
t1 = [0.2e1 * m(6) * (t38 * t8 + t39 * t7 + t116) + 0.2e1 * m(4) * (-0.1e1 + t105) * t127 + 0.2e1 * (-t14 * t42 + t116 - t127) * m(5); t39 * t10 - t70 * t12 + t15 * t23 + t7 * t25 + t8 * t26 + t38 * t9 + t41 * t5 - (mrSges(4,1) * t65 + mrSges(4,2) * t69) * t100 + (-t14 * t46 + t15 * t47 - t31 * t41 + t32 * t42) * mrSges(5,3) + ((mrSges(4,3) * t105 - mrSges(3,2)) * t70 + (-mrSges(4,1) * t69 + mrSges(4,2) * t65 - mrSges(3,1) + t35) * t66) * qJD(2) + (m(6) * (t77 * t102 + (-t126 * qJD(5) - t63 * t7 - t67 * t8) * t69) + (-t100 * t65 - t104 * t69) * m(5)) * pkin(2); -(-t18 * t63 + t19 * t67) * t31 + (-0.2e1 * Ifges(4,4) * t65 + (t35 - t78) * t132) * t102 + (((2 * Ifges(5,2)) + Ifges(6,3)) * t32 - t84 * t31 + t107) * t46 + (-0.2e1 * Ifges(5,1) * t31 - t63 * t3 + t67 * t4 + (Ifges(6,5) * t67 + t84) * t32 + (-t67 * t18 - t63 * t19 - t46 * t80) * qJD(5)) * t47 + ((qJD(5) * t79 - t10 * t63 - t67 * t9 - t12) * t132 + 0.2e1 * (Ifges(4,4) * t69 + (-Ifges(4,2) + Ifges(4,1) + (-m(6) * t106 - m(5)) * t71) * t65) * qJD(3)) * t69; -t128 * mrSges(4,2) + (-t101 * t69 - t92) * mrSges(4,1) + (((t126 * qJD(4) - t15) * t122 + (-qJD(4) * t42 - t15) * t123) * t68 + ((-t38 * t95 - t39 * t96 + t130 + t99) * t122 + (t14 + t99) * t123) * t64) * t132 + t72; t73 + (Ifges(4,5) * t69 - Ifges(4,6) * t65) * qJD(3) + ((mrSges(5,3) * t31 - t5 + (-mrSges(5,3) * t46 - t79) * qJD(4)) * t68 + (qJD(4) * t23 + t10 * t67 - t63 * t9 + t78 * qJD(5) + (-t32 + t98) * mrSges(5,3)) * t64) * pkin(2); -0.2e1 * pkin(2) * t109 + (t87 * t132 + (0.2e1 * m(6) * (-0.1e1 + t106) * t71 * t64 + t129 * t132) * t68) * qJD(4) + t74; t72; t73; (-t109 + (t129 * t68 + t87) * qJD(4)) * pkin(2) + t74; t74; mrSges(6,1) * t8 - mrSges(6,2) * t7; -Ifges(6,5) * t89 - t76 * Ifges(6,6) + ((-t102 * t63 + t67 * t94) * mrSges(6,2) + (t102 * t67 + t63 * t94) * mrSges(6,1)) * pkin(2) + t107; ((t64 * t96 - t67 * t97) * mrSges(6,2) + (-t63 * t97 - t64 * t95) * mrSges(6,1)) * pkin(2) + t49; t49; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
