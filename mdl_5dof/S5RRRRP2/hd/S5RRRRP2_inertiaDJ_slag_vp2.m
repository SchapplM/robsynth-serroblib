% Calculate time derivative of joint inertia matrix for
% S5RRRRP2
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:06
% EndTime: 2020-01-03 12:11:09
% DurationCPUTime: 1.16s
% Computational Cost: add. (1701->203), mult. (3841->256), div. (0->0), fcn. (3035->6), ass. (0->100)
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t156 = t96 ^ 2 + t99 ^ 2;
t98 = cos(qJ(4));
t153 = pkin(3) * t98;
t152 = -mrSges(5,1) - mrSges(6,1);
t97 = sin(qJ(2));
t85 = pkin(1) * t97 + pkin(7);
t137 = -pkin(8) - t85;
t68 = t137 * t96;
t92 = t99 * pkin(8);
t69 = t85 * t99 + t92;
t95 = sin(qJ(4));
t30 = t95 * t68 + t98 * t69;
t141 = -pkin(8) - pkin(7);
t83 = t141 * t96;
t84 = pkin(7) * t99 + t92;
t47 = t95 * t83 + t98 * t84;
t151 = pkin(3) * qJD(4);
t100 = cos(qJ(2));
t150 = t156 * t100;
t149 = qJD(3) + qJD(4);
t81 = -mrSges(4,1) * t99 + mrSges(4,2) * t96;
t148 = (t81 - mrSges(3,1)) * t97 + (t156 * mrSges(4,3) - mrSges(3,2)) * t100;
t120 = qJD(4) * t98;
t121 = qJD(4) * t95;
t122 = qJD(3) * t99;
t72 = t95 * t99 + t96 * t98;
t40 = t149 * t72;
t71 = -t95 * t96 + t98 * t99;
t147 = Ifges(4,5) * t122 + (mrSges(5,3) + mrSges(6,3)) * (t120 * t71 + t121 * t72 - t40 * t95) * pkin(3);
t146 = 2 * m(5);
t145 = 2 * m(6);
t39 = t149 * t71;
t14 = t40 * mrSges(6,1) + t39 * mrSges(6,2);
t144 = 0.2e1 * t14;
t15 = mrSges(5,1) * t40 + mrSges(5,2) * t39;
t143 = 0.2e1 * t15;
t41 = -mrSges(6,1) * t71 + mrSges(6,2) * t72;
t142 = 0.2e1 * t41;
t42 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t140 = pkin(3) * t42;
t139 = t71 * pkin(4);
t138 = t99 * pkin(3);
t136 = Ifges(4,4) * t96;
t134 = Ifges(4,6) * t96;
t133 = t39 * mrSges(6,3);
t123 = qJD(3) * t96;
t90 = pkin(3) * t123;
t128 = pkin(1) * qJD(2);
t91 = t97 * t128;
t79 = t91 + t90;
t132 = t79 * t42;
t129 = -t40 * qJ(5) + t71 * qJD(5);
t127 = qJ(5) * t72;
t119 = 0.2e1 * mrSges(5,3);
t118 = 0.2e1 * mrSges(6,3);
t116 = t98 * t39 * mrSges(5,3);
t115 = t100 * t128;
t88 = -pkin(2) - t138;
t87 = -pkin(1) * t100 - pkin(2);
t28 = pkin(4) * t40 + t90;
t114 = qJD(3) * t141;
t113 = (-mrSges(5,2) - mrSges(6,2)) * t98;
t29 = t98 * t68 - t69 * t95;
t46 = t98 * t83 - t84 * t95;
t112 = qJD(3) * t137;
t111 = -(Ifges(5,6) + Ifges(6,6)) * t40 + (Ifges(5,5) + Ifges(6,5)) * t39;
t104 = -qJ(5) * t39 - qJD(5) * t72;
t48 = t96 * t112 + t99 * t115;
t49 = t99 * t112 - t96 * t115;
t9 = -qJD(4) * t30 - t48 * t95 + t98 * t49;
t3 = t104 + t9;
t107 = m(6) * t3 - t133;
t77 = t96 * t114;
t78 = t99 * t114;
t22 = -qJD(4) * t47 - t77 * t95 + t98 * t78;
t6 = t104 + t22;
t106 = m(6) * t6 - t133;
t105 = mrSges(4,1) * t96 + mrSges(4,2) * t99;
t80 = t87 - t138;
t8 = t68 * t120 - t69 * t121 + t98 * t48 + t95 * t49;
t21 = t83 * t120 - t84 * t121 + t98 * t77 + t95 * t78;
t2 = t8 + t129;
t103 = t9 * mrSges(5,1) + t3 * mrSges(6,1) - t8 * mrSges(5,2) - t2 * mrSges(6,2) + t111;
t102 = (Ifges(4,1) * t99 - t136) * t123 + (0.2e1 * Ifges(4,4) * t99 + (Ifges(4,1) - Ifges(4,2)) * t96) * t122 + 0.2e1 * (Ifges(6,4) + Ifges(5,4)) * (t39 * t71 - t40 * t72) - 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t71 * t40 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t72 * t39;
t5 = t21 + t129;
t101 = t22 * mrSges(5,1) + t6 * mrSges(6,1) - t21 * mrSges(5,2) - t5 * mrSges(6,2) + t111;
t86 = pkin(4) + t153;
t82 = Ifges(4,2) * t99 + t136;
t76 = t105 * qJD(3);
t66 = t71 * qJ(5);
t51 = t88 - t139;
t50 = t80 - t139;
t27 = t66 + t47;
t26 = t46 - t127;
t25 = t28 + t91;
t24 = t66 + t30;
t23 = t29 - t127;
t1 = [0.2e1 * t87 * t76 + 0.2e1 * t132 + t80 * t143 + t50 * t144 + t25 * t142 - t82 * t123 + (-t29 * t39 - t30 * t40 + t8 * t71 - t9 * t72) * t119 + 0.2e1 * (m(4) * (t150 * t85 + t87 * t97) + t148) * t128 + (t2 * t71 - t23 * t39 - t24 * t40 - t3 * t72) * t118 + t102 + (t29 * t9 + t30 * t8 + t79 * t80) * t146 + (t2 * t24 + t23 * t3 + t25 * t50) * t145; ((-t22 - t9) * t72 + (t21 + t8) * t71 - (t30 + t47) * t40 + (-t29 - t46) * t39) * mrSges(5,3) + ((-t3 - t6) * t72 + (t2 + t5) * t71 - (t24 + t27) * t40 + (-t23 - t26) * t39) * mrSges(6,3) + t132 + m(5) * (t21 * t30 + t22 * t29 + t46 * t9 + t47 * t8 + t79 * t88 + t80 * t90) + m(6) * (t2 * t27 + t23 * t6 + t24 * t5 + t25 * t51 + t26 * t3 + t28 * t50) + (t87 - pkin(2)) * t76 + (t28 + t25) * t41 + (m(4) * (-pkin(2) * t97 + pkin(7) * t150) + t148) * t128 + (t88 + t80) * t15 + (t50 + t51) * t14 + t102 + (-t82 + t140) * t123; (-t82 + 0.2e1 * t140) * t123 + (t21 * t47 + t22 * t46 + t88 * t90) * t146 + (t26 * t6 + t27 * t5 + t28 * t51) * t145 + t88 * t143 - 0.2e1 * pkin(2) * t76 + t51 * t144 + t28 * t142 + (-t26 * t39 - t27 * t40 + t5 * t71 - t6 * t72) * t118 + (t21 * t71 - t22 * t72 - t46 * t39 - t47 * t40) * t119 + t102; t107 * t86 + (-t116 + m(6) * (t24 * t120 - t23 * t121 + t2 * t95) + m(5) * (t30 * t120 - t29 * t121 + t8 * t95 + t9 * t98)) * pkin(3) + (t81 * t85 - t134) * qJD(3) - t105 * t115 + t103 + t147; t106 * t86 + (-t116 + m(6) * (t27 * t120 - t26 * t121 + t5 * t95) + m(5) * (t47 * t120 - t46 * t121 + t21 * t95 + t22 * t98)) * pkin(3) + (t81 * pkin(7) - t134) * qJD(3) + t101 + t147; 0.2e1 * (t113 + ((-t86 + t153) * m(6) + t152) * t95) * t151; t107 * pkin(4) + t103; t106 * pkin(4) + t101; (t113 + (-m(6) * pkin(4) + t152) * t95) * t151; 0; m(6) * t25 + t14; m(6) * t28 + t14; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
