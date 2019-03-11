% Calculate time derivative of joint inertia matrix for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:37:00
% EndTime: 2019-03-09 03:37:04
% DurationCPUTime: 2.25s
% Computational Cost: add. (4746->362), mult. (10067->547), div. (0->0), fcn. (9582->10), ass. (0->148)
t111 = sin(qJ(5));
t144 = qJD(5) * t111;
t107 = sin(pkin(11));
t108 = cos(pkin(11));
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t89 = t107 * t115 + t108 * t112;
t136 = t89 * t144;
t114 = cos(qJ(5));
t88 = t107 * t112 - t108 * t115;
t83 = t88 * qJD(3);
t147 = t114 * t83;
t117 = t136 + t147;
t143 = qJD(5) * t114;
t118 = -t111 * t83 + t143 * t89;
t31 = t118 * mrSges(6,1) - t117 * mrSges(6,2);
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t120 = t110 * t111 - t113 * t114;
t169 = qJD(5) + qJD(6);
t91 = t110 * t114 + t111 * t113;
t68 = t169 * t91;
t18 = t120 * t83 - t68 * t89;
t56 = t120 * t89;
t19 = t169 * t56 + t91 * t83;
t6 = -mrSges(7,1) * t19 + t18 * mrSges(7,2);
t172 = t31 + t6;
t67 = t169 * t120;
t34 = mrSges(7,1) * t68 - t67 * mrSges(7,2);
t126 = mrSges(6,1) * t111 + mrSges(6,2) * t114;
t92 = t126 * qJD(5);
t171 = t34 + t92;
t82 = t89 * qJD(3);
t170 = -Ifges(6,5) * t147 + Ifges(6,3) * t82;
t168 = 2 * m(7);
t167 = -2 * mrSges(5,3);
t102 = sin(pkin(10)) * pkin(1) + pkin(7);
t145 = qJ(4) + t102;
t128 = qJD(3) * t145;
t116 = -qJD(4) * t112 - t115 * t128;
t74 = qJD(4) * t115 - t112 * t128;
t43 = t107 * t74 - t108 * t116;
t166 = 0.2e1 * t43;
t129 = t145 * t112;
t84 = t145 * t115;
t57 = t107 * t84 + t108 * t129;
t165 = 0.2e1 * t57;
t135 = -cos(pkin(10)) * pkin(1) - pkin(2);
t96 = -pkin(3) * t115 + t135;
t164 = 0.2e1 * t96;
t163 = m(5) * pkin(3);
t162 = m(7) * pkin(5);
t160 = -t89 / 0.2e1;
t151 = Ifges(6,4) * t111;
t98 = Ifges(6,2) * t114 + t151;
t159 = -t98 / 0.2e1;
t158 = t43 * t57;
t157 = t82 * Ifges(6,5);
t156 = t82 * Ifges(6,6);
t155 = t88 * Ifges(6,6);
t69 = t88 * t82;
t101 = pkin(3) * t107 + pkin(8);
t154 = pkin(9) + t101;
t58 = -t107 * t129 + t108 * t84;
t53 = t114 * t58;
t54 = pkin(4) * t88 - pkin(8) * t89 + t96;
t28 = t111 * t54 + t53;
t153 = -Ifges(7,5) * t67 - Ifges(7,6) * t68;
t97 = -mrSges(6,1) * t114 + mrSges(6,2) * t111;
t152 = t97 - mrSges(5,1);
t150 = Ifges(6,4) * t114;
t149 = Ifges(6,6) * t111;
t148 = t111 * t89;
t146 = t114 * t89;
t142 = qJD(6) * t110;
t141 = qJD(6) * t113;
t140 = 0.2e1 * t115;
t139 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t82;
t138 = pkin(3) * qJD(3) * t112;
t137 = pkin(5) * t144;
t103 = -pkin(3) * t108 - pkin(4);
t77 = t83 * mrSges(5,2);
t134 = t82 * mrSges(5,1) - t77;
t133 = -(2 * Ifges(5,4)) - t149;
t44 = t107 * t116 + t108 * t74;
t52 = pkin(4) * t82 + pkin(8) * t83 + t138;
t132 = -t111 * t44 + t114 * t52;
t27 = -t111 * t58 + t114 * t54;
t131 = qJD(5) * t154;
t130 = (t111 ^ 2 + t114 ^ 2) * t83;
t127 = t43 * t88 + t57 * t82;
t125 = Ifges(6,1) * t114 - t151;
t124 = -Ifges(6,2) * t111 + t150;
t20 = pkin(5) * t88 - pkin(9) * t146 + t27;
t21 = -pkin(9) * t148 + t28;
t9 = -t110 * t21 + t113 * t20;
t10 = t110 * t20 + t113 * t21;
t85 = t154 * t111;
t86 = t154 * t114;
t60 = -t110 * t86 - t113 * t85;
t61 = -t110 * t85 + t113 * t86;
t123 = -t111 * t27 + t114 * t28;
t62 = -mrSges(6,2) * t88 - mrSges(6,3) * t148;
t63 = mrSges(6,1) * t88 - mrSges(6,3) * t146;
t122 = -t111 * t63 + t114 * t62;
t80 = t111 * t131;
t81 = t114 * t131;
t32 = qJD(6) * t60 - t110 * t81 - t113 * t80;
t33 = -qJD(6) * t61 + t110 * t80 - t113 * t81;
t121 = mrSges(7,1) * t33 - t32 * mrSges(7,2) + t153;
t7 = pkin(9) * t147 + pkin(5) * t82 + (-t53 + (pkin(9) * t89 - t54) * t111) * qJD(5) + t132;
t11 = t111 * t52 + t114 * t44 + t143 * t54 - t144 * t58;
t8 = -pkin(9) * t118 + t11;
t2 = qJD(6) * t9 + t110 * t7 + t113 * t8;
t3 = -qJD(6) * t10 - t110 * t8 + t113 * t7;
t119 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t139;
t104 = Ifges(6,5) * t143;
t99 = Ifges(6,1) * t111 + t150;
t95 = -pkin(5) * t114 + t103;
t94 = t125 * qJD(5);
t93 = t124 * qJD(5);
t87 = (-mrSges(7,1) * t110 - mrSges(7,2) * t113) * qJD(6) * pkin(5);
t72 = Ifges(7,1) * t91 - Ifges(7,4) * t120;
t71 = Ifges(7,4) * t91 - Ifges(7,2) * t120;
t70 = mrSges(7,1) * t120 + mrSges(7,2) * t91;
t59 = t126 * t89;
t55 = t91 * t89;
t46 = Ifges(6,5) * t88 + t125 * t89;
t45 = t124 * t89 + t155;
t42 = pkin(5) * t148 + t57;
t40 = mrSges(7,1) * t88 + mrSges(7,3) * t56;
t39 = -mrSges(7,2) * t88 - mrSges(7,3) * t55;
t38 = -mrSges(6,2) * t82 - mrSges(6,3) * t118;
t37 = mrSges(6,1) * t82 + mrSges(6,3) * t117;
t36 = -Ifges(7,1) * t67 - Ifges(7,4) * t68;
t35 = -Ifges(7,4) * t67 - Ifges(7,2) * t68;
t29 = mrSges(7,1) * t55 - mrSges(7,2) * t56;
t26 = pkin(5) * t118 + t43;
t25 = -Ifges(7,1) * t56 - Ifges(7,4) * t55 + Ifges(7,5) * t88;
t24 = -Ifges(7,4) * t56 - Ifges(7,2) * t55 + Ifges(7,6) * t88;
t23 = -Ifges(6,1) * t117 - Ifges(6,4) * t118 + t157;
t22 = -Ifges(6,4) * t117 - Ifges(6,2) * t118 + t156;
t14 = -mrSges(7,2) * t82 + mrSges(7,3) * t19;
t13 = mrSges(7,1) * t82 - mrSges(7,3) * t18;
t12 = -qJD(5) * t28 + t132;
t5 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t82 * Ifges(7,5);
t4 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t82 * Ifges(7,6);
t1 = [(t58 * t167 + (Ifges(6,5) * t114 + t133) * t89 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t88 - Ifges(7,5) * t56 - Ifges(7,6) * t55) * t82 + (mrSges(5,3) * t166 - 0.2e1 * Ifges(5,1) * t83 - t111 * t22 + t114 * t23 + (-t114 * t45 - t111 * t46 + t88 * (-Ifges(6,5) * t111 - Ifges(6,6) * t114)) * qJD(5)) * t89 + (t10 * t2 + t26 * t42 + t3 * t9) * t168 + t134 * t164 + t31 * t165 + t59 * t166 - (mrSges(5,3) * t165 - t111 * t45 + t114 * t46) * t83 + ((mrSges(4,2) * t135 + Ifges(4,4) * t115) * t140 + (0.2e1 * pkin(3) * (mrSges(5,1) * t88 + mrSges(5,2) * t89) + 0.2e1 * t135 * mrSges(4,1) + t163 * t164 - 0.2e1 * Ifges(4,4) * t112 + (Ifges(4,1) - Ifges(4,2)) * t140) * t112) * qJD(3) + 0.2e1 * m(5) * (t44 * t58 + t158) + 0.2e1 * m(6) * (t11 * t28 + t12 * t27 + t158) + (-t133 * t83 + t44 * t167 + t139 + t170) * t88 + 0.2e1 * t11 * t62 + 0.2e1 * t12 * t63 - t55 * t4 - t56 * t5 + 0.2e1 * t27 * t37 + 0.2e1 * t28 * t38 + 0.2e1 * t2 * t39 + 0.2e1 * t3 * t40 + 0.2e1 * t42 * t6 + t19 * t24 + t18 * t25 + 0.2e1 * t26 * t29 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14; -t55 * t13 - t56 * t14 + t18 * t39 + t19 * t40 + t172 * t88 - t122 * t83 + (t29 + t59) * t82 + (-t111 * t37 + t114 * t38 + (-t111 * t62 - t114 * t63) * qJD(5)) * t89 + m(7) * (t10 * t18 + t19 * t9 - t2 * t56 + t26 * t88 - t3 * t55 + t42 * t82) + m(6) * (-t123 * t83 + (t11 * t114 - t111 * t12 + (-t111 * t28 - t114 * t27) * qJD(5)) * t89 + t127) + m(5) * (t44 * t89 - t58 * t83 + t127); 0.2e1 * m(7) * (-t18 * t56 - t19 * t55 + t69) + 0.2e1 * m(6) * (-t130 * t89 + t69) + 0.2e1 * m(5) * (-t83 * t89 + t69); (t157 / 0.2e1 + t23 / 0.2e1 - t12 * mrSges(6,3) + t93 * t160 - t83 * t159 + (t99 * t160 - t28 * mrSges(6,3) + pkin(5) * t29 - t155 / 0.2e1 - t45 / 0.2e1 + t42 * t162) * qJD(5) + (-m(6) * t12 - t37 + (-m(6) * t28 - t62) * qJD(5)) * t101) * t111 + (t156 / 0.2e1 + t22 / 0.2e1 + t11 * mrSges(6,3) + t89 * t94 / 0.2e1 - t83 * t99 / 0.2e1 + (t89 * t159 - t27 * mrSges(6,3) + t46 / 0.2e1) * qJD(5) + (m(6) * (-t27 * qJD(5) + t11) + t38 - qJD(5) * t63) * t101) * t114 + (m(5) * (t107 * t44 - t108 * t43) + (-t107 * t82 + t108 * t83) * mrSges(5,3)) * pkin(3) + (Ifges(4,5) * t115 - Ifges(4,6) * t112 + (-mrSges(4,1) * t115 + mrSges(4,2) * t112) * t102) * qJD(3) + (m(6) * t103 + t152) * t43 + t91 * t5 / 0.2e1 + t57 * t92 + t95 * t6 + t103 * t31 - Ifges(5,6) * t82 - Ifges(5,5) * t83 - t67 * t25 / 0.2e1 - t68 * t24 / 0.2e1 + t26 * t70 + t19 * t71 / 0.2e1 + t18 * t72 / 0.2e1 + t60 * t13 + t61 * t14 - t55 * t35 / 0.2e1 - t56 * t36 / 0.2e1 + t32 * t39 + t33 * t40 + t42 * t34 - t44 * mrSges(5,2) + m(7) * (t10 * t32 + t2 * t61 + t26 * t95 + t3 * t60 + t33 * t9) + (t153 + t104) * t88 / 0.2e1 + t82 * (Ifges(7,5) * t91 - Ifges(7,6) * t120) / 0.2e1 - t120 * t4 / 0.2e1 + (-t10 * t68 - t120 * t2 - t3 * t91 + t9 * t67) * mrSges(7,3); t77 + t171 * t88 + (-mrSges(4,1) * t112 - mrSges(4,2) * t115) * qJD(3) - mrSges(6,3) * t130 + (t70 + t152) * t82 + m(7) * (t137 * t88 + t18 * t61 + t19 * t60 - t32 * t56 - t33 * t55 + t82 * t95) + m(6) * (-t101 * t130 + t103 * t82) + (-t107 * t83 - t108 * t82) * t163 + (-t120 * t18 - t19 * t91 - t55 * t67 + t56 * t68) * mrSges(7,3); -t67 * t72 + t91 * t36 - t68 * t71 - t120 * t35 + (t137 * t95 + t32 * t61 + t33 * t60) * t168 + 0.2e1 * t70 * t137 + 0.2e1 * t95 * t34 + 0.2e1 * t103 * t92 + t111 * t94 - t98 * t144 + (qJD(5) * t99 + t93) * t114 + 0.2e1 * (-t120 * t32 - t33 * t91 + t60 * t67 - t61 * t68) * mrSges(7,3); m(5) * t138 + t111 * t38 + t114 * t37 - t120 * t13 + t91 * t14 - t67 * t39 - t68 * t40 + t122 * qJD(5) + m(7) * (-t10 * t67 - t120 * t3 + t2 * t91 - t68 * t9) + m(6) * (qJD(5) * t123 + t11 * t111 + t114 * t12) + t134; m(7) * (-t120 * t19 + t18 * t91 + t55 * t68 + t56 * t67); m(7) * (-t120 * t33 + t32 * t91 - t60 * t68 - t61 * t67); (t120 * t68 - t67 * t91) * t168; -Ifges(6,5) * t136 + t12 * mrSges(6,1) - t11 * mrSges(6,2) - t118 * Ifges(6,6) + (m(7) * (t10 * t141 + t110 * t2 + t113 * t3 - t142 * t9) + t39 * t141 + t110 * t14 - t40 * t142 + t113 * t13) * pkin(5) + t119 + t170; (t110 * t18 + t113 * t19 + (t110 * t55 - t113 * t56) * qJD(6)) * t162 - t172; t104 + (t101 * t97 - t149) * qJD(5) + (m(7) * (t110 * t32 + t113 * t33 + (-t110 * t60 + t113 * t61) * qJD(6)) + (-t110 * t68 + t113 * t67 + (t110 * t91 - t113 * t120) * qJD(6)) * mrSges(7,3)) * pkin(5) + t121; (-t110 * t67 - t113 * t68 + (t110 * t120 + t113 * t91) * qJD(6)) * t162 - t171; 0.2e1 * t87; t119; -t6; t121; -t34; t87; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
