% Calculate time derivative of joint inertia matrix for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:59
% EndTime: 2019-03-09 03:31:03
% DurationCPUTime: 2.06s
% Computational Cost: add. (1263->300), mult. (2580->406), div. (0->0), fcn. (1659->4), ass. (0->137)
t169 = Ifges(7,4) + Ifges(6,5);
t159 = 2 * qJD(2);
t81 = sin(qJ(5));
t83 = cos(qJ(5));
t137 = -t81 ^ 2 - t83 ^ 2;
t130 = qJD(5) * t83;
t84 = cos(qJ(3));
t133 = qJD(3) * t84;
t82 = sin(qJ(3));
t89 = t82 * t130 + t81 * t133;
t135 = qJD(3) * t82;
t124 = mrSges(6,1) * t135;
t21 = -mrSges(6,3) * t89 - t124;
t73 = mrSges(7,1) * t135;
t22 = mrSges(7,2) * t89 + t73;
t144 = t21 - t22;
t123 = mrSges(6,2) * t135;
t131 = qJD(5) * t82;
t119 = t81 * t131;
t120 = t83 * t133;
t90 = -t119 + t120;
t20 = mrSges(6,3) * t90 + t123;
t122 = mrSges(7,3) * t135;
t23 = mrSges(7,2) * t90 - t122;
t145 = t20 + t23;
t149 = t82 * t83;
t48 = -mrSges(6,2) * t84 + mrSges(6,3) * t149;
t49 = mrSges(7,2) * t149 + mrSges(7,3) * t84;
t140 = t48 + t49;
t151 = t81 * t82;
t46 = mrSges(6,1) * t84 - mrSges(6,3) * t151;
t47 = -mrSges(7,1) * t84 + mrSges(7,2) * t151;
t141 = t46 - t47;
t88 = t140 * t83 - t141 * t81;
t168 = t88 * qJD(5) + t144 * t83 + t145 * t81;
t167 = 2 * qJ(4);
t51 = t82 * pkin(3) - qJ(4) * t84 + qJ(2);
t34 = pkin(8) * t82 + t51;
t86 = -pkin(1) - pkin(7);
t158 = pkin(4) - t86;
t53 = t158 * t84;
t165 = t83 * t34 + t81 * t53;
t163 = m(7) * qJ(6) + mrSges(7,3);
t127 = qJ(4) * qJD(3);
t113 = pkin(3) * t133 + t82 * t127 + qJD(2);
t19 = (pkin(8) * qJD(3) - qJD(4)) * t84 + t113;
t38 = t158 * t135;
t4 = -qJD(5) * t165 - t19 * t81 - t38 * t83;
t29 = -qJD(4) * t84 + t113;
t161 = -0.2e1 * t29;
t160 = 0.2e1 * t51;
t157 = Ifges(6,4) * t81;
t156 = Ifges(6,4) * t83;
t155 = Ifges(7,5) * t81;
t154 = Ifges(7,5) * t83;
t153 = Ifges(6,6) * t84;
t152 = Ifges(7,6) * t84;
t85 = -pkin(3) - pkin(8);
t150 = t81 * t85;
t148 = t83 * t85;
t147 = mrSges(5,3) - mrSges(4,2);
t146 = -Ifges(7,2) - Ifges(6,3);
t101 = -Ifges(7,3) * t83 + t155;
t24 = t101 * t82 + t152;
t102 = Ifges(6,2) * t83 + t157;
t25 = t102 * t82 + t153;
t143 = t24 - t25;
t103 = Ifges(7,1) * t81 - t154;
t26 = t84 * Ifges(7,4) + t103 * t82;
t104 = Ifges(6,1) * t81 + t156;
t27 = t84 * Ifges(6,5) + t104 * t82;
t142 = t26 + t27;
t139 = t137 * t85 * t135;
t138 = qJD(4) * t82 + t84 * t127;
t136 = qJD(3) * t81;
t134 = qJD(3) * t83;
t132 = qJD(5) * t81;
t129 = qJD(5) * t84;
t128 = qJD(6) * t81;
t126 = 0.2e1 * t84;
t125 = pkin(5) * t135;
t116 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t58 = Ifges(7,3) * t81 + t154;
t59 = -Ifges(6,2) * t81 + t156;
t115 = t58 / 0.2e1 - t59 / 0.2e1;
t60 = Ifges(7,1) * t83 + t155;
t61 = Ifges(6,1) * t83 - t157;
t114 = t60 / 0.2e1 + t61 / 0.2e1;
t112 = m(5) * t86 - mrSges(5,1);
t99 = pkin(5) * t81 - qJ(6) * t83;
t50 = qJ(4) + t99;
t56 = t81 * mrSges(7,1) - t83 * mrSges(7,3);
t111 = m(7) * t50 + t56;
t110 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(5,2);
t3 = t53 * t130 - t132 * t34 + t83 * t19 - t81 * t38;
t91 = qJ(6) * t135 - qJD(6) * t84;
t1 = t3 - t91;
t2 = t125 - t4;
t109 = t1 * t81 - t2 * t83;
t108 = t3 * t81 + t4 * t83;
t107 = (m(7) * t85 - mrSges(7,2)) * t81;
t106 = mrSges(6,1) * t83 - mrSges(6,2) * t81;
t57 = t81 * mrSges(6,1) + t83 * mrSges(6,2);
t105 = mrSges(7,1) * t83 + mrSges(7,3) * t81;
t100 = pkin(5) * t83 + qJ(6) * t81;
t10 = qJ(6) * t84 + t165;
t14 = -t34 * t81 + t53 * t83;
t11 = -pkin(5) * t84 - t14;
t98 = t10 * t83 + t11 * t81;
t97 = -t14 * t81 + t165 * t83;
t94 = (m(6) / 0.2e1 + m(7) / 0.2e1) * t135;
t93 = Ifges(6,6) * t120 + Ifges(7,6) * t119 + t169 * t89;
t92 = -pkin(4) - t100;
t87 = -m(7) * t99 - t56 - t57;
t75 = t82 * t86;
t74 = Ifges(7,6) * t130;
t70 = t86 * t133;
t52 = -pkin(4) * t82 + t75;
t45 = t104 * qJD(5);
t44 = t103 * qJD(5);
t43 = t102 * qJD(5);
t42 = t101 * qJD(5);
t41 = t106 * qJD(5);
t40 = t105 * qJD(5);
t39 = -pkin(4) * t133 + t70;
t36 = t106 * t82;
t35 = t105 * t82;
t28 = qJD(5) * t100 - qJD(6) * t83 + qJD(4);
t18 = t82 * t92 + t75;
t13 = -mrSges(6,1) * t90 + mrSges(6,2) * t89;
t12 = -mrSges(7,1) * t90 - mrSges(7,3) * t89;
t9 = t61 * t131 + (-t82 * Ifges(6,5) + t104 * t84) * qJD(3);
t8 = t60 * t131 + (-t82 * Ifges(7,4) + t103 * t84) * qJD(3);
t7 = t59 * t131 + (-t82 * Ifges(6,6) + t102 * t84) * qJD(3);
t6 = t58 * t131 + (-t82 * Ifges(7,6) + t101 * t84) * qJD(3);
t5 = t70 + (qJD(5) * t99 - t128) * t82 + t92 * t133;
t15 = [m(5) * t29 * t160 + 0.2e1 * t1 * t49 + 0.2e1 * t10 * t23 + 0.2e1 * t11 * t22 + 0.2e1 * t18 * t12 + 0.2e1 * t52 * t13 + 0.2e1 * t14 * t21 + 0.2e1 * t165 * t20 + 0.2e1 * t2 * t47 + 0.2e1 * t3 * t48 - 0.2e1 * t5 * t35 - 0.2e1 * t39 * t36 + 0.2e1 * t4 * t46 + 0.2e1 * m(6) * (t14 * t4 + t165 * t3 + t39 * t52) + 0.2e1 * m(7) * (t1 * t10 + t11 * t2 + t18 * t5) + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t159 + ((mrSges(4,2) * t159) + mrSges(5,3) * t161 + (0.2e1 * qJ(2) * mrSges(4,1) - 0.2e1 * t51 * mrSges(5,2) + (-Ifges(4,4) - Ifges(5,6)) * t126 + (-t143 - t152) * t83 + t142 * t81) * qJD(3) + t93) * t84 + (mrSges(4,1) * t159 + mrSges(5,2) * t161 + (-t6 + t7) * t83 + (t8 + t9) * t81 + (t142 * t83 + (t143 - t153) * t81) * qJD(5) + (mrSges(5,3) * t160 - 0.2e1 * qJ(2) * mrSges(4,2) + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + (-Ifges(6,6) + Ifges(7,6)) * t83 - t169 * t81) * t82 + (-Ifges(4,1) - Ifges(5,2) + Ifges(5,3) + Ifges(4,2) + t146) * t126) * qJD(3)) * t82; (t12 + t13 + (t140 * t81 + t141 * t83) * qJD(3) + m(7) * (t10 * t136 - t11 * t134 + t5) + m(6) * (t134 * t14 + t136 * t165 + t39)) * t82 + ((-t35 - t36) * qJD(3) + m(7) * (qJD(3) * t18 - t10 * t130 - t11 * t132 - t109) + m(6) * (qJD(3) * t52 - t130 * t165 + t132 * t14 - t108) - t168) * t84; 0.4e1 * (0.1e1 + t137) * t84 * t94; t84 * t74 / 0.2e1 + t18 * t40 + t50 * t12 + t52 * t41 + t5 * t56 + t39 * t57 - t28 * t35 - qJD(4) * t36 + qJ(4) * t13 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t144 * t85) * t83 + (t6 / 0.2e1 - t7 / 0.2e1 - t1 * mrSges(7,2) - t3 * mrSges(6,3) + t145 * t85) * t81 + m(7) * (t1 * t150 - t148 * t2 + t18 * t28 + t5 * t50) + m(6) * (qJ(4) * t39 + qJD(4) * t52 + t148 * t4 + t150 * t3) + ((t42 / 0.2e1 - t43 / 0.2e1) * t83 + (-t44 / 0.2e1 - t45 / 0.2e1) * t81 + t112 * qJD(4)) * t82 + ((-t10 * mrSges(7,2) - t165 * mrSges(6,3) - t153 / 0.2e1 + t24 / 0.2e1 - t25 / 0.2e1 + t114 * t82) * t83 + (-t26 / 0.2e1 - t27 / 0.2e1 - t11 * mrSges(7,2) + t14 * mrSges(6,3) + t116 * t84 + t115 * t82) * t81 + (m(6) * t97 + m(7) * t98 + t88) * t85) * qJD(5) + ((pkin(3) * mrSges(5,1) + Ifges(5,4) - Ifges(4,5) + t116 * t83 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t81 + t110 * t86) * t82 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) - t115 * t83 + t114 * t81 + (m(5) * qJ(4) + t147) * t86) * t84) * qJD(3); (t40 + t41) * t82 + m(6) * (t138 - t139) + m(7) * (t28 * t82 - t139) + m(5) * t138 + ((t111 + t57 + t147) * t84 + (t110 + (mrSges(7,2) + mrSges(6,3)) * t137) * t82) * qJD(3); t41 * t167 + 0.2e1 * t40 * t50 + (-t44 - t45) * t83 + (-t42 + t43) * t81 + ((t58 - t59) * t83 + (-t60 - t61) * t81) * qJD(5) + 0.2e1 * t111 * t28 + (0.2e1 * mrSges(5,3) + 0.2e1 * t57 + (m(5) + m(6)) * t167) * qJD(4); t112 * t135 + m(7) * (qJD(5) * t98 + t109) + m(6) * (qJD(5) * t97 + t108) + t168; m(5) * t135 - 0.2e1 * t137 * t94; 0; 0; -Ifges(6,6) * t119 - pkin(5) * t22 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t10) + qJD(6) * t49 + qJ(6) * t23 + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) + (t146 * t82 - t83 * t152) * qJD(3) + t93; (m(7) * t125 + t124 + t73) * t83 + (m(7) * t91 + t122 - t123) * t81 + ((mrSges(6,2) - t163) * t83 + (m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1)) * t81) * t129; t74 + qJD(6) * t107 + ((-qJ(6) * mrSges(7,2) - Ifges(6,6)) * t83 + (pkin(5) * mrSges(7,2) - t169) * t81 + t87 * t85) * qJD(5); m(7) * t128 + qJD(5) * t87; 0.2e1 * t163 * qJD(6); m(7) * t2 + t22; (-t129 * t81 - t134 * t82) * m(7); qJD(5) * t107; m(7) * t132; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
