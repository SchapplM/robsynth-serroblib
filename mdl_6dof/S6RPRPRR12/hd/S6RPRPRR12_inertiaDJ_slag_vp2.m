% Calculate time derivative of joint inertia matrix for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:14
% EndTime: 2019-03-09 04:18:21
% DurationCPUTime: 2.58s
% Computational Cost: add. (2737->351), mult. (5537->523), div. (0->0), fcn. (4371->6), ass. (0->156)
t111 = sin(qJ(5));
t112 = sin(qJ(3));
t114 = cos(qJ(5));
t154 = qJD(5) * t114;
t115 = cos(qJ(3));
t157 = qJD(3) * t115;
t121 = t111 * t157 + t112 * t154;
t117 = -pkin(1) - pkin(7);
t179 = pkin(4) - t117;
t88 = mrSges(6,1) * t111 + mrSges(6,2) * t114;
t195 = mrSges(5,3) + t88;
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t162 = t113 * t114;
t125 = t110 * t111 - t162;
t59 = t125 * t112;
t194 = -t112 * pkin(3) + qJ(4) * t115;
t193 = t112 * mrSges(4,1) + t115 * mrSges(4,2);
t156 = qJD(5) * t111;
t159 = qJD(3) * t112;
t145 = pkin(3) * t157 + qJ(4) * t159 + qJD(2);
t52 = (pkin(8) * qJD(3) - qJD(4)) * t115 + t145;
t83 = qJ(2) - t194;
t67 = pkin(8) * t112 + t83;
t74 = t179 * t159;
t87 = t179 * t115;
t11 = -t111 * t74 + t114 * t52 + t87 * t154 - t156 * t67;
t140 = -t111 * t52 - t114 * t74;
t70 = t111 * t87;
t37 = t114 * t67 + t70;
t12 = -qJD(5) * t37 + t140;
t133 = t11 * t111 + t114 * t12;
t192 = qJD(5) + qJD(6);
t146 = t114 * t157;
t191 = Ifges(6,5) * t121 + Ifges(6,6) * t146;
t81 = -mrSges(6,3) * t111 * t112 + mrSges(6,1) * t115;
t163 = t112 * t114;
t82 = -mrSges(6,2) * t115 + mrSges(6,3) * t163;
t128 = t111 * t81 - t114 * t82;
t155 = qJD(5) * t112;
t122 = -t111 * t155 + t146;
t53 = mrSges(6,2) * t159 + mrSges(6,3) * t122;
t54 = -mrSges(6,1) * t159 - mrSges(6,3) * t121;
t190 = t128 * qJD(5) - t111 * t53 - t114 * t54;
t126 = t110 * t114 + t113 * t111;
t152 = qJD(6) * t110;
t41 = -t110 * t156 - t111 * t152 + t162 * t192;
t42 = t192 * t126;
t189 = qJD(6) * (-t110 * t125 - t113 * t126) - t110 * t41 + t113 * t42;
t188 = t115 * t192;
t187 = 2 * m(7);
t57 = -qJD(4) * t115 + t145;
t186 = -0.2e1 * t57;
t185 = 0.2e1 * t83;
t184 = m(7) * pkin(5);
t183 = t111 / 0.2e1;
t182 = t114 / 0.2e1;
t181 = t115 / 0.2e1;
t180 = -Ifges(6,3) - Ifges(7,3);
t116 = -pkin(3) - pkin(8);
t178 = pkin(9) - t116;
t124 = qJD(3) * t126;
t24 = t115 * t124 - t192 * t59;
t123 = qJD(3) * t125;
t26 = -t112 * t42 - t115 * t123;
t177 = Ifges(7,5) * t24 + Ifges(7,6) * t26;
t176 = -Ifges(7,5) * t42 - Ifges(7,6) * t41;
t175 = Ifges(6,4) * t111;
t174 = Ifges(6,4) * t114;
t136 = Ifges(6,1) * t111 + t174;
t56 = t115 * Ifges(6,5) + t112 * t136;
t172 = t111 * t56;
t90 = Ifges(6,1) * t114 - t175;
t171 = t111 * t90;
t135 = Ifges(6,2) * t114 + t175;
t165 = t115 * Ifges(6,6);
t55 = t112 * t135 + t165;
t168 = t114 * t55;
t89 = -Ifges(6,2) * t111 + t174;
t167 = t114 * t89;
t161 = qJ(4) * t157 + qJD(4) * t112;
t160 = t111 ^ 2 + t114 ^ 2;
t158 = qJD(3) * t114;
t153 = qJD(5) * t115;
t151 = qJD(6) * t113;
t150 = 2 * mrSges(7,3);
t149 = 0.2e1 * t115;
t96 = t112 * t157;
t144 = m(5) * t117 - mrSges(5,1);
t25 = t112 * t124 + t125 * t188;
t27 = -t112 * t123 + t126 * t188;
t143 = t27 * mrSges(7,1) - t25 * mrSges(7,2);
t142 = -t42 * mrSges(7,1) - t41 * mrSges(7,2);
t141 = -pkin(9) * t112 - t67;
t86 = t178 * t114;
t75 = t179 * t157;
t139 = -t125 * t42 - t126 * t41;
t7 = (-pkin(9) * t111 * t115 - pkin(5) * t112) * qJD(3) + (t114 * t141 - t70) * qJD(5) + t140;
t8 = pkin(9) * t122 + t11;
t71 = t114 * t87;
t30 = pkin(5) * t115 + t111 * t141 + t71;
t34 = pkin(9) * t163 + t37;
t9 = -t110 * t34 + t113 * t30;
t2 = qJD(6) * t9 + t110 * t7 + t113 * t8;
t10 = t110 * t30 + t113 * t34;
t3 = -qJD(6) * t10 - t110 * t8 + t113 * t7;
t138 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t177;
t137 = mrSges(6,1) * t114 - mrSges(6,2) * t111;
t134 = -Ifges(6,5) * t111 - Ifges(6,6) * t114;
t84 = t178 * t111;
t48 = -t110 * t86 - t113 * t84;
t47 = t110 * t84 - t113 * t86;
t36 = -t111 * t67 + t71;
t130 = t111 * t36 - t114 * t37;
t72 = t178 * t156;
t73 = qJD(5) * t86;
t16 = qJD(6) * t47 + t110 * t72 - t113 * t73;
t17 = -qJD(6) * t48 + t110 * t73 + t113 * t72;
t127 = t17 * mrSges(7,1) - t16 * mrSges(7,2) + t176;
t120 = t10 * t41 - t125 * t3 + t126 * t2 - t42 * t9;
t119 = t125 * t17 - t126 * t16 - t41 * t48 + t42 * t47;
t60 = t125 * t115;
t62 = t126 * t115;
t118 = t125 * t27 - t126 * t25 + t41 * t62 + t42 * t60;
t104 = t112 * t117;
t99 = pkin(5) * t111 + qJ(4);
t94 = pkin(5) * t154 + qJD(4);
t85 = -pkin(4) * t112 + t104;
t80 = t136 * qJD(5);
t79 = t135 * qJD(5);
t78 = t137 * qJD(5);
t69 = (-mrSges(7,1) * t110 - mrSges(7,2) * t113) * qJD(6) * pkin(5);
t68 = t137 * t112;
t63 = t104 + (-pkin(5) * t114 - pkin(4)) * t112;
t61 = t126 * t112;
t51 = mrSges(7,1) * t115 - mrSges(7,3) * t61;
t50 = -mrSges(7,2) * t115 - mrSges(7,3) * t59;
t46 = -Ifges(7,1) * t125 - Ifges(7,4) * t126;
t45 = -Ifges(7,4) * t125 - Ifges(7,2) * t126;
t44 = mrSges(7,1) * t126 - mrSges(7,2) * t125;
t43 = -pkin(5) * t122 - t75;
t35 = -mrSges(6,1) * t122 + mrSges(6,2) * t121;
t33 = mrSges(7,1) * t59 + mrSges(7,2) * t61;
t32 = t90 * t155 + (-t112 * Ifges(6,5) + t115 * t136) * qJD(3);
t31 = t89 * t155 + (-t112 * Ifges(6,6) + t115 * t135) * qJD(3);
t29 = Ifges(7,1) * t61 - Ifges(7,4) * t59 + Ifges(7,5) * t115;
t28 = Ifges(7,4) * t61 - Ifges(7,2) * t59 + Ifges(7,6) * t115;
t20 = -Ifges(7,1) * t42 - Ifges(7,4) * t41;
t19 = -Ifges(7,4) * t42 - Ifges(7,2) * t41;
t18 = mrSges(7,1) * t41 - mrSges(7,2) * t42;
t14 = mrSges(7,2) * t159 + mrSges(7,3) * t26;
t13 = -mrSges(7,1) * t159 - mrSges(7,3) * t24;
t6 = -mrSges(7,1) * t26 + mrSges(7,2) * t24;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t26 - Ifges(7,5) * t159;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t26 - Ifges(7,6) * t159;
t1 = [-t59 * t4 + t61 * t5 + 0.2e1 * t63 * t6 + 0.2e1 * t43 * t33 + 0.2e1 * t2 * t50 + 0.2e1 * t3 * t51 + 0.2e1 * t37 * t53 + 0.2e1 * t36 * t54 + t26 * t28 + t24 * t29 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + m(5) * t57 * t185 + 0.2e1 * t75 * t68 + 0.2e1 * t12 * t81 + 0.2e1 * t11 * t82 + 0.2e1 * t85 * t35 + (t10 * t2 + t3 * t9 + t43 * t63) * t187 + 0.2e1 * m(6) * (t11 * t37 + t12 * t36 - t75 * t85) + (mrSges(5,3) * t186 + t177 + t191) * t115 + (mrSges(5,2) * t186 + t111 * t32 + t114 * t31 + (t114 * t56 + (-t55 - t165) * t111) * qJD(5)) * t112 + 0.2e1 * (mrSges(3,3) + (m(3) + m(4)) * qJ(2) + t193) * qJD(2) + ((0.2e1 * qJ(2) * mrSges(4,1) - 0.2e1 * t83 * mrSges(5,2) + t172 + t168 + (-Ifges(4,4) - Ifges(5,6)) * t149) * t115 + (-0.2e1 * qJ(2) * mrSges(4,2) + mrSges(5,3) * t185 - Ifges(7,5) * t61 + Ifges(7,6) * t59 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t134) * t112 + (-Ifges(4,1) + Ifges(4,2) - Ifges(5,2) + Ifges(5,3) + t180) * t149) * t112) * qJD(3); m(7) * (t10 * t25 - t2 * t62 + t27 * t9 + t3 * t60) + t25 * t50 - t62 * t14 + t27 * t51 + t60 * t13 + (t35 + t6 + (t111 * t82 + t114 * t81) * qJD(3) + m(7) * t43 + m(6) * (qJD(3) * t111 * t37 + t158 * t36 - t75)) * t112 + (m(6) * (-t154 * t37 + t156 * t36 - t133) + (m(6) * t85 + m(7) * t63 + t33 - t68) * qJD(3) + t190) * t115; 0.2e1 * m(6) * (-t160 + 0.1e1) * t96 + 0.2e1 * m(7) * (-t25 * t62 + t27 * t60 + t96); -t59 * t19 / 0.2e1 + t61 * t20 / 0.2e1 + t63 * t18 - t42 * t29 / 0.2e1 + t43 * t44 + t26 * t45 / 0.2e1 + t24 * t46 / 0.2e1 + t47 * t13 + t48 * t14 + t16 * t50 + t17 * t51 + qJ(4) * t35 - t41 * t28 / 0.2e1 - t120 * mrSges(7,3) + (t144 * qJD(4) - t79 * t182 - t80 * t183) * t112 - t126 * t4 / 0.2e1 + ((pkin(3) * mrSges(5,1) + Ifges(7,5) * t125 / 0.2e1 + Ifges(7,6) * t126 / 0.2e1 + Ifges(5,4) - Ifges(4,5) - Ifges(6,5) * t114 / 0.2e1 + Ifges(6,6) * t183) * t112 + (-qJ(4) * mrSges(5,1) + t171 / 0.2e1 + t167 / 0.2e1 + Ifges(5,5) - Ifges(4,6)) * t115 + (m(5) * t194 + t112 * mrSges(5,2) + t115 * mrSges(5,3) - t193) * t117) * qJD(3) - t125 * t5 / 0.2e1 + t176 * t181 + m(6) * (-qJ(4) * t75 + qJD(4) * t85 + t116 * t133) + m(7) * (t10 * t16 + t17 * t9 + t2 * t48 + t3 * t47 + t43 * t99 + t63 * t94) - qJD(4) * t68 + t85 * t78 - t75 * t88 + t94 * t33 + t99 * t6 + (-t172 / 0.2e1 - t168 / 0.2e1 + t134 * t181 + (-t111 * t89 / 0.2e1 + t90 * t182) * t112 + t130 * mrSges(6,3) + (-m(6) * t130 - t128) * t116) * qJD(5) + (-t31 / 0.2e1 - t11 * mrSges(6,3) + t116 * t53) * t111 + (-t12 * mrSges(6,3) + t116 * t54 + t32 / 0.2e1) * t114; (t18 + t78) * t112 + ((-mrSges(4,2) + t44 + t195) * t115 + (-mrSges(6,3) * t160 - mrSges(4,1) + mrSges(5,2)) * t112) * qJD(3) + m(6) * (t116 * t159 * t160 + t161) + m(7) * (t112 * t94 + t157 * t99 - t16 * t62 + t17 * t60 + t25 * t48 + t27 * t47) + m(5) * (-pkin(3) * t159 + t161) + t118 * mrSges(7,3); -t41 * t45 - t126 * t19 - t42 * t46 - t125 * t20 + (t16 * t48 + t17 * t47 + t94 * t99) * t187 + 0.2e1 * t94 * t44 + 0.2e1 * t99 * t18 + 0.2e1 * qJ(4) * t78 - t114 * t80 + t111 * t79 + (-t167 - t171) * qJD(5) + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t195) * qJD(4) + t119 * t150; -t125 * t13 + t126 * t14 + t41 * t50 - t42 * t51 + t144 * t159 + m(7) * t120 + m(6) * (-qJD(5) * t130 + t133) - t190; -m(7) * t118 + (m(6) * t160 + m(5)) * t159; -m(7) * t119 + t139 * t150; -0.2e1 * m(7) * t139; t12 * mrSges(6,1) - t11 * mrSges(6,2) + (-Ifges(6,6) * t156 + qJD(3) * t180) * t112 + (m(7) * (t10 * t151 + t110 * t2 + t113 * t3 - t152 * t9) + t50 * t151 + t110 * t14 - t51 * t152 + t113 * t13) * pkin(5) + t138 + t191; (-t111 * t159 + t114 * t153) * mrSges(6,2) + (t111 * t153 + t112 * t158) * mrSges(6,1) + (t110 * t25 + t113 * t27 + (-t110 * t60 - t113 * t62) * qJD(6)) * t184 + t143; ((-mrSges(6,2) * t116 - Ifges(6,6)) * t114 + (-mrSges(6,1) * t116 - Ifges(6,5)) * t111) * qJD(5) + (m(7) * (t110 * t16 + t113 * t17 + (-t110 * t47 + t113 * t48) * qJD(6)) + t189 * mrSges(7,3)) * pkin(5) + t127; -t88 * qJD(5) - t184 * t189 + t142; 0.2e1 * t69; -Ifges(7,3) * t159 + t138; t143; t127; t142; t69; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
