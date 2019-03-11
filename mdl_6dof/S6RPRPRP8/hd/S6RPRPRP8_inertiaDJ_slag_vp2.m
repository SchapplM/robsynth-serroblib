% Calculate time derivative of joint inertia matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:32
% EndTime: 2019-03-09 03:24:39
% DurationCPUTime: 3.47s
% Computational Cost: add. (2377->295), mult. (4873->418), div. (0->0), fcn. (4223->6), ass. (0->136)
t208 = Ifges(6,1) + Ifges(7,1);
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t210 = -t90 ^ 2 - t92 ^ 2;
t200 = Ifges(7,4) + Ifges(6,5);
t142 = cos(pkin(9));
t176 = cos(qJ(3));
t89 = sin(pkin(9));
t91 = sin(qJ(3));
t61 = t142 * t91 + t176 * t89;
t209 = t200 * t61;
t207 = -Ifges(6,6) + Ifges(7,6);
t172 = Ifges(7,5) * t90;
t174 = Ifges(6,4) * t90;
t206 = t208 * t92 + t172 - t174;
t139 = qJD(5) * t92;
t57 = t61 * qJD(3);
t154 = t90 * t57;
t111 = t142 * t176;
t60 = t89 * t91 - t111;
t99 = t60 * t139 + t154;
t140 = qJD(5) * t90;
t127 = t60 * t140;
t151 = t92 * t57;
t98 = t127 - t151;
t141 = qJD(3) * t91;
t56 = -qJD(3) * t111 + t141 * t89;
t20 = mrSges(6,2) * t56 + mrSges(6,3) * t99;
t21 = mrSges(7,2) * t99 - mrSges(7,3) * t56;
t149 = t20 + t21;
t18 = -mrSges(6,1) * t56 - mrSges(6,3) * t98;
t19 = t56 * mrSges(7,1) + mrSges(7,2) * t98;
t150 = t18 - t19;
t204 = t149 * t92 - t150 * t90;
t203 = (Ifges(6,4) - Ifges(7,5)) * t99 + t208 * t98 - t200 * t56;
t55 = -pkin(5) * t140 + qJ(6) * t139 + qJD(6) * t90;
t202 = m(7) * t55;
t161 = t57 * t60;
t201 = -Ifges(4,1) + Ifges(4,2);
t198 = -t206 * t60 + t209;
t155 = t90 * mrSges(6,3);
t36 = -mrSges(6,2) * t61 + t155 * t60;
t156 = t90 * mrSges(7,2);
t39 = mrSges(7,3) * t61 + t156 * t60;
t148 = t36 + t39;
t152 = t92 * mrSges(6,3);
t37 = mrSges(6,1) * t61 + t152 * t60;
t153 = t92 * mrSges(7,2);
t38 = -mrSges(7,1) * t61 - t153 * t60;
t147 = t37 - t38;
t80 = pkin(3) * t91 + qJ(2);
t35 = pkin(4) * t61 + pkin(8) * t60 + t80;
t93 = -pkin(1) - pkin(7);
t143 = qJ(4) - t93;
t68 = t143 * t91;
t97 = t143 * t176;
t42 = -t142 * t68 - t89 * t97;
t197 = t35 * t90 + t42 * t92;
t196 = t206 * qJD(5);
t171 = Ifges(7,5) * t92;
t173 = Ifges(6,4) * t92;
t195 = t208 * t90 - t171 + t173;
t15 = t35 * t92 - t42 * t90;
t194 = -t139 * t15 - t140 * t197;
t11 = qJ(6) * t61 + t197;
t12 = -pkin(5) * t61 - t15;
t193 = -t11 * t140 + t12 * t139;
t119 = qJD(3) * t176;
t192 = -mrSges(4,1) * t141 - mrSges(4,2) * t119;
t191 = t139 * t200 + t140 * t207;
t190 = t142 * t57 + t56 * t89;
t104 = -pkin(5) * t92 - qJ(6) * t90;
t79 = -pkin(3) * t142 - pkin(4);
t58 = t104 + t79;
t69 = -mrSges(7,1) * t92 - mrSges(7,3) * t90;
t188 = -m(7) * t58 - t69;
t51 = -qJD(3) * t97 - qJD(4) * t91;
t96 = -qJD(4) * t176 + t141 * t143;
t29 = t142 * t51 + t89 * t96;
t75 = pkin(3) * t119 + qJD(2);
t30 = -pkin(4) * t56 + pkin(8) * t57 + t75;
t4 = -qJD(5) * t197 - t29 * t90 + t30 * t92;
t70 = -mrSges(6,1) * t92 + mrSges(6,2) * t90;
t187 = m(6) * t79 - mrSges(5,1) + t70;
t186 = Ifges(6,6) * t92 + t200 * t90;
t185 = t147 * t90 - t148 * t92;
t184 = 2 * m(5);
t183 = -2 * mrSges(5,3);
t182 = 0.2e1 * qJD(2);
t181 = m(5) * pkin(3);
t28 = -t142 * t96 + t51 * t89;
t41 = t142 * t97 - t68 * t89;
t169 = t28 * t41;
t168 = t41 * t57;
t167 = t56 * t61;
t163 = t56 * t90;
t162 = t56 * t92;
t160 = t60 * t90;
t159 = t60 * t92;
t145 = -Ifges(6,5) * t151 - Ifges(6,3) * t56;
t138 = qJD(6) * t92;
t137 = t60 * t183;
t136 = 0.2e1 * t56 * mrSges(5,3);
t120 = -t56 * mrSges(5,1) - mrSges(5,2) * t57;
t117 = -Ifges(7,4) * t151 - Ifges(7,2) * t56 - Ifges(7,6) * t99;
t3 = t139 * t35 - t140 * t42 + t29 * t92 + t30 * t90;
t1 = -qJ(6) * t56 + qJD(6) * t61 + t3;
t2 = pkin(5) * t56 - t4;
t113 = t1 * t92 + t2 * t90;
t112 = t3 * t92 - t4 * t90;
t110 = t90 * mrSges(6,1) + t92 * mrSges(6,2);
t109 = t90 * mrSges(7,1) - t92 * mrSges(7,3);
t106 = -Ifges(6,2) * t90 + t173;
t105 = Ifges(7,3) * t90 + t171;
t103 = pkin(5) * t90 - qJ(6) * t92;
t102 = t28 * t60 + t168;
t95 = qJD(5) * t104 + t138;
t94 = m(7) * t138 + (m(7) * t104 + t69 + t70) * qJD(5);
t78 = pkin(3) * t89 + pkin(8);
t72 = Ifges(6,2) * t92 + t174;
t71 = -Ifges(7,3) * t92 + t172;
t65 = t106 * qJD(5);
t64 = t105 * qJD(5);
t63 = t110 * qJD(5);
t62 = t109 * qJD(5);
t34 = t110 * t60;
t33 = t109 * t60;
t23 = t61 * Ifges(6,6) - t106 * t60;
t22 = Ifges(7,6) * t61 - t105 * t60;
t17 = -t103 * t60 + t41;
t14 = -mrSges(6,1) * t99 + mrSges(6,2) * t98;
t13 = -mrSges(7,1) * t99 - mrSges(7,3) * t98;
t7 = Ifges(6,4) * t98 + Ifges(6,2) * t99 - t56 * Ifges(6,6);
t6 = Ifges(7,5) * t98 - t56 * Ifges(7,6) - Ifges(7,3) * t99;
t5 = -t103 * t57 + t60 * t95 + t28;
t8 = [t61 * t117 + 0.2e1 * t75 * (mrSges(5,1) * t61 - mrSges(5,2) * t60) - 0.2e1 * t5 * t33 + 0.2e1 * t3 * t36 + 0.2e1 * t4 * t37 + 0.2e1 * t2 * t38 + 0.2e1 * t1 * t39 + 0.2e1 * t41 * t14 + 0.2e1 * t12 * t19 + 0.2e1 * t11 * t21 + 0.2e1 * t17 * t13 + 0.2e1 * t15 * t18 - t203 * t159 + 0.2e1 * Ifges(5,1) * t161 - t198 * t151 + (-0.2e1 * Ifges(4,4) * t176 + t201 * t91) * t119 + (0.2e1 * Ifges(4,4) * t91 + t176 * t201) * t141 + (-t6 + t7) * t160 + t61 * (Ifges(6,6) * t99 + t145) + 0.2e1 * m(7) * (t1 * t11 + t12 * t2 + t17 * t5) + (-0.2e1 * t34 + t137) * t28 + (t91 * mrSges(4,1) + mrSges(4,2) * t176 + mrSges(3,3)) * t182 + 0.2e1 * t80 * t120 + t99 * (-t22 + t23) + 0.2e1 * (-t56 * t60 + t57 * t61) * Ifges(5,4) + (t29 * t42 + t75 * t80 + t169) * t184 + 0.2e1 * t197 * t20 + 0.2e1 * m(6) * (t15 * t4 + t197 * t3 + t169) + ((t200 * t92 + t207 * t90) * t60 + (-(2 * Ifges(5,2)) - Ifges(7,2) - Ifges(6,3)) * t61) * t56 + (t198 + t209) * t127 + (t29 * t61 + t168) * t183 + t42 * t136 + ((m(4) + m(3)) * t182 + 0.2e1 * (mrSges(4,1) * t176 - mrSges(4,2) * t91) * qJD(3)) * qJ(2); (t13 + t14) * t60 + (-t33 - t34 + t137) * t57 + t185 * t56 + m(6) * (t15 * t163 - t162 * t197 + t102) + m(7) * (-t11 * t162 - t12 * t163 + t17 * t57 + t5 * t60) + m(5) * (-t42 * t56 + t102) + (t136 + (-t147 * t92 - t148 * t90) * qJD(5) + m(6) * (t112 + t194) + m(7) * (t113 + t193) + m(5) * t29 + t204) * t61; (t161 - t167) * t184 + 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t167 * t210 + t161); -t92 * t6 / 0.2e1 + t92 * t7 / 0.2e1 + t5 * t69 + t79 * t14 + t17 * t62 + t41 * t63 - Ifges(5,5) * t57 + t58 * t13 + t55 * t33 + Ifges(5,6) * t56 + t203 * t90 / 0.2e1 + (t60 * t72 + t198) * t139 / 0.2e1 + m(7) * (-t17 * t55 + t5 * t58) - t195 * t151 / 0.2e1 - (qJD(5) * t71 + t196) * t159 / 0.2e1 + t191 * t61 / 0.2e1 + t192 * t93 + t193 * mrSges(7,2) + t194 * mrSges(6,3) + (t181 * t89 - mrSges(5,2)) * t29 - t4 * t155 - t23 * t140 / 0.2e1 - Ifges(4,5) * t141 + (-t71 / 0.2e1 + t72 / 0.2e1) * t154 + (t65 / 0.2e1 - t64 / 0.2e1) * t160 - Ifges(4,6) * t119 + t190 * mrSges(5,3) * pkin(3) + (m(6) * ((-t15 * t92 - t197 * t90) * qJD(5) + t112) - t147 * t139 - t148 * t140 + m(7) * ((-t11 * t90 + t12 * t92) * qJD(5) + t113) + t204) * t78 + t3 * t152 + t1 * t153 + t2 * t156 + (t195 * t60 + t22) * t140 / 0.2e1 - (-Ifges(7,6) * t92 + t186) * t56 / 0.2e1 + (-t142 * t181 + t187) * t28; -t190 * t181 + (t62 + t63 - t202) * t60 + (t187 - t188) * t57 + t192 + (mrSges(5,2) + ((m(6) + m(7)) * t78 + mrSges(6,3) + mrSges(7,2)) * t210) * t56; 0.2e1 * t58 * t62 + 0.2e1 * t63 * t79 + (-t64 + t65) * t92 + t196 * t90 + 0.2e1 * t188 * t55 + (t195 * t92 + (t71 - t72) * t90) * qJD(5); t150 * t92 + t149 * t90 - t185 * qJD(5) + m(6) * (t3 * t90 + t4 * t92 + (-t15 * t90 + t197 * t92) * qJD(5)) + m(7) * (t1 * t90 - t2 * t92 + (t11 * t92 + t12 * t90) * qJD(5)) + m(5) * t75 + t120; 0; 0; 0; Ifges(6,6) * t154 - pkin(5) * t19 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t11) + qJD(6) * t39 + qJ(6) * t21 + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) + t186 * t60 * qJD(5) + t117 + t145; (m(7) * t103 + t109 + t110) * t56 + t94 * t61; mrSges(7,2) * t95 + t78 * t94 + t191; t202 + ((-mrSges(6,2) + mrSges(7,3)) * t92 + (-mrSges(6,1) - mrSges(7,1)) * t90) * qJD(5); 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t19; (t139 * t61 - t163) * m(7); (m(7) * t78 + mrSges(7,2)) * t139; m(7) * t140; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
