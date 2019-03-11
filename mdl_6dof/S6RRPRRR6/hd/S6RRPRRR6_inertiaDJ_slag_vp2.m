% Calculate time derivative of joint inertia matrix for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:33
% EndTime: 2019-03-09 13:50:45
% DurationCPUTime: 4.84s
% Computational Cost: add. (6963->392), mult. (14422->576), div. (0->0), fcn. (13493->8), ass. (0->175)
t254 = qJD(2) - qJD(4);
t133 = sin(qJ(6));
t137 = cos(qJ(6));
t113 = -t137 * mrSges(7,1) + mrSges(7,2) * t133;
t268 = mrSges(6,1) - t113;
t131 = t133 ^ 2;
t132 = t137 ^ 2;
t188 = t131 + t132;
t135 = sin(qJ(4));
t139 = cos(qJ(4));
t140 = cos(qJ(2));
t237 = pkin(7) - pkin(8);
t176 = t237 * t140;
t136 = sin(qJ(2));
t260 = t237 * t136;
t267 = t135 * t176 - t139 * t260;
t266 = m(4) * pkin(7) + mrSges(4,2);
t160 = mrSges(7,1) * t133 + mrSges(7,2) * t137;
t105 = t160 * qJD(6);
t218 = mrSges(7,3) * t132;
t219 = mrSges(7,3) * t131;
t253 = qJD(4) + qJD(5);
t134 = sin(qJ(5));
t138 = cos(qJ(5));
t99 = t134 * t135 - t138 * t139;
t68 = t253 * t99;
t102 = t134 * t139 + t135 * t138;
t69 = t253 * t102;
t149 = t99 * t105 - (-mrSges(6,2) + t218 + t219) * t68 - t268 * t69;
t265 = t149 - (mrSges(5,1) * t135 + mrSges(5,2) * t139) * qJD(4);
t74 = t135 * t260 + t139 * t176;
t154 = t135 * t136 + t139 * t140;
t155 = t135 * t140 - t136 * t139;
t156 = t134 * t155 - t138 * t154;
t62 = -t134 * t154 - t138 * t155;
t259 = -t140 * pkin(2) - t136 * qJ(3);
t112 = -pkin(1) + t259;
t95 = t140 * pkin(3) - t112;
t72 = pkin(4) * t154 + t95;
t32 = -pkin(5) * t156 - pkin(10) * t62 + t72;
t147 = pkin(9) * t155 - t267;
t57 = -pkin(9) * t154 + t74;
t37 = t134 * t147 + t138 * t57;
t17 = t133 * t32 + t137 * t37;
t190 = t17 * qJD(6);
t47 = t254 * t74;
t71 = t254 * t154;
t144 = -t71 * pkin(9) + t47;
t209 = t134 * t57;
t46 = t254 * t267;
t70 = t254 * t155;
t261 = -t70 * pkin(9) + qJD(5) * t147 + t46;
t11 = -qJD(5) * t209 + t134 * t144 + t261 * t138;
t29 = qJD(5) * t156 - t134 * t70 + t138 * t71;
t30 = qJD(5) * t62 + t134 * t71 + t138 * t70;
t141 = -pkin(2) - pkin(3);
t187 = qJD(2) * t136;
t186 = qJD(2) * t140;
t189 = qJ(3) * t186 + t136 * qJD(3);
t76 = t141 * t187 + t189;
t48 = pkin(4) * t70 + t76;
t13 = pkin(5) * t30 - pkin(10) * t29 + t48;
t3 = -t11 * t133 + t13 * t137 - t190;
t264 = t3 + t190;
t182 = qJD(6) * t137;
t183 = qJD(6) * t133;
t107 = Ifges(7,4) * t182 - Ifges(7,2) * t183;
t108 = Ifges(7,1) * t182 - Ifges(7,4) * t183;
t263 = t137 * t107 + t133 * t108;
t262 = mrSges(7,3) * t188;
t110 = -qJ(3) * t135 + t139 * t141;
t258 = t188 * t138;
t119 = pkin(4) * t134 + pkin(10);
t257 = t188 * t119;
t16 = -t133 * t37 + t137 * t32;
t256 = -t133 * t16 + t137 * t17;
t252 = t47 * mrSges(5,1) - t46 * mrSges(5,2) + Ifges(5,5) * t71 - Ifges(5,6) * t70;
t106 = Ifges(7,5) * t182 - Ifges(7,6) * t183;
t114 = Ifges(7,5) * t133 + Ifges(7,6) * t137;
t184 = qJD(5) * t138;
t12 = t261 * t134 - t138 * t144 + t57 * t184;
t200 = qJD(6) * t16;
t2 = t11 * t137 + t13 * t133 + t200;
t231 = t137 * t2;
t36 = -t138 * t147 + t209;
t251 = -t11 * mrSges(6,2) + mrSges(7,3) * t231 + Ifges(6,5) * t29 + t36 * t105 - t156 * t106 / 0.2e1 + (-Ifges(6,6) + t114 / 0.2e1) * t30 - t268 * t12;
t250 = 2 * m(5);
t249 = 2 * m(6);
t248 = 2 * m(7);
t247 = -0.2e1 * pkin(1);
t246 = 0.2e1 * t12;
t245 = 0.2e1 * t48;
t244 = 0.2e1 * t112;
t243 = m(6) / 0.2e1;
t241 = -t29 / 0.2e1;
t240 = t29 / 0.2e1;
t234 = pkin(5) * t105;
t233 = t12 * t36;
t232 = t133 * t3;
t109 = -pkin(4) + t110;
t111 = t139 * qJ(3) + t135 * t141;
t67 = t134 * t109 + t138 * t111;
t78 = t139 * qJD(3) + t110 * qJD(4);
t79 = -t135 * qJD(3) - qJD(4) * t111;
t42 = qJD(5) * t67 + t134 * t78 - t138 * t79;
t230 = t36 * t42;
t66 = t109 * t138 - t111 * t134;
t41 = qJD(5) * t66 + t134 * t79 + t138 * t78;
t229 = t41 * mrSges(6,2);
t228 = t42 * t99;
t225 = t69 * t99;
t224 = t78 * mrSges(5,2);
t223 = t79 * mrSges(5,1);
t222 = qJD(6) / 0.2e1;
t221 = Ifges(4,5) - Ifges(3,4);
t206 = t137 * t29;
t220 = Ifges(7,5) * t206 + Ifges(7,3) * t30;
t217 = Ifges(7,4) * t133;
t216 = Ifges(7,4) * t137;
t215 = pkin(4) * qJD(5);
t213 = t133 * t29;
t212 = t133 * t62;
t211 = t134 * mrSges(6,1);
t208 = t134 * t99;
t205 = t137 * t62;
t204 = t138 * mrSges(6,2);
t35 = t42 * t113;
t64 = pkin(5) - t66;
t54 = t64 * t105;
t65 = -pkin(10) + t67;
t199 = qJD(6) * t65;
t120 = -pkin(4) * t138 - pkin(5);
t196 = t120 * t105;
t115 = Ifges(7,2) * t137 + t217;
t193 = t133 * t115;
t192 = t134 * t113;
t116 = Ifges(7,1) * t133 + t216;
t191 = t137 * t116;
t181 = 0.2e1 * t140;
t180 = pkin(4) * t184;
t175 = t62 * t183;
t151 = t175 - t206;
t152 = t182 * t62 + t213;
t8 = mrSges(7,1) * t152 - mrSges(7,2) * t151;
t179 = m(7) * t12 + t8;
t174 = pkin(10) * t188;
t173 = t115 * t222;
t171 = -t182 / 0.2e1;
t170 = -Ifges(7,6) * t133 - (2 * Ifges(6,4));
t169 = t188 * t65;
t168 = t102 * t188;
t164 = t138 * t262;
t163 = t12 * t99 + t36 * t69;
t162 = -t140 * mrSges(4,1) - t136 * mrSges(4,3);
t14 = mrSges(7,1) * t30 + mrSges(7,3) * t151;
t15 = -mrSges(7,2) * t30 - mrSges(7,3) * t152;
t159 = -t133 * t14 + t137 * t15;
t44 = mrSges(7,2) * t156 - mrSges(7,3) * t212;
t45 = -mrSges(7,1) * t156 - mrSges(7,3) * t205;
t158 = -t133 * t45 + t137 * t44;
t97 = t116 * t182;
t150 = -t115 * t183 + t263 + t97;
t148 = -t232 + (-t133 * t17 - t137 * t16) * qJD(6);
t38 = t41 * t219;
t39 = t41 * t218;
t40 = t42 * mrSges(6,1);
t145 = t35 + t38 + t39 - t40 + t54 + 0.2e1 * t133 * t173 - t97 / 0.2e1 + t116 * t171 - t229 - t263;
t143 = m(7) * (-t16 * t182 - t17 * t183 + t231 - t232) - t44 * t183 - t45 * t182 + t159;
t23 = -Ifges(7,6) * t156 + (-Ifges(7,2) * t133 + t216) * t62;
t24 = -Ifges(7,5) * t156 + (Ifges(7,1) * t137 - t217) * t62;
t6 = -Ifges(7,4) * t151 - Ifges(7,2) * t152 + Ifges(7,6) * t30;
t7 = -Ifges(7,1) * t151 - Ifges(7,4) * t152 + Ifges(7,5) * t30;
t142 = t148 * mrSges(7,3) + t193 * t241 + t191 * t240 + t24 * t182 / 0.2e1 + t133 * t7 / 0.2e1 - t107 * t212 / 0.2e1 + t137 * t6 / 0.2e1 + t108 * t205 / 0.2e1 + t62 * t115 * t171 + t251 - (t62 * t116 + t23) * t183 / 0.2e1;
t43 = t160 * t62;
t1 = [t24 * t206 - t23 * t213 + 0.2e1 * t16 * t14 + 0.2e1 * t17 * t15 + 0.2e1 * t36 * t8 + t43 * t246 + 0.2e1 * t2 * t44 + 0.2e1 * t3 * t45 + 0.2e1 * t72 * (mrSges(6,1) * t30 + mrSges(6,2) * t29) + 0.2e1 * t95 * (mrSges(5,1) * t70 + mrSges(5,2) * t71) + 0.2e1 * t154 * Ifges(5,2) * t70 - 0.2e1 * t155 * t71 * Ifges(5,1) + 0.2e1 * t76 * (mrSges(5,1) * t154 - mrSges(5,2) * t155) + (-t267 * t47 + t46 * t74 + t76 * t95) * t250 + (t11 * t37 + t48 * t72 + t233) * t249 + (t16 * t3 + t17 * t2 + t233) * t248 - (mrSges(6,1) * t245 - 0.2e1 * t11 * mrSges(6,3) + ((2 * Ifges(6,2)) + Ifges(7,3)) * t30 + t170 * t29 + t220) * t156 + (mrSges(6,2) * t245 + mrSges(6,3) * t246 + 0.2e1 * Ifges(6,1) * t29 - t133 * t6 + t137 * t7 + (Ifges(7,5) * t137 + t170) * t30 + (t114 * t156 - t133 * t24 - t137 * t23) * qJD(6)) * t62 + (m(4) * t244 + 0.2e1 * t162) * (pkin(2) * t187 - t189) + 0.2e1 * (t29 * t36 - t30 * t37) * mrSges(6,3) + 0.2e1 * (-t154 * t71 + t155 * t70) * Ifges(5,4) + 0.2e1 * (-t154 * t46 + t155 * t47 + t267 * t71 - t70 * t74) * mrSges(5,3) + ((mrSges(3,2) * t247 - 0.2e1 * t112 * mrSges(4,3) - t181 * t221) * t140 + (mrSges(3,1) * t247 + mrSges(4,1) * t244 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t181 + 0.2e1 * t221 * t136) * t136) * qJD(2); -t251 + (-t110 * t71 - t111 * t70 - t154 * t78 + t155 * t79) * mrSges(5,3) + (t156 * t41 - t66 * t29 - t67 * t30 + t42 * t62) * mrSges(6,3) + m(5) * (t110 * t47 + t111 * t46 - t267 * t79 + t74 * t78) + (t23 * t222 + t115 * t240 - t7 / 0.2e1 - t44 * t199 + (t107 / 0.2e1 + t116 * t222) * t62 + t264 * mrSges(7,3) + (-m(7) * t16 - t45) * t41 + (-m(7) * t264 - t14) * t65) * t133 + t42 * t43 + t64 * t8 + m(6) * (t11 * t67 - t12 * t66 + t37 * t41 + t230) + m(7) * (t12 * t64 + t230) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t140 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t136 + (m(4) * t259 - t140 * mrSges(3,1) + t136 * mrSges(3,2) + t162) * pkin(7)) * qJD(2) + (-qJD(6) * t24 / 0.2e1 + t116 * t241 + m(7) * (-t16 * t199 + t17 * t41 + t2 * t65) + t65 * t15 + t41 * t44 - t6 / 0.2e1 + mrSges(7,3) * t200 - t45 * t199 + (-t108 / 0.2e1 + t173) * t62) * t137 - t252 + t266 * qJD(3) * t140; -0.2e1 * t223 + 0.2e1 * t224 + 0.2e1 * t229 - 0.2e1 * t54 - 0.2e1 * t35 - 0.2e1 * t38 - 0.2e1 * t39 + 0.2e1 * t40 + (t191 - t193) * qJD(6) + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) + (t169 * t41 + t42 * t64) * t248 + (t41 * t67 - t42 * t66) * t249 + (t110 * t79 + t111 * t78) * t250 + t263; t69 * t43 + t99 * t8 - t158 * t68 + t266 * t186 + ((-t133 * t44 - t137 * t45) * qJD(6) + t159) * t102 + m(7) * (-t256 * t68 + (t148 + t231) * t102 + t163) + m(6) * (t102 * t11 - t37 * t68 + t163) + m(5) * (t135 * t46 + t139 * t47 + (t135 * t267 + t139 * t74) * qJD(4)) + (-t102 * t30 - t156 * t68 + t29 * t99 + t62 * t69) * mrSges(6,3) + (-t135 * t70 - t139 * t71 + (-t135 * t155 - t139 * t154) * qJD(4)) * mrSges(5,3); m(7) * (t168 * t41 - t169 * t68 + t64 * t69 + t228) + m(6) * (t102 * t41 - t66 * t69 - t67 * t68 + t228) + m(5) * (t135 * t78 + t139 * t79 + (-t110 * t135 + t111 * t139) * qJD(4)) - t265; 0.2e1 * m(6) * (-t102 * t68 + t225) + 0.2e1 * m(7) * (-t168 * t68 + t225); t142 + (m(6) * (t11 * t134 - t12 * t138) + (-t134 * t30 - t138 * t29) * mrSges(6,3) + ((m(6) * t37 + m(7) * t256 + mrSges(6,3) * t156 + t158) * t138 + (t62 * mrSges(6,3) + t43 + (m(6) + m(7)) * t36) * t134) * qJD(5)) * pkin(4) + t179 * t120 + t143 * t119 + t252; t145 + m(7) * (t120 * t42 + t257 * t41) + (m(6) * (t134 * t41 - t138 * t42) + (t268 * t134 + (mrSges(6,2) - t262) * t138 + m(6) * (-t134 * t66 + t138 * t67) + m(7) * (t134 * t64 + t258 * t65)) * qJD(5)) * pkin(4) - t224 + t223 - t196; m(7) * (t120 * t69 - t257 * t68) + 0.2e1 * ((-t134 * t68 - t138 * t69) * t243 + ((t102 * t138 + t208) * t243 + m(7) * (t258 * t102 + t208) / 0.2e1) * qJD(5)) * pkin(4) + t265; 0.2e1 * t196 + (-0.2e1 * t204 - 0.2e1 * t211 + 0.2e1 * t192 + (t258 * t119 + t120 * t134) * t248 + 0.2e1 * t164) * t215 + t150; -pkin(5) * t179 + pkin(10) * t143 + t142; t145 + m(7) * (-pkin(5) * t42 + t174 * t41) + t234; m(7) * (-pkin(5) * t69 - t174 * t68) + t149; (t120 - pkin(5)) * t105 + (t192 + m(7) * (-pkin(5) * t134 + t258 * pkin(10)) - t204 - t211 + t164) * t215 + t150; t150 - 0.2e1 * t234; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t175 - Ifges(7,6) * t152 + t220; (-t137 * t41 + t183 * t65) * mrSges(7,2) + (-t133 * t41 - t182 * t65) * mrSges(7,1) - t106; (t102 * t183 + t137 * t68) * mrSges(7,2) + (-t102 * t182 + t133 * t68) * mrSges(7,1); (t119 * t183 - t137 * t180) * mrSges(7,2) + (-t119 * t182 - t133 * t180) * mrSges(7,1) + t106; pkin(10) * qJD(6) * t113 + t106; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
