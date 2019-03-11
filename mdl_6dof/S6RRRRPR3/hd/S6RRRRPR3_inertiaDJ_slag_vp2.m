% Calculate time derivative of joint inertia matrix for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:02:03
% EndTime: 2019-03-09 22:02:13
% DurationCPUTime: 4.52s
% Computational Cost: add. (7935->348), mult. (17216->514), div. (0->0), fcn. (16650->8), ass. (0->179)
t115 = sin(qJ(6));
t117 = cos(qJ(6));
t199 = t115 ^ 2 + t117 ^ 2;
t172 = t199 * mrSges(7,3);
t269 = mrSges(5,3) + mrSges(6,1);
t196 = qJD(6) * t117;
t116 = sin(qJ(4));
t235 = sin(qJ(3));
t236 = sin(qJ(2));
t238 = cos(qJ(3));
t239 = cos(qJ(2));
t141 = t235 * t236 - t238 * t239;
t137 = t141 * qJD(3);
t124 = -qJD(2) * t141 - t137;
t257 = qJD(2) + qJD(3);
t92 = t235 * t239 + t236 * t238;
t125 = t257 * t92;
t138 = t116 * t141;
t237 = cos(qJ(4));
t175 = qJD(4) * t237;
t41 = -qJD(4) * t138 + t116 * t124 + t125 * t237 + t175 * t92;
t135 = t237 * t141;
t73 = t116 * t92 + t135;
t152 = t115 * t41 + t73 * t196;
t215 = Ifges(7,4) * t117;
t101 = -Ifges(7,2) * t115 + t215;
t216 = Ifges(7,4) * t115;
t102 = Ifges(7,1) * t117 - t216;
t268 = t101 * t117 + t102 * t115;
t197 = qJD(6) * t115;
t262 = t117 * t41 - t73 * t197;
t74 = t237 * t92 - t138;
t110 = -pkin(2) * t239 - pkin(1);
t80 = pkin(3) * t141 + t110;
t132 = -t74 * qJ(5) + t80;
t240 = pkin(4) + pkin(10);
t23 = t240 * t73 + t132;
t264 = -pkin(8) - pkin(7);
t148 = t264 * t236;
t143 = t238 * t148;
t149 = t264 * t239;
t78 = t235 * t149 + t143;
t130 = -t92 * pkin(9) + t78;
t61 = t237 * t130;
t93 = t235 * t148;
t79 = -t238 * t149 + t93;
t62 = -pkin(9) * t141 + t79;
t46 = t116 * t62 - t61;
t24 = pkin(5) * t74 + t46;
t18 = -t115 * t23 + t117 * t24;
t19 = t115 * t24 + t117 * t23;
t209 = t115 * t73;
t49 = mrSges(7,1) * t74 - mrSges(7,3) * t209;
t205 = t117 * t73;
t50 = -mrSges(7,2) * t74 + mrSges(7,3) * t205;
t267 = m(7) * (t115 * t19 + t117 * t18) + t117 * t49 + t115 * t50 + t74 * t269;
t198 = qJD(4) * t116;
t40 = qJD(4) * t135 + t116 * t125 - t124 * t237 + t198 * t92;
t266 = -0.2e1 * t40;
t265 = 0.2e1 * t73;
t263 = Ifges(3,1) - Ifges(3,2);
t145 = qJD(2) * t149;
t52 = qJD(2) * t143 + qJD(3) * t78 + t145 * t235;
t119 = -pkin(9) * t125 + t52;
t53 = -qJD(2) * t93 - qJD(3) * t79 + t145 * t238;
t120 = -pkin(9) * t124 + t53;
t126 = t116 * t130;
t16 = qJD(4) * t126 + t116 * t119 - t237 * t120 + t175 * t62;
t10 = -t40 * pkin(5) + t16;
t203 = qJD(6) * t19;
t177 = qJD(2) * t236;
t111 = pkin(2) * t177;
t64 = pkin(3) * t125 + t111;
t12 = t41 * pkin(4) + t40 * qJ(5) - t74 * qJD(5) + t64;
t8 = t41 * pkin(10) + t12;
t2 = t10 * t117 - t115 * t8 - t203;
t1 = qJD(6) * t18 + t10 * t115 + t117 * t8;
t230 = t1 * t115;
t261 = t2 * t117 + t230;
t20 = mrSges(7,2) * t40 + mrSges(7,3) * t262;
t21 = -mrSges(7,1) * t40 - mrSges(7,3) * t152;
t131 = m(7) * ((-t115 * t18 + t117 * t19) * qJD(6) + t261) + t117 * t21 + t115 * t20;
t247 = t50 * t196 - t49 * t197 + t131;
t228 = t40 * mrSges(6,1);
t255 = m(6) * t16 - t228;
t190 = t238 * pkin(2);
t165 = t190 + pkin(3);
t150 = t237 * t165;
t168 = qJD(3) * t190;
t176 = qJD(3) * t235;
t68 = -(qJD(4) * t235 + t176) * pkin(2) * t116 + qJD(4) * t150 + t237 * t168;
t188 = pkin(3) * t198;
t107 = mrSges(6,2) * t188;
t166 = pkin(3) * t175;
t104 = t166 + qJD(5);
t213 = t104 * mrSges(6,3);
t231 = pkin(3) * t116;
t108 = qJ(5) + t231;
t162 = mrSges(7,1) * t117 - mrSges(7,2) * t115;
t95 = t162 * qJD(6);
t81 = t108 * t95;
t100 = mrSges(7,1) * t115 + mrSges(7,2) * t117;
t84 = t104 * t100;
t254 = t107 + t213 + t81 + t84;
t15 = -qJD(4) * t61 - t116 * t120 - t237 * t119 + t198 * t62;
t159 = Ifges(7,5) * t115 + Ifges(7,6) * t117;
t180 = t205 / 0.2e1;
t47 = t237 * t62 + t126;
t25 = -t73 * pkin(5) + t47;
t160 = Ifges(7,2) * t117 + t216;
t31 = Ifges(7,6) * t74 + t160 * t73;
t161 = Ifges(7,1) * t115 + t215;
t32 = Ifges(7,5) * t74 + t161 * t73;
t5 = Ifges(7,4) * t152 + Ifges(7,2) * t262 - t40 * Ifges(7,6);
t6 = Ifges(7,1) * t152 + Ifges(7,4) * t262 - t40 * Ifges(7,5);
t9 = -pkin(5) * t41 - t15;
t96 = t160 * qJD(6);
t97 = t161 * qJD(6);
t251 = t18 * mrSges(7,3) * t197 + t25 * t95 - t31 * t196 / 0.2e1 - t115 * t5 / 0.2e1 + t117 * t6 / 0.2e1 - t97 * t209 / 0.2e1 - t96 * t180 + t9 * t100 + (-Ifges(7,5) * t117 / 0.2e1 + Ifges(7,6) * t115 / 0.2e1 + Ifges(6,4) - Ifges(5,5)) * t40 - (t101 * t73 + t32) * t197 / 0.2e1 + (mrSges(6,2) - mrSges(5,1)) * t16 + (t102 * t180 - t74 * t159 / 0.2e1) * qJD(6) + (mrSges(5,2) - mrSges(6,3)) * t15 + (Ifges(6,5) - Ifges(5,6) + t268 / 0.2e1) * t41;
t221 = t68 * mrSges(5,2);
t67 = -qJD(5) - t68;
t223 = t67 * mrSges(6,3);
t58 = t67 * t100;
t164 = t237 * t235;
t69 = t188 + (qJD(3) + qJD(4)) * pkin(2) * (t116 * t238 + t164);
t65 = t69 * mrSges(6,2);
t66 = t69 * mrSges(5,1);
t88 = pkin(2) * t164 + t116 * t165;
t83 = qJ(5) + t88;
t77 = t83 * t95;
t250 = t65 + t77 - t221 - t223 - t58 - t66;
t11 = -mrSges(7,1) * t262 + mrSges(7,2) * t152;
t249 = -m(6) * t15 + m(7) * t9 - t41 * mrSges(6,1) + t11;
t48 = t162 * t73;
t248 = m(6) * t47 + m(7) * t25 - t73 * mrSges(6,1) - t48;
t245 = 2 * m(5);
t244 = 0.2e1 * m(6);
t243 = 0.2e1 * m(7);
t242 = m(5) * pkin(3);
t233 = Ifges(4,4) * t92;
t227 = t40 * mrSges(5,3);
t225 = t41 * mrSges(5,3);
t224 = t46 * t69;
t222 = t67 * t83;
t218 = t73 * mrSges(5,3);
t217 = t92 * mrSges(4,3);
t214 = pkin(3) * qJD(4);
t200 = t104 * t108;
t191 = pkin(2) * t235;
t189 = t237 * pkin(3);
t187 = mrSges(5,2) * t237;
t179 = m(7) * t199;
t178 = qJD(2) * t239;
t173 = t115 * t96 - t117 * t97;
t171 = t199 * t69;
t167 = pkin(2) * t176;
t109 = -t189 - pkin(4);
t163 = -t15 * t47 + t16 * t46;
t158 = t104 * t83 - t108 * t67;
t156 = -0.2e1 * t172;
t155 = -qJ(5) * t67 + qJD(5) * t83;
t154 = qJ(5) * t104 + qJD(5) * t108;
t153 = t199 * t188;
t146 = t152 * Ifges(7,5) + t262 * Ifges(7,6) - Ifges(7,3) * t40;
t87 = -t116 * t191 + t150;
t85 = -pkin(4) - t87;
t140 = -t268 * qJD(6) + t173;
t139 = Ifges(4,4) * t141;
t114 = qJD(5) * mrSges(6,3);
t91 = qJ(5) * t95;
t98 = qJD(5) * t100;
t134 = t114 + t91 + t98 + t140;
t123 = mrSges(4,3) * t125;
t122 = mrSges(4,3) * t124;
t121 = -t52 * mrSges(4,2) + t53 * mrSges(4,1) - Ifges(4,6) * t125 + Ifges(4,5) * t124 + (-t19 * t196 - t261) * mrSges(7,3) + t251;
t106 = -pkin(10) + t109;
t82 = -pkin(10) + t85;
t45 = t73 * pkin(4) + t132;
t3 = [(Ifges(6,2) + Ifges(5,1)) * t74 * t266 + (Ifges(6,3) + Ifges(5,2)) * t41 * t265 - t40 * (Ifges(7,3) * t74 + t159 * t73) + t269 * (t15 * t265 + 0.2e1 * t16 * t74 + t46 * t266 - 0.2e1 * t41 * t47) - 0.2e1 * pkin(1) * (mrSges(3,1) * t236 + mrSges(3,2) * t239) * qJD(2) - 0.2e1 * t78 * t122 - 0.2e1 * t79 * t123 - 0.2e1 * t52 * t141 * mrSges(4,3) + 0.2e1 * (Ifges(6,6) + Ifges(5,4)) * (t40 * t73 - t41 * t74) + t124 * (Ifges(4,1) * t92 - t139) + (t64 * t80 + t163) * t245 + (-0.2e1 * Ifges(3,4) * t236 + t239 * t263) * t177 + (0.2e1 * Ifges(3,4) * t239 + t236 * t263) * t178 + t262 * t31 + t152 * t32 + (t1 * t19 + t18 * t2 + t25 * t9) * t243 + (t12 * t45 + t163) * t244 + t5 * t205 + t6 * t209 + (t92 * (-Ifges(4,1) * t141 - t233) + 0.2e1 * t110 * (mrSges(4,1) * t92 - mrSges(4,2) * t141) - t141 * (-Ifges(4,2) * t92 - t139)) * t257 + t74 * t146 + 0.2e1 * (mrSges(4,1) * t141 + t92 * mrSges(4,2)) * t111 + 0.2e1 * m(4) * (t110 * t111 + t79 * t52 + t78 * t53) - 0.2e1 * t53 * t217 + 0.2e1 * t19 * t20 + 0.2e1 * t18 * t21 + 0.2e1 * t25 * t11 - t125 * (-Ifges(4,2) * t141 + t233) + 0.2e1 * t45 * (-mrSges(6,2) * t41 + mrSges(6,3) * t40) - 0.2e1 * t9 * t48 + 0.2e1 * t2 * t49 + 0.2e1 * t1 * t50 + 0.2e1 * t64 * (mrSges(5,1) * t73 + mrSges(5,2) * t74) + 0.2e1 * t12 * (-mrSges(6,2) * t73 - mrSges(6,3) * t74) + 0.2e1 * t80 * (mrSges(5,1) * t41 - mrSges(5,2) * t40); m(4) * (t238 * t53 + t235 * t52 + (-t235 * t78 + t238 * t79) * qJD(3)) * pkin(2) + t121 + t87 * t227 + t167 * t217 - t123 * t191 + Ifges(3,5) * t178 - Ifges(3,6) * t177 - t68 * t218 + m(6) * (t16 * t85 + t224) + m(5) * (-t15 * t88 - t16 * t87 + t47 * t68 + t224) - t88 * t225 - t85 * t228 + t249 * t83 - t248 * t67 + (-mrSges(4,3) * t137 - t122) * t190 + (-mrSges(3,1) * t178 + mrSges(3,2) * t177) * pkin(7) + t247 * t82 + t267 * t69; -0.2e1 * t221 - 0.2e1 * t223 - 0.2e1 * t58 + 0.2e1 * t65 - 0.2e1 * t66 + 0.2e1 * t77 + 0.2e1 * (-mrSges(4,1) * t235 - mrSges(4,2) * t238) * qJD(3) * pkin(2) + t69 * t156 + (t171 * t82 - t222) * t243 + (t69 * t85 - t222) * t244 + (t68 * t88 - t69 * t87) * t245 + t140; (-t237 * t16 - t116 * t15 + (t116 * t46 + t237 * t47) * qJD(4)) * t242 + t121 - t225 * t231 + t189 * t227 - t166 * t218 + t255 * t109 + t249 * t108 + t248 * t104 + t247 * t106 + (m(6) * t46 + t267) * t188; t250 + (-t237 * t69 + t116 * t68 + (-t116 * t87 + t237 * t88) * qJD(4)) * t242 + t173 - mrSges(5,2) * t166 - mrSges(4,1) * t167 - mrSges(4,2) * t168 + m(7) * (t106 * t171 + t153 * t82 + t158) - mrSges(5,1) * t188 + m(6) * (t109 * t69 + t188 * t85 + t158) - t101 * t196 - t102 * t197 + t254 + (-t188 - t69) * t172; 0.2e1 * t213 + 0.2e1 * t107 + 0.2e1 * t81 + 0.2e1 * t84 + (-0.2e1 * t187 + (-0.2e1 * mrSges(5,1) + t156) * t116) * t214 + (t106 * t153 + t200) * t243 + (t109 * t188 + t200) * t244 + t140; -t247 * t240 + (pkin(4) * t40 - qJ(5) * t41 - qJD(5) * t73) * mrSges(6,1) + m(6) * (-pkin(4) * t16 - qJ(5) * t15 + qJD(5) * t47) + m(7) * (qJ(5) * t9 + qJD(5) * t25) + qJ(5) * t11 + (-t230 + (-t2 - t203) * t117) * mrSges(7,3) - qJD(5) * t48 + t251; -t69 * t172 + m(7) * (-t171 * t240 + t155) + m(6) * (-pkin(4) * t69 + t155) + t134 + t250; (-t187 + (-mrSges(5,1) - t172) * t116) * t214 + m(7) * (-t153 * t240 + t154) + m(6) * (-pkin(4) * t188 + t154) + t134 + t254; 0.2e1 * t114 + 0.2e1 * t91 + 0.2e1 * t98 + 0.2e1 * (m(6) + m(7)) * qJD(5) * qJ(5) + t140; (-t115 * t49 + t117 * t50) * qJD(6) + t131 + t255; 0.2e1 * (t179 / 0.2e1 + m(6) / 0.2e1) * t69; (t179 + m(6)) * t188; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t146; t162 * t69 + ((-mrSges(7,2) * t82 - Ifges(7,6)) * t117 + (-mrSges(7,1) * t82 - Ifges(7,5)) * t115) * qJD(6); t162 * t188 + (-t100 * t106 - t159) * qJD(6); ((mrSges(7,2) * t240 - Ifges(7,6)) * t117 + (mrSges(7,1) * t240 - Ifges(7,5)) * t115) * qJD(6); -t100 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
