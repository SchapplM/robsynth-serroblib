% Calculate time derivative of joint inertia matrix for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:11
% EndTime: 2019-12-31 22:08:25
% DurationCPUTime: 5.27s
% Computational Cost: add. (3778->525), mult. (10458->743), div. (0->0), fcn. (9218->8), ass. (0->217)
t188 = cos(qJ(4));
t271 = Ifges(5,6) + Ifges(6,6);
t274 = t271 * t188;
t185 = sin(qJ(4));
t273 = (Ifges(5,5) + Ifges(6,5)) * t185;
t189 = cos(qJ(3));
t228 = qJD(3) * t189;
t209 = t188 * t228;
t186 = sin(qJ(3));
t226 = qJD(4) * t186;
t212 = t185 * t226;
t192 = t209 - t212;
t210 = t185 * t228;
t225 = qJD(4) * t188;
t191 = t186 * t225 + t210;
t184 = cos(pkin(5));
t183 = sin(pkin(5));
t187 = sin(qJ(2));
t242 = t183 * t187;
t121 = t184 * t186 + t189 * t242;
t190 = cos(qJ(2));
t241 = t183 * t190;
t193 = -t121 * t188 + t185 * t241;
t230 = qJD(2) * t187;
t214 = t183 * t230;
t120 = -t184 * t189 + t186 * t242;
t231 = qJD(2) * t183;
t213 = t190 * t231;
t84 = -t120 * qJD(3) + t189 * t213;
t43 = t193 * qJD(4) - t84 * t185 + t188 * t214;
t85 = -t121 * t185 - t188 * t241;
t44 = t85 * qJD(4) + t185 * t214 + t84 * t188;
t83 = t121 * qJD(3) + t186 * t213;
t5 = Ifges(6,5) * t44 + Ifges(6,6) * t43 + Ifges(6,3) * t83;
t6 = Ifges(5,5) * t44 + Ifges(5,6) * t43 + Ifges(5,3) * t83;
t270 = t5 + t6;
t266 = -m(5) * pkin(9) - mrSges(5,3);
t14 = -mrSges(5,1) * t43 + mrSges(5,2) * t44;
t125 = t184 * t187 * pkin(1) + pkin(7) * t241;
t107 = pkin(8) * t184 + t125;
t108 = (-pkin(2) * t190 - pkin(8) * t187 - pkin(1)) * t183;
t115 = (pkin(2) * t187 - pkin(8) * t190) * t231;
t169 = pkin(7) * t242;
t256 = pkin(1) * t190;
t124 = t184 * t256 - t169;
t116 = t124 * qJD(2);
t229 = qJD(3) * t186;
t26 = -t107 * t228 - t108 * t229 + t115 * t189 - t186 * t116;
t24 = -pkin(3) * t214 - t26;
t265 = -m(5) * t24 - t14;
t148 = -pkin(3) * t189 - pkin(9) * t186 - pkin(2);
t238 = t188 * t189;
t173 = pkin(8) * t238;
t105 = t185 * t148 + t173;
t264 = 0.2e1 * m(5);
t263 = 2 * m(6);
t262 = 0.2e1 * pkin(8);
t261 = -2 * mrSges(3,3);
t260 = -2 * mrSges(6,3);
t258 = m(6) * pkin(4);
t255 = pkin(8) * t185;
t253 = -qJ(5) - pkin(9);
t106 = t169 + (-pkin(2) - t256) * t184;
t55 = t120 * pkin(3) - t121 * pkin(9) + t106;
t65 = t189 * t107 + t186 * t108;
t57 = -pkin(9) * t241 + t65;
t17 = t185 * t55 + t188 * t57;
t252 = Ifges(4,4) * t186;
t251 = Ifges(4,4) * t189;
t250 = Ifges(5,4) * t185;
t249 = Ifges(5,4) * t188;
t248 = Ifges(6,4) * t185;
t247 = Ifges(6,4) * t188;
t246 = t116 * mrSges(3,2);
t117 = t125 * qJD(2);
t245 = t117 * mrSges(3,1);
t244 = t117 * mrSges(4,1);
t243 = t117 * mrSges(4,2);
t240 = t185 * t186;
t239 = t186 * t188;
t194 = -Ifges(6,2) * t185 + t247;
t111 = -Ifges(6,6) * t189 + t194 * t186;
t195 = -Ifges(5,2) * t185 + t249;
t112 = -Ifges(5,6) * t189 + t195 * t186;
t237 = -t111 - t112;
t196 = Ifges(6,1) * t188 - t248;
t113 = -Ifges(6,5) * t189 + t196 * t186;
t197 = Ifges(5,1) * t188 - t250;
t114 = -Ifges(5,5) * t189 + t197 * t186;
t236 = t113 + t114;
t146 = (pkin(3) * t186 - pkin(9) * t189) * qJD(3);
t235 = t185 * t146 + t148 * t225;
t234 = t188 * t146 + t229 * t255;
t233 = Ifges(6,5) * t209 + Ifges(6,3) * t229;
t232 = Ifges(5,5) * t209 + Ifges(5,3) * t229;
t227 = qJD(4) * t185;
t131 = mrSges(6,1) * t227 + mrSges(6,2) * t225;
t224 = qJD(5) * t188;
t7 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t83;
t8 = Ifges(5,4) * t44 + Ifges(5,2) * t43 + Ifges(5,6) * t83;
t223 = -t7 / 0.2e1 - t8 / 0.2e1;
t10 = Ifges(5,1) * t44 + Ifges(5,4) * t43 + Ifges(5,5) * t83;
t9 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t83;
t221 = t9 / 0.2e1 + t10 / 0.2e1;
t220 = Ifges(4,5) * t84 - Ifges(4,6) * t83 + Ifges(4,3) * t214;
t219 = Ifges(4,6) * t241;
t32 = -Ifges(6,4) * t193 + Ifges(6,2) * t85 + Ifges(6,6) * t120;
t33 = -Ifges(5,4) * t193 + Ifges(5,2) * t85 + Ifges(5,6) * t120;
t218 = -t32 / 0.2e1 - t33 / 0.2e1;
t34 = -Ifges(6,1) * t193 + Ifges(6,4) * t85 + Ifges(6,5) * t120;
t35 = -Ifges(5,1) * t193 + Ifges(5,4) * t85 + Ifges(5,5) * t120;
t217 = t34 / 0.2e1 + t35 / 0.2e1;
t155 = Ifges(6,2) * t188 + t248;
t72 = -t155 * t226 + (Ifges(6,6) * t186 + t194 * t189) * qJD(3);
t156 = Ifges(5,2) * t188 + t250;
t73 = -t156 * t226 + (Ifges(5,6) * t186 + t195 * t189) * qJD(3);
t216 = t72 / 0.2e1 + t73 / 0.2e1;
t158 = Ifges(6,1) * t185 + t247;
t74 = -t158 * t226 + (Ifges(6,5) * t186 + t196 * t189) * qJD(3);
t159 = Ifges(5,1) * t185 + t249;
t75 = -t159 * t226 + (Ifges(5,5) * t186 + t197 * t189) * qJD(3);
t215 = t74 / 0.2e1 + t75 / 0.2e1;
t208 = t111 / 0.2e1 + t112 / 0.2e1;
t207 = t113 / 0.2e1 + t114 / 0.2e1;
t180 = Ifges(6,5) * t225;
t181 = Ifges(5,5) * t225;
t206 = t180 / 0.2e1 + t181 / 0.2e1 - t271 * t227 / 0.2e1;
t136 = t194 * qJD(4);
t137 = t195 * qJD(4);
t205 = t136 / 0.2e1 + t137 / 0.2e1;
t139 = t196 * qJD(4);
t140 = t197 * qJD(4);
t204 = t139 / 0.2e1 + t140 / 0.2e1;
t203 = t273 / 0.2e1 + t274 / 0.2e1;
t202 = -t156 / 0.2e1 - t155 / 0.2e1;
t201 = t158 / 0.2e1 + t159 / 0.2e1;
t13 = -t43 * mrSges(6,1) + t44 * mrSges(6,2);
t16 = -t185 * t57 + t188 * t55;
t200 = qJD(4) * t253;
t64 = -t186 * t107 + t108 * t189;
t199 = t214 / 0.2e1;
t56 = pkin(3) * t241 - t64;
t198 = mrSges(5,1) * t185 + mrSges(5,2) * t188;
t25 = -t107 * t229 + t108 * t228 + t186 * t115 + t189 * t116;
t23 = pkin(9) * t214 + t25;
t36 = t83 * pkin(3) - t84 * pkin(9) + t117;
t3 = t185 * t36 + t188 * t23 + t55 * t225 - t57 * t227;
t76 = t191 * mrSges(6,1) + t192 * mrSges(6,2);
t4 = -t17 * qJD(4) - t185 * t23 + t188 * t36;
t182 = Ifges(4,5) * t228;
t175 = -pkin(4) * t188 - pkin(3);
t162 = Ifges(3,5) * t213;
t160 = Ifges(4,1) * t186 + t251;
t157 = Ifges(4,2) * t189 + t252;
t152 = t253 * t188;
t151 = -mrSges(5,1) * t188 + mrSges(5,2) * t185;
t150 = -mrSges(6,1) * t188 + mrSges(6,2) * t185;
t149 = t253 * t185;
t147 = (pkin(4) * t185 + pkin(8)) * t186;
t145 = -mrSges(5,1) * t189 - mrSges(5,3) * t239;
t144 = -mrSges(6,1) * t189 - mrSges(6,3) * t239;
t143 = mrSges(5,2) * t189 - mrSges(5,3) * t240;
t142 = mrSges(6,2) * t189 - mrSges(6,3) * t240;
t141 = (Ifges(4,1) * t189 - t252) * qJD(3);
t138 = (-Ifges(4,2) * t186 + t251) * qJD(3);
t133 = (mrSges(4,1) * t186 + mrSges(4,2) * t189) * qJD(3);
t132 = t198 * qJD(4);
t130 = t188 * t148;
t127 = t198 * t186;
t126 = (mrSges(6,1) * t185 + mrSges(6,2) * t188) * t186;
t119 = -qJD(5) * t185 + t188 * t200;
t118 = t185 * t200 + t224;
t110 = -Ifges(5,3) * t189 + (Ifges(5,5) * t188 - Ifges(5,6) * t185) * t186;
t109 = -Ifges(6,3) * t189 + (Ifges(6,5) * t188 - Ifges(6,6) * t185) * t186;
t104 = -t189 * t255 + t130;
t100 = pkin(4) * t191 + pkin(8) * t228;
t96 = -mrSges(5,2) * t229 - mrSges(5,3) * t191;
t95 = -mrSges(6,2) * t229 - mrSges(6,3) * t191;
t94 = mrSges(5,1) * t229 - mrSges(5,3) * t192;
t93 = mrSges(6,1) * t229 - mrSges(6,3) * t192;
t89 = -mrSges(4,1) * t241 - mrSges(4,3) * t121;
t88 = mrSges(4,2) * t241 - mrSges(4,3) * t120;
t87 = -qJ(5) * t240 + t105;
t78 = -qJ(5) * t239 + t130 + (-pkin(4) - t255) * t189;
t77 = mrSges(5,1) * t191 + mrSges(5,2) * t192;
t71 = -Ifges(5,5) * t212 - Ifges(5,6) * t191 + t232;
t70 = -Ifges(6,5) * t212 - Ifges(6,6) * t191 + t233;
t69 = mrSges(4,1) * t214 - mrSges(4,3) * t84;
t68 = -mrSges(4,2) * t214 - mrSges(4,3) * t83;
t67 = Ifges(4,1) * t121 - Ifges(4,4) * t120 - Ifges(4,5) * t241;
t66 = Ifges(4,4) * t121 - Ifges(4,2) * t120 - t219;
t63 = -qJD(4) * t105 + t234;
t62 = (-t188 * t229 - t189 * t227) * pkin(8) + t235;
t61 = mrSges(5,1) * t120 + mrSges(5,3) * t193;
t60 = mrSges(6,1) * t120 + mrSges(6,3) * t193;
t59 = -mrSges(5,2) * t120 + mrSges(5,3) * t85;
t58 = -mrSges(6,2) * t120 + mrSges(6,3) * t85;
t50 = -mrSges(5,1) * t85 - mrSges(5,2) * t193;
t49 = -mrSges(6,1) * t85 - mrSges(6,2) * t193;
t48 = mrSges(4,1) * t83 + mrSges(4,2) * t84;
t47 = Ifges(4,1) * t84 - Ifges(4,4) * t83 + Ifges(4,5) * t214;
t46 = Ifges(4,4) * t84 - Ifges(4,2) * t83 + Ifges(4,6) * t214;
t45 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t239 + (-qJD(5) * t186 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t189) * t185 + t235;
t37 = -t186 * t224 + (pkin(4) * t186 - qJ(5) * t238) * qJD(3) + (-t173 + (qJ(5) * t186 - t148) * t185) * qJD(4) + t234;
t31 = -Ifges(5,5) * t193 + Ifges(5,6) * t85 + Ifges(5,3) * t120;
t30 = -Ifges(6,5) * t193 + Ifges(6,6) * t85 + Ifges(6,3) * t120;
t27 = -pkin(4) * t85 + t56;
t21 = mrSges(5,1) * t83 - mrSges(5,3) * t44;
t20 = mrSges(6,1) * t83 - mrSges(6,3) * t44;
t19 = -mrSges(5,2) * t83 + mrSges(5,3) * t43;
t18 = -mrSges(6,2) * t83 + mrSges(6,3) * t43;
t15 = qJ(5) * t85 + t17;
t12 = pkin(4) * t120 + qJ(5) * t193 + t16;
t11 = -pkin(4) * t43 + t24;
t2 = qJ(5) * t43 + qJD(5) * t85 + t3;
t1 = pkin(4) * t83 - qJ(5) * t44 + qJD(5) * t193 + t4;
t22 = [(-t190 * t220 + 0.2e1 * (t116 * t190 + t117 * t187) * mrSges(3,3) + ((t124 * t261 + Ifges(3,5) * t184 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t190) * t183) * t190 + (t125 * t261 + Ifges(4,5) * t121 - 0.2e1 * Ifges(3,6) * t184 - Ifges(4,6) * t120 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t187 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t190) * t183) * t187) * qJD(2)) * t183 - (t9 + t10) * t193 + 0.2e1 * m(4) * (t106 * t117 + t25 * t65 + t26 * t64) + 0.2e1 * m(3) * (t116 * t125 - t117 * t124) + (t7 + t8) * t85 + (t30 + t31 - t66) * t83 + (t34 + t35) * t44 + (t32 + t33) * t43 + (t1 * t12 + t11 * t27 + t15 * t2) * t263 + (t16 * t4 + t17 * t3 + t24 * t56) * t264 + (-t46 + 0.2e1 * t244 + t270) * t120 + 0.2e1 * t15 * t18 + 0.2e1 * t17 * t19 + 0.2e1 * t12 * t20 + 0.2e1 * t16 * t21 + 0.2e1 * t27 * t13 + 0.2e1 * t11 * t49 + 0.2e1 * t24 * t50 + 0.2e1 * t56 * t14 + 0.2e1 * t2 * t58 + 0.2e1 * t3 * t59 + 0.2e1 * t1 * t60 + 0.2e1 * t4 * t61 + 0.2e1 * t65 * t68 + 0.2e1 * t64 * t69 + t84 * t67 + 0.2e1 * t25 * t88 + 0.2e1 * t26 * t89 + 0.2e1 * t106 * t48 + (t47 + 0.2e1 * t243) * t121 + (t162 - 0.2e1 * t245 - 0.2e1 * t246) * t184; (-m(4) * t117 - t48) * pkin(2) - t215 * t193 + t162 + m(5) * (t104 * t4 + t105 * t3 + t16 * t63 + t17 * t62) - t246 - t245 + m(6) * (t1 * t78 + t100 * t27 + t11 * t147 + t12 * t37 + t15 * t45 + t2 * t87) + (t109 / 0.2e1 + t110 / 0.2e1 - t157 / 0.2e1) * t83 + (Ifges(4,5) * t199 - t26 * mrSges(4,3) + t47 / 0.2e1 + t243 + t221 * t188 + t223 * t185 + (-t185 * t217 + t188 * t218) * qJD(4) + (-t65 * mrSges(4,3) + t30 / 0.2e1 + t31 / 0.2e1 - t66 / 0.2e1 + t219 / 0.2e1) * qJD(3) + (-qJD(3) * t88 - t69 + m(4) * (-qJD(3) * t65 - t26) - t265) * pkin(8)) * t186 + (t70 / 0.2e1 + t71 / 0.2e1 - t138 / 0.2e1) * t120 + t45 * t58 + t37 * t60 + t62 * t59 + t63 * t61 + t27 * t76 + t56 * t77 + t78 * t20 + t87 * t18 + t12 * t93 + t16 * t94 + t15 * t95 + t17 * t96 + t100 * t49 + t104 * t21 + t105 * t19 + t207 * t44 + t208 * t43 + t11 * t126 + t24 * t127 + t106 * t133 + t121 * t141 / 0.2e1 + t2 * t142 + t3 * t143 + t1 * t144 + t4 * t145 + t147 * t13 + t84 * t160 / 0.2e1 + t216 * t85 + (-t190 * t182 / 0.2e1 - Ifges(3,6) * t230) * t183 + (Ifges(4,6) * t199 + t25 * mrSges(4,3) - t5 / 0.2e1 - t6 / 0.2e1 + t46 / 0.2e1 - t244 + (m(4) * t25 + t68) * pkin(8) + (-t64 * mrSges(4,3) + t67 / 0.2e1 + t217 * t188 + t218 * t185 + (-m(4) * t64 + m(5) * t56 + t50 - t89) * pkin(8)) * qJD(3)) * t189; -0.2e1 * pkin(2) * t133 + 0.2e1 * t100 * t126 + 0.2e1 * t104 * t94 + 0.2e1 * t105 * t96 + 0.2e1 * t45 * t142 + 0.2e1 * t62 * t143 + 0.2e1 * t37 * t144 + 0.2e1 * t63 * t145 + 0.2e1 * t147 * t76 + 0.2e1 * t78 * t93 + 0.2e1 * t87 * t95 + (t100 * t147 + t37 * t78 + t45 * t87) * t263 + (t104 * t63 + t105 * t62) * t264 + (t138 - t70 - t71 + (t127 * t262 + t185 * t237 + t188 * t236 + t160) * qJD(3)) * t189 + (t77 * t262 + t141 + (t74 + t75) * t188 + (-t72 - t73) * t185 + (-t185 * t236 + t188 * t237) * qJD(4) + (pkin(8) ^ 2 * t189 * t264 + t109 + t110 - t157) * qJD(3)) * t186; m(6) * (t1 * t149 + t11 * t175 + t118 * t15 + t119 * t12 - t152 * t2) + t220 - t25 * mrSges(4,2) + t26 * mrSges(4,1) + t201 * t44 - t202 * t43 + t203 * t83 - t204 * t193 + t205 * t85 + t206 * t120 + t118 * t58 + t119 * t60 + t27 * t131 + t56 * t132 + t149 * t20 + t11 * t150 + t24 * t151 - t152 * t18 + t175 * t13 + (t2 * mrSges(6,3) + t3 * mrSges(5,3) + (-t16 * mrSges(5,3) - t12 * mrSges(6,3) + t217) * qJD(4) + (-qJD(4) * t61 + t19 + m(5) * (-qJD(4) * t16 + t3)) * pkin(9) - t223) * t188 + (-t1 * mrSges(6,3) - t4 * mrSges(5,3) + (-m(5) * t4 - t21) * pkin(9) + (-t15 * mrSges(6,3) + pkin(4) * t49 - pkin(9) * t59 + t17 * t266 + t27 * t258 + t218) * qJD(4) + t221) * t185 + t265 * pkin(3); t182 + m(6) * (t100 * t175 + t118 * t87 + t119 * t78 + t149 * t37 - t152 * t45) + t186 * pkin(8) * t132 - pkin(3) * t77 + t118 * t142 + t119 * t144 + t147 * t131 + t149 * t93 + t100 * t150 - t152 * t95 + t175 * t76 - t206 * t189 + ((-Ifges(4,6) + t203) * t186 + (t186 * mrSges(4,2) + (-m(5) * pkin(3) - mrSges(4,1) + t151) * t189) * pkin(8)) * qJD(3) + (-t37 * mrSges(6,3) - t63 * mrSges(5,3) - t205 * t186 + t202 * t228 + (-m(5) * t63 - t94) * pkin(9) + (-t87 * mrSges(6,3) + pkin(4) * t126 - pkin(9) * t143 + t105 * t266 + t147 * t258 - t201 * t186 - t208) * qJD(4) + t215) * t185 + (t45 * mrSges(6,3) + t62 * mrSges(5,3) + t204 * t186 + t201 * t228 + (m(5) * t62 + t96) * pkin(9) + (-t78 * mrSges(6,3) - t104 * mrSges(5,3) + t202 * t186 + (-m(5) * t104 - t145) * pkin(9) + t207) * qJD(4) + t216) * t188; (-t118 * t152 + t119 * t149) * t263 - 0.2e1 * pkin(3) * t132 + 0.2e1 * t175 * t131 + (t119 * t260 + t139 + t140 + (-t152 * t260 - t155 - t156 + 0.2e1 * (m(6) * t175 + t150) * pkin(4)) * qJD(4)) * t185 + (0.2e1 * t118 * mrSges(6,3) + t136 + t137 + (t149 * t260 + t158 + t159) * qJD(4)) * t188; mrSges(5,1) * t4 + mrSges(6,1) * t1 - mrSges(5,2) * t3 - mrSges(6,2) * t2 + (m(6) * t1 + t20) * pkin(4) + t270; mrSges(5,1) * t63 + mrSges(6,1) * t37 - mrSges(5,2) * t62 - mrSges(6,2) * t45 - t271 * t210 + (m(6) * t37 + t93) * pkin(4) + (-t273 - t274) * t226 + t232 + t233; -mrSges(6,2) * t118 + t180 + t181 + (mrSges(6,1) + t258) * t119 + ((-mrSges(5,1) * pkin(9) - (mrSges(6,3) * pkin(4))) * t188 + (mrSges(5,2) * pkin(9) - t271) * t185) * qJD(4); 0; m(6) * t11 + t13; m(6) * t100 + t76; t227 * t258 + t131; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t22(1), t22(2), t22(4), t22(7), t22(11); t22(2), t22(3), t22(5), t22(8), t22(12); t22(4), t22(5), t22(6), t22(9), t22(13); t22(7), t22(8), t22(9), t22(10), t22(14); t22(11), t22(12), t22(13), t22(14), t22(15);];
Mq = res;
