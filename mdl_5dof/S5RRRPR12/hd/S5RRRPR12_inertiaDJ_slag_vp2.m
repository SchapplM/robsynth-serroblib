% Calculate time derivative of joint inertia matrix for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:07
% EndTime: 2019-12-31 21:37:18
% DurationCPUTime: 3.90s
% Computational Cost: add. (4832->497), mult. (13219->758), div. (0->0), fcn. (12329->10), ass. (0->216)
t196 = sin(pkin(10));
t198 = cos(pkin(10));
t264 = mrSges(5,1) * t196 + mrSges(5,2) * t198;
t263 = 2 * pkin(8);
t200 = sin(qJ(5));
t203 = cos(qJ(5));
t206 = t196 * t200 - t198 * t203;
t244 = -t206 / 0.2e1;
t164 = t196 * t203 + t198 * t200;
t243 = t164 / 0.2e1;
t241 = t198 / 0.2e1;
t197 = sin(pkin(5));
t202 = sin(qJ(2));
t220 = qJD(2) * t202;
t212 = t197 * t220;
t199 = cos(pkin(5));
t205 = cos(qJ(2));
t225 = t197 * t205;
t159 = pkin(1) * t199 * t202 + pkin(7) * t225;
t140 = pkin(8) * t199 + t159;
t141 = (-pkin(2) * t205 - pkin(8) * t202 - pkin(1)) * t197;
t221 = qJD(2) * t197;
t145 = (pkin(2) * t202 - pkin(8) * t205) * t221;
t226 = t197 * t202;
t185 = pkin(7) * t226;
t240 = pkin(1) * t205;
t158 = t199 * t240 - t185;
t146 = t158 * qJD(2);
t201 = sin(qJ(3));
t204 = cos(qJ(3));
t218 = qJD(3) * t204;
t219 = qJD(3) * t201;
t45 = -t140 * t218 - t141 * t219 + t145 * t204 - t146 * t201;
t38 = -pkin(3) * t212 - t45;
t153 = -t199 * t204 + t201 * t226;
t211 = t205 * t221;
t112 = -qJD(3) * t153 + t204 * t211;
t82 = -t112 * t196 + t198 * t212;
t83 = t112 * t198 + t196 * t212;
t43 = -mrSges(5,1) * t82 + mrSges(5,2) * t83;
t262 = -m(5) * t38 - t43;
t151 = t206 * qJD(5);
t261 = 0.2e1 * m(5);
t260 = 2 * m(6);
t259 = -2 * mrSges(3,3);
t258 = m(5) * pkin(8);
t154 = t199 * t201 + t204 * t226;
t109 = -t154 * t196 - t198 * t225;
t110 = t154 * t198 - t196 * t225;
t59 = t109 * t203 - t110 * t200;
t257 = t59 / 0.2e1;
t60 = t109 * t200 + t110 * t203;
t256 = t60 / 0.2e1;
t255 = t82 / 0.2e1;
t254 = t83 / 0.2e1;
t102 = Ifges(6,4) * t164 - Ifges(6,2) * t206;
t252 = t102 / 0.2e1;
t103 = Ifges(6,1) * t164 - Ifges(6,4) * t206;
t251 = t103 / 0.2e1;
t233 = Ifges(5,4) * t198;
t208 = -Ifges(5,2) * t196 + t233;
t124 = (Ifges(5,6) * t201 + t204 * t208) * qJD(3);
t250 = t124 / 0.2e1;
t234 = Ifges(5,4) * t196;
t209 = Ifges(5,1) * t198 - t234;
t125 = (Ifges(5,5) * t201 + t204 * t209) * qJD(3);
t249 = t125 / 0.2e1;
t142 = t164 * t201;
t248 = -t142 / 0.2e1;
t143 = t206 * t201;
t247 = -t143 / 0.2e1;
t246 = -t151 / 0.2e1;
t152 = t164 * qJD(5);
t245 = -t152 / 0.2e1;
t242 = -t196 / 0.2e1;
t239 = pkin(9) + qJ(4);
t44 = -t140 * t219 + t141 * t218 + t145 * t201 + t146 * t204;
t35 = (qJ(4) * t220 - qJD(4) * t205) * t197 + t44;
t111 = qJD(3) * t154 + t201 * t211;
t147 = t159 * qJD(2);
t39 = pkin(3) * t111 - qJ(4) * t112 - qJD(4) * t154 + t147;
t12 = t196 * t39 + t198 * t35;
t139 = t185 + (-pkin(2) - t240) * t199;
t66 = t153 * pkin(3) - t154 * qJ(4) + t139;
t74 = t140 * t204 + t141 * t201;
t67 = -qJ(4) * t225 + t74;
t33 = t196 * t66 + t198 * t67;
t236 = Ifges(4,4) * t201;
t235 = Ifges(4,4) * t204;
t232 = t146 * mrSges(3,2);
t231 = t147 * mrSges(3,1);
t230 = t147 * mrSges(4,1);
t229 = t147 * mrSges(4,2);
t228 = t196 * t201;
t227 = t196 * t204;
t224 = t198 * t201;
t223 = t198 * t204;
t96 = -Ifges(6,5) * t151 - Ifges(6,6) * t152;
t144 = t264 * t218;
t150 = -qJD(4) * t201 + (pkin(3) * t201 - qJ(4) * t204) * qJD(3);
t216 = pkin(8) * t219;
t107 = t150 * t198 + t196 * t216;
t171 = -pkin(3) * t204 - qJ(4) * t201 - pkin(2);
t130 = pkin(8) * t223 + t171 * t196;
t19 = qJD(5) * t59 + t200 * t82 + t203 * t83;
t20 = -qJD(5) * t60 - t200 * t83 + t203 * t82;
t3 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t111;
t89 = -t152 * t201 - t206 * t218;
t90 = t151 * t201 - t164 * t218;
t40 = Ifges(6,5) * t89 + Ifges(6,6) * t90 + Ifges(6,3) * t219;
t215 = Ifges(4,6) * t225;
t214 = Ifges(4,5) * t112 - Ifges(4,6) * t111 + Ifges(4,3) * t212;
t213 = pkin(4) * t196 + pkin(8);
t210 = Ifges(5,5) * t196 / 0.2e1 + Ifges(5,6) * t241 + Ifges(6,5) * t243 + Ifges(6,6) * t244;
t6 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t46 = -mrSges(6,1) * t90 + mrSges(6,2) * t89;
t11 = -t196 * t35 + t198 * t39;
t32 = -t196 * t67 + t198 * t66;
t73 = -t140 * t201 + t141 * t204;
t68 = pkin(3) * t225 - t73;
t207 = Ifges(5,5) * t198 - Ifges(5,6) * t196;
t18 = pkin(4) * t153 - pkin(9) * t110 + t32;
t26 = pkin(9) * t109 + t33;
t7 = t18 * t203 - t200 * t26;
t8 = t18 * t200 + t203 * t26;
t113 = -pkin(9) * t228 + t130;
t162 = t198 * t171;
t99 = -pkin(9) * t224 + t162 + (-pkin(8) * t196 - pkin(4)) * t204;
t58 = t113 * t203 + t200 * t99;
t57 = -t113 * t200 + t203 * t99;
t172 = t239 * t196;
t174 = t239 * t198;
t114 = -t172 * t203 - t174 * t200;
t115 = -t172 * t200 + t174 * t203;
t193 = Ifges(4,5) * t218;
t191 = -pkin(4) * t198 - pkin(3);
t183 = Ifges(3,5) * t211;
t179 = Ifges(4,1) * t201 + t235;
t178 = Ifges(4,2) * t204 + t236;
t177 = Ifges(5,1) * t196 + t233;
t176 = Ifges(5,2) * t198 + t234;
t173 = -mrSges(5,1) * t198 + mrSges(5,2) * t196;
t170 = t213 * t201;
t169 = (Ifges(4,1) * t204 - t236) * qJD(3);
t168 = (-Ifges(4,2) * t201 + t235) * qJD(3);
t167 = (mrSges(4,1) * t201 + mrSges(4,2) * t204) * qJD(3);
t166 = -mrSges(5,1) * t204 - mrSges(5,3) * t224;
t165 = mrSges(5,2) * t204 - mrSges(5,3) * t228;
t160 = t213 * t218;
t157 = (mrSges(5,1) * t201 - mrSges(5,3) * t223) * qJD(3);
t156 = (-mrSges(5,2) * t201 - mrSges(5,3) * t227) * qJD(3);
t155 = t264 * t201;
t138 = -Ifges(5,5) * t204 + t201 * t209;
t137 = -Ifges(5,6) * t204 + t201 * t208;
t136 = -Ifges(5,3) * t204 + t201 * t207;
t134 = t196 * t150;
t129 = -pkin(8) * t227 + t162;
t123 = (Ifges(5,3) * t201 + t204 * t207) * qJD(3);
t119 = -mrSges(6,1) * t204 + mrSges(6,3) * t143;
t118 = mrSges(6,2) * t204 - mrSges(6,3) * t142;
t117 = -mrSges(4,1) * t225 - mrSges(4,3) * t154;
t116 = mrSges(4,2) * t225 - mrSges(4,3) * t153;
t108 = -t198 * t216 + t134;
t100 = mrSges(6,1) * t206 + mrSges(6,2) * t164;
t98 = -Ifges(6,1) * t151 - Ifges(6,4) * t152;
t97 = -Ifges(6,4) * t151 - Ifges(6,2) * t152;
t95 = mrSges(6,1) * t152 - mrSges(6,2) * t151;
t94 = t134 + (-pkin(8) * t224 - pkin(9) * t227) * qJD(3);
t93 = mrSges(4,1) * t212 - mrSges(4,3) * t112;
t92 = -mrSges(4,2) * t212 - mrSges(4,3) * t111;
t91 = mrSges(6,1) * t142 - mrSges(6,2) * t143;
t88 = Ifges(4,1) * t154 - Ifges(4,4) * t153 - Ifges(4,5) * t225;
t87 = Ifges(4,4) * t154 - Ifges(4,2) * t153 - t215;
t81 = (pkin(4) * t201 - pkin(9) * t223) * qJD(3) + t107;
t80 = -Ifges(6,1) * t143 - Ifges(6,4) * t142 - Ifges(6,5) * t204;
t79 = -Ifges(6,4) * t143 - Ifges(6,2) * t142 - Ifges(6,6) * t204;
t78 = -Ifges(6,5) * t143 - Ifges(6,6) * t142 - Ifges(6,3) * t204;
t76 = -qJD(4) * t164 - qJD(5) * t115;
t75 = -qJD(4) * t206 + qJD(5) * t114;
t72 = mrSges(5,1) * t153 - mrSges(5,3) * t110;
t71 = -mrSges(5,2) * t153 + mrSges(5,3) * t109;
t70 = -mrSges(6,2) * t219 + mrSges(6,3) * t90;
t69 = mrSges(6,1) * t219 - mrSges(6,3) * t89;
t62 = mrSges(4,1) * t111 + mrSges(4,2) * t112;
t61 = -mrSges(5,1) * t109 + mrSges(5,2) * t110;
t56 = Ifges(4,1) * t112 - Ifges(4,4) * t111 + Ifges(4,5) * t212;
t55 = Ifges(4,4) * t112 - Ifges(4,2) * t111 + Ifges(4,6) * t212;
t54 = mrSges(5,1) * t111 - mrSges(5,3) * t83;
t53 = -mrSges(5,2) * t111 + mrSges(5,3) * t82;
t52 = Ifges(5,1) * t110 + Ifges(5,4) * t109 + Ifges(5,5) * t153;
t51 = Ifges(5,4) * t110 + Ifges(5,2) * t109 + Ifges(5,6) * t153;
t50 = Ifges(5,5) * t110 + Ifges(5,6) * t109 + Ifges(5,3) * t153;
t49 = -pkin(4) * t109 + t68;
t48 = mrSges(6,1) * t153 - mrSges(6,3) * t60;
t47 = -mrSges(6,2) * t153 + mrSges(6,3) * t59;
t42 = Ifges(6,1) * t89 + Ifges(6,4) * t90 + Ifges(6,5) * t219;
t41 = Ifges(6,4) * t89 + Ifges(6,2) * t90 + Ifges(6,6) * t219;
t31 = Ifges(5,1) * t83 + Ifges(5,4) * t82 + Ifges(5,5) * t111;
t30 = Ifges(5,4) * t83 + Ifges(5,2) * t82 + Ifges(5,6) * t111;
t29 = Ifges(5,5) * t83 + Ifges(5,6) * t82 + Ifges(5,3) * t111;
t28 = -mrSges(6,1) * t59 + mrSges(6,2) * t60;
t27 = -pkin(4) * t82 + t38;
t25 = Ifges(6,1) * t60 + Ifges(6,4) * t59 + Ifges(6,5) * t153;
t24 = Ifges(6,4) * t60 + Ifges(6,2) * t59 + Ifges(6,6) * t153;
t23 = Ifges(6,5) * t60 + Ifges(6,6) * t59 + Ifges(6,3) * t153;
t22 = -qJD(5) * t58 - t200 * t94 + t203 * t81;
t21 = qJD(5) * t57 + t200 * t81 + t203 * t94;
t14 = -mrSges(6,2) * t111 + mrSges(6,3) * t20;
t13 = mrSges(6,1) * t111 - mrSges(6,3) * t19;
t10 = pkin(9) * t82 + t12;
t9 = pkin(4) * t111 - pkin(9) * t83 + t11;
t5 = Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t111;
t4 = Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t111;
t2 = -qJD(5) * t8 - t10 * t200 + t203 * t9;
t1 = qJD(5) * t7 + t10 * t203 + t200 * t9;
t15 = [0.2e1 * m(4) * (t139 * t147 + t44 * t74 + t45 * t73) + 0.2e1 * m(3) * (t146 * t159 - t147 * t158) + (t50 + t23 - t87) * t111 + 0.2e1 * t139 * t62 + 0.2e1 * t44 * t116 + 0.2e1 * t45 * t117 + t109 * t30 + t110 * t31 + t112 * t88 + 0.2e1 * t74 * t92 + 0.2e1 * t73 * t93 + t82 * t51 + t83 * t52 + 0.2e1 * t12 * t71 + 0.2e1 * t11 * t72 + 0.2e1 * t38 * t61 + 0.2e1 * t68 * t43 + t59 * t4 + t60 * t5 + 0.2e1 * t1 * t47 + 0.2e1 * t2 * t48 + 0.2e1 * t49 * t6 + 0.2e1 * t33 * t53 + 0.2e1 * t32 * t54 + t20 * t24 + t19 * t25 + 0.2e1 * t27 * t28 + 0.2e1 * t7 * t13 + 0.2e1 * t8 * t14 + (-t205 * t214 + 0.2e1 * (t146 * t205 + t147 * t202) * mrSges(3,3) + ((t158 * t259 + Ifges(3,5) * t199 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t205) * t197) * t205 + (t159 * t259 + Ifges(4,5) * t154 - 0.2e1 * Ifges(3,6) * t199 - Ifges(4,6) * t153 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t202 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t205) * t197) * t202) * qJD(2)) * t197 + (t56 + 0.2e1 * t229) * t154 + (t29 + t3 - t55 + 0.2e1 * t230) * t153 + (t183 - 0.2e1 * t231 - 0.2e1 * t232) * t199 + (t1 * t8 + t2 * t7 + t27 * t49) * t260 + (t11 * t32 + t12 * t33 + t38 * t68) * t261; (-m(4) * t147 - t62) * pkin(2) + t183 + m(5) * (t107 * t32 + t108 * t33 + t11 * t129 + t12 * t130) - t231 - t232 + m(6) * (t1 * t58 + t160 * t49 + t170 * t27 + t2 * t57 + t21 * t8 + t22 * t7) + t139 * t167 + t154 * t169 / 0.2e1 + t170 * t6 + t112 * t179 / 0.2e1 + t33 * t156 + t32 * t157 + t160 * t28 + t12 * t165 + t11 * t166 + t38 * t155 + t68 * t144 + t129 * t54 + t130 * t53 + t1 * t118 + t2 * t119 + t107 * t72 + t108 * t71 + t90 * t24 / 0.2e1 + t27 * t91 + t19 * t80 / 0.2e1 + t89 * t25 / 0.2e1 + t20 * t79 / 0.2e1 + t7 * t69 + t8 * t70 + t57 * t13 + t58 * t14 + t21 * t47 + t22 * t48 + t49 * t46 + (t204 * t92 + (t43 - t93) * t201 + (-t201 * t116 + (-t117 + t61) * t204) * qJD(3) + m(5) * (t201 * t38 + t218 * t68) + m(4) * (-t201 * t45 + t204 * t44 - t218 * t73 - t219 * t74)) * pkin(8) + (-t205 * t193 / 0.2e1 + (Ifges(4,5) * t201 / 0.2e1 + Ifges(4,6) * t204 / 0.2e1 - Ifges(3,6)) * t220) * t197 + (-t230 + t55 / 0.2e1 - t29 / 0.2e1 - t3 / 0.2e1 + t44 * mrSges(4,3)) * t204 + ((t52 * t241 + t51 * t242 - t73 * mrSges(4,3) + t88 / 0.2e1) * t204 + (-t74 * mrSges(4,3) + t23 / 0.2e1 + t50 / 0.2e1 - t87 / 0.2e1 + t215 / 0.2e1) * t201) * qJD(3) + (t229 + t56 / 0.2e1 - t45 * mrSges(4,3) + t30 * t242 + t31 * t241) * t201 + t5 * t247 + t4 * t248 + t110 * t249 + t109 * t250 + t138 * t254 + t137 * t255 + t42 * t256 + t41 * t257 + (-t178 / 0.2e1 + t136 / 0.2e1 + t78 / 0.2e1) * t111 + (-t168 / 0.2e1 + t123 / 0.2e1 + t40 / 0.2e1) * t153; -0.2e1 * pkin(2) * t167 + 0.2e1 * t107 * t166 + 0.2e1 * t108 * t165 + 0.2e1 * t21 * t118 + 0.2e1 * t22 * t119 + 0.2e1 * t129 * t157 + 0.2e1 * t130 * t156 - t142 * t41 - t143 * t42 + 0.2e1 * t160 * t91 + 0.2e1 * t170 * t46 + 0.2e1 * t57 * t69 + 0.2e1 * t58 * t70 + t90 * t79 + t89 * t80 + (t160 * t170 + t21 * t58 + t22 * t57) * t260 + (t107 * t129 + t108 * t130) * t261 + (t168 - t123 - t40) * t204 + (-t196 * t124 + t198 * t125 + t144 * t263 + t169) * t201 + ((t136 - t178 + t78) * t201 + (-t196 * t137 + t198 * t138 + t179 + (t201 * t258 + t155) * t263) * t204) * qJD(3); (t30 / 0.2e1 + t12 * mrSges(5,3) + qJ(4) * t53 + qJD(4) * t71 + m(5) * (qJ(4) * t12 + qJD(4) * t33)) * t198 + (t31 / 0.2e1 - t11 * mrSges(5,3) - qJ(4) * t54 - qJD(4) * t72 + m(5) * (-qJ(4) * t11 - qJD(4) * t32)) * t196 + t262 * pkin(3) + t210 * t111 + (-t1 * t206 + t151 * t7 - t152 * t8 - t164 * t2) * mrSges(6,3) + t214 + t191 * t6 + t38 * t173 + t153 * t96 / 0.2e1 + t114 * t13 + t115 * t14 + t27 * t100 + t49 * t95 + m(6) * (t1 * t115 + t114 * t2 + t191 * t27 + t7 * t76 + t75 * t8) + t75 * t47 + t76 * t48 + t45 * mrSges(4,1) - t44 * mrSges(4,2) + t5 * t243 + t4 * t244 + t24 * t245 + t25 * t246 + t19 * t251 + t20 * t252 + t177 * t254 + t176 * t255 + t98 * t256 + t97 * t257; t193 - t204 * t96 / 0.2e1 + t191 * t46 + t170 * t95 + t160 * t100 + t41 * t244 + t42 * t243 + t80 * t246 + t79 * t245 + t98 * t247 - pkin(3) * t144 + t97 * t248 + t114 * t69 + t115 * t70 + t75 * t118 + t76 * t119 + t90 * t252 + t89 * t251 + m(6) * (t114 * t22 + t115 * t21 + t160 * t191 + t57 * t76 + t58 * t75) + (t250 + t108 * mrSges(5,3) + qJ(4) * t156 + qJD(4) * t165 + m(5) * (qJ(4) * t108 + qJD(4) * t130)) * t198 + (t249 - t107 * mrSges(5,3) - qJ(4) * t157 - qJD(4) * t166 + m(5) * (-qJ(4) * t107 - qJD(4) * t129)) * t196 + (t151 * t57 - t152 * t58 - t164 * t22 - t206 * t21) * mrSges(6,3) + ((mrSges(4,2) * pkin(8) - Ifges(4,6) + t210) * t201 + (t177 * t241 + t176 * t242 + (-m(5) * pkin(3) - mrSges(4,1) + t173) * pkin(8)) * t204) * qJD(3); (t114 * t76 + t115 * t75) * t260 - t151 * t103 + t164 * t98 - t152 * t102 - t206 * t97 + 0.2e1 * t191 * t95 + 0.2e1 * (t114 * t151 - t115 * t152 - t164 * t76 - t206 * t75) * mrSges(6,3) + (qJ(4) * t261 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t196 ^ 2 + t198 ^ 2); m(6) * t27 - t262 + t6; m(6) * t160 + t218 * t258 + t144 + t46; t95; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t22 - mrSges(6,2) * t21 + t40; mrSges(6,1) * t76 - mrSges(6,2) * t75 + t96; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
