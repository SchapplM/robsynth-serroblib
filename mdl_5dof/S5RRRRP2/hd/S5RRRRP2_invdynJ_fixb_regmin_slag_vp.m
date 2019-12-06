% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:15
% EndTime: 2019-12-05 18:48:21
% DurationCPUTime: 1.68s
% Computational Cost: add. (2495->256), mult. (3783->314), div. (0->0), fcn. (2445->12), ass. (0->176)
t242 = cos(qJ(4));
t147 = sin(qJ(4));
t151 = cos(qJ(3));
t142 = qJD(1) + qJD(2);
t149 = sin(qJ(2));
t241 = pkin(1) * t149;
t204 = qJD(1) * t241;
t92 = pkin(7) * t142 + t204;
t191 = pkin(8) * t142 + t92;
t55 = t191 * t151;
t48 = t147 * t55;
t148 = sin(qJ(3));
t54 = t191 * t148;
t51 = qJD(3) * pkin(3) - t54;
t189 = t242 * t51 - t48;
t197 = t242 * t148;
t84 = t147 * t151 + t197;
t63 = t84 * t142;
t222 = t63 * qJ(5);
t252 = t222 - t189;
t145 = qJ(3) + qJ(4);
t131 = sin(t145);
t133 = cos(t145);
t146 = qJ(1) + qJ(2);
t132 = sin(t146);
t134 = cos(t146);
t215 = -g(2) * t132 + g(3) * t134;
t251 = -g(1) * t133 + t215 * t131;
t250 = g(2) * t134 + g(3) * t132;
t190 = qJD(4) * t242;
t196 = t242 * t151;
t249 = -qJD(3) * t196 - t151 * t190;
t136 = t151 * pkin(3);
t234 = pkin(2) + t136;
t152 = cos(qJ(2));
t211 = qJD(1) * t152;
t203 = pkin(1) * t211;
t244 = -pkin(7) - pkin(8);
t108 = t244 * t148;
t135 = t151 * pkin(8);
t109 = pkin(7) * t151 + t135;
t225 = t147 * t108 + t242 * t109;
t198 = qJD(3) * t244;
t89 = t148 * t198;
t90 = t151 * t198;
t248 = -t225 * qJD(4) - t147 * t89 + t84 * t203 + t242 * t90;
t218 = t147 * t148;
t168 = t196 - t218;
t206 = qJD(4) * t147;
t247 = -t108 * t190 + t109 * t206 - t147 * t90 + t168 * t203 - t242 * t89;
t208 = qJD(3) * t148;
t129 = pkin(3) * t208;
t246 = t129 - t204;
t141 = qJD(3) + qJD(4);
t245 = t63 ^ 2;
t243 = pkin(4) * t168;
t240 = pkin(1) * t152;
t139 = qJDD(1) + qJDD(2);
t239 = pkin(2) * t139;
t238 = pkin(2) * t142;
t200 = t142 * t218;
t61 = -t142 * t196 + t200;
t65 = -t142 * t234 - t203;
t35 = pkin(4) * t61 + qJD(5) + t65;
t236 = t35 * t63;
t235 = t63 * t61;
t124 = pkin(7) + t241;
t233 = -pkin(8) - t124;
t41 = t141 * t84;
t229 = -t41 * qJ(5) + qJD(5) * t168;
t232 = t229 - t247;
t174 = t141 * t218;
t40 = t174 + t249;
t172 = t40 * qJ(5) - t84 * qJD(5);
t231 = t172 + t248;
t16 = pkin(4) * t141 - t252;
t230 = t16 + t252;
t228 = -t242 * t54 - t48;
t207 = qJD(3) * t151;
t210 = qJD(2) * t149;
t130 = pkin(1) * t210;
t213 = -qJD(1) * t130 + qJDD(1) * t240;
t66 = -t213 - t239;
t93 = -t203 - t238;
t227 = t66 * t148 + t93 * t207;
t80 = t233 * t148;
t81 = t124 * t151 + t135;
t226 = t147 * t80 + t242 * t81;
t224 = qJ(5) * t84;
t223 = t61 * qJ(5);
t221 = t132 * t133;
t220 = t133 * t134;
t219 = t142 * t148;
t217 = t148 * t139;
t216 = t151 * t139;
t214 = pkin(4) * t133 + t136;
t143 = t148 ^ 2;
t212 = -t151 ^ 2 + t143;
t209 = qJD(2) * t152;
t205 = qJDD(1) * t149;
t202 = pkin(1) * t209;
t201 = t250 * t151 + t93 * t208;
t50 = t242 * t55;
t126 = -pkin(2) - t240;
t195 = t142 * t210;
t194 = t142 * t208;
t193 = t142 * t207;
t192 = pkin(4) * t41 + t129;
t188 = t147 * t54 - t50;
t187 = -t147 * t81 + t242 * t80;
t186 = qJD(3) * t233;
t185 = t242 * t108 - t109 * t147;
t140 = -qJ(5) + t244;
t88 = pkin(2) + t214;
t184 = t132 * t140 - t134 * t88;
t38 = pkin(3) * t194 - t139 * t234 - t213;
t183 = g(2) * t220 + g(3) * t221 - t168 * t38 + t65 * t41;
t182 = -t139 * t197 + t249 * t142 - t147 * t216;
t181 = t142 * t204;
t180 = -t213 - t250;
t103 = t126 - t136;
t175 = -t139 * t196 + t147 * t217;
t170 = -t147 * t51 - t50;
t18 = -t170 - t223;
t138 = qJDD(3) + qJDD(4);
t67 = pkin(7) * t139 + (qJD(1) * t209 + t205) * pkin(1);
t25 = -t92 * t207 + qJDD(3) * pkin(3) - t148 * t67 + (-t193 - t217) * pkin(8);
t26 = -t92 * t208 + t151 * t67 + (-t194 + t216) * pkin(8);
t160 = t170 * qJD(4) - t147 * t26 + t242 * t25;
t21 = t142 * t174 + t182;
t3 = t138 * pkin(4) + t21 * qJ(5) - t63 * qJD(5) + t160;
t158 = t147 * t25 + t51 * t190 - t55 * t206 + t242 * t26;
t22 = t41 * t142 + t175;
t4 = -t22 * qJ(5) - t61 * qJD(5) + t158;
t173 = t16 * t40 + t168 * t4 - t18 * t41 - t3 * t84 - t215;
t171 = -t132 * t88 - t134 * t140;
t52 = t148 * t186 + t151 * t202;
t53 = -t148 * t202 + t151 * t186;
t169 = t147 * t53 + t80 * t190 - t81 * t206 + t242 * t52;
t166 = -t142 * t93 + t215 - t67;
t165 = -t131 * t250 + t38 * t84 - t65 * t40;
t154 = qJD(3) ^ 2;
t164 = -pkin(7) * t154 + t181 + t239;
t163 = -pkin(1) * t195 - t124 * t154 - t126 * t139;
t162 = -pkin(7) * qJDD(3) + (t203 - t238) * qJD(3);
t161 = -qJDD(3) * t124 + (t126 * t142 - t202) * qJD(3);
t159 = -t226 * qJD(4) - t147 * t52 + t242 * t53;
t13 = pkin(4) * t22 + qJDD(5) + t38;
t156 = g(1) * t131 - g(2) * t221 + g(3) * t220 + t65 * t61 - t158;
t155 = -t65 * t63 + t160 + t251;
t153 = cos(qJ(1));
t150 = sin(qJ(1));
t137 = t142 ^ 2;
t125 = t242 * pkin(3) + pkin(4);
t102 = qJDD(3) * t151 - t148 * t154;
t101 = qJDD(3) * t148 + t151 * t154;
t91 = t130 + t129;
t79 = t168 * qJ(5);
t68 = t139 * t143 + 0.2e1 * t148 * t193;
t60 = t61 ^ 2;
t42 = -0.2e1 * t212 * t142 * qJD(3) + 0.2e1 * t148 * t216;
t37 = t79 + t225;
t36 = t185 - t224;
t30 = t79 + t226;
t29 = t187 - t224;
t28 = t138 * t168 - t141 * t41;
t27 = t138 * t84 - t141 * t40;
t24 = -t60 + t245;
t20 = -t222 + t228;
t19 = t188 + t223;
t14 = -t182 + (-t200 + t61) * t141;
t8 = -t21 * t84 - t40 * t63;
t7 = t159 + t172;
t6 = t169 + t229;
t5 = -t168 * t21 - t22 * t84 + t40 * t61 - t41 * t63;
t1 = [qJDD(1), g(2) * t153 + g(3) * t150, -g(2) * t150 + g(3) * t153, t139, (t139 * t152 - t195) * pkin(1) - t180, ((-qJDD(1) - t139) * t149 + (-qJD(1) - t142) * t209) * pkin(1) + t215, t68, t42, t101, t102, 0, t161 * t148 + (t163 - t66) * t151 + t201, t161 * t151 + (-t163 - t250) * t148 + t227, t8, t5, t27, t28, 0, t103 * t22 + t187 * t138 + t159 * t141 + t91 * t61 + t183, -t103 * t21 - t226 * t138 - t169 * t141 + t91 * t63 + t165, t21 * t29 - t22 * t30 - t6 * t61 - t63 * t7 + t173, t4 * t30 + t18 * t6 + t3 * t29 + t16 * t7 + t13 * (t103 - t243) + t35 * (t130 + t192) - g(2) * (-pkin(1) * t153 + t184) - g(3) * (-pkin(1) * t150 + t171); 0, 0, 0, t139, -t180 + t181, (-t205 + (-qJD(2) + t142) * t211) * pkin(1) + t215, t68, t42, t101, t102, 0, t162 * t148 + (t164 - t66) * t151 + t201, t162 * t151 + (-t164 - t250) * t148 + t227, t8, t5, t27, t28, 0, t185 * t138 + t248 * t141 - t22 * t234 + t246 * t61 + t183, -t225 * t138 + t247 * t141 + t21 * t234 + t246 * t63 + t165, t21 * t36 - t22 * t37 - t231 * t63 - t232 * t61 + t173, t4 * t37 + t3 * t36 + t13 * (-t234 - t243) - g(2) * t184 - g(3) * t171 + (t192 - t204) * t35 + t232 * t18 + t231 * t16; 0, 0, 0, 0, 0, 0, -t148 * t137 * t151, t212 * t137, t217, t216, qJDD(3), -g(1) * t151 + t166 * t148, g(1) * t148 + t166 * t151, t235, t24, t14, -t175, t138, -t188 * t141 + (t242 * t138 - t141 * t206 - t61 * t219) * pkin(3) + t155, t228 * t141 + (-t147 * t138 - t141 * t190 - t63 * t219) * pkin(3) + t156, t125 * t21 + (t18 + t19) * t63 + (-t16 + t20) * t61 + (-t147 * t22 + (t147 * t63 - t242 * t61) * qJD(4)) * pkin(3), t3 * t125 - t18 * t20 - t16 * t19 - pkin(4) * t236 - g(1) * t214 - t215 * (-pkin(3) * t148 - pkin(4) * t131) + (-t35 * t219 + t4 * t147 + (-t147 * t16 + t242 * t18) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t24, t14, -t175, t138, -t170 * t141 + t155, t189 * t141 + t156, pkin(4) * t21 - t230 * t61, t230 * t18 + (-t236 + t251 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 - t245, t16 * t63 + t18 * t61 + t13 - t250;];
tau_reg = t1;
