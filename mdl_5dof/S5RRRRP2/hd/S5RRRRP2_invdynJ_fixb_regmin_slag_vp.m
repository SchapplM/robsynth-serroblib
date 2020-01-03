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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:48
% EndTime: 2020-01-03 12:11:52
% DurationCPUTime: 1.54s
% Computational Cost: add. (2495->259), mult. (3783->318), div. (0->0), fcn. (2445->12), ass. (0->179)
t146 = qJ(1) + qJ(2);
t132 = sin(t146);
t134 = cos(t146);
t176 = g(2) * t134 + g(3) * t132;
t149 = sin(qJ(2));
t208 = qJD(2) * t149;
t130 = pkin(1) * t208;
t152 = cos(qJ(2));
t243 = pkin(1) * t152;
t211 = -qJD(1) * t130 + qJDD(1) * t243;
t139 = qJDD(1) + qJDD(2);
t242 = pkin(2) * t139;
t66 = -t211 - t242;
t256 = t66 + t176;
t245 = cos(qJ(4));
t147 = sin(qJ(4));
t151 = cos(qJ(3));
t142 = qJD(1) + qJD(2);
t244 = pkin(1) * t149;
t202 = qJD(1) * t244;
t94 = pkin(7) * t142 + t202;
t190 = pkin(8) * t142 + t94;
t55 = t190 * t151;
t48 = t147 * t55;
t148 = sin(qJ(3));
t54 = t190 * t148;
t51 = qJD(3) * pkin(3) - t54;
t188 = t245 * t51 - t48;
t196 = t245 * t148;
t86 = t147 * t151 + t196;
t63 = t86 * t142;
t223 = t63 * qJ(5);
t255 = t223 - t188;
t189 = qJD(4) * t245;
t195 = t245 * t151;
t254 = -qJD(3) * t195 - t151 * t189;
t136 = t151 * pkin(3);
t234 = pkin(2) + t136;
t209 = qJD(1) * t152;
t201 = pkin(1) * t209;
t247 = -pkin(7) - pkin(8);
t110 = t247 * t148;
t135 = t151 * pkin(8);
t111 = pkin(7) * t151 + t135;
t226 = t147 * t110 + t245 * t111;
t197 = qJD(3) * t247;
t91 = t148 * t197;
t92 = t151 * t197;
t253 = -t226 * qJD(4) - t147 * t91 + t86 * t201 + t245 * t92;
t216 = t147 * t148;
t168 = t195 - t216;
t204 = qJD(4) * t147;
t252 = -t110 * t189 + t111 * t204 - t147 * t92 + t168 * t201 - t245 * t91;
t251 = g(2) * t132 - g(3) * t134;
t206 = qJD(3) * t148;
t129 = pkin(3) * t206;
t250 = t129 - t202;
t154 = qJD(3) ^ 2;
t249 = pkin(7) * t154 - t242;
t141 = qJD(3) + qJD(4);
t248 = t63 ^ 2;
t246 = pkin(4) * t168;
t241 = pkin(2) * t142;
t145 = qJ(3) + qJ(4);
t133 = cos(t145);
t239 = g(1) * t133;
t199 = t142 * t216;
t61 = -t142 * t195 + t199;
t65 = -t142 * t234 - t201;
t35 = pkin(4) * t61 + qJD(5) + t65;
t236 = t35 * t63;
t235 = t63 * t61;
t124 = pkin(7) + t244;
t233 = -pkin(8) - t124;
t41 = t141 * t86;
t229 = -t41 * qJ(5) + qJD(5) * t168;
t232 = t229 - t252;
t173 = t141 * t216;
t40 = t173 + t254;
t171 = t40 * qJ(5) - t86 * qJD(5);
t231 = t171 + t253;
t16 = pkin(4) * t141 - t255;
t230 = t16 + t255;
t228 = -t245 * t54 - t48;
t82 = t233 * t148;
t83 = t124 * t151 + t135;
t227 = t147 * t82 + t245 * t83;
t225 = qJ(5) * t86;
t224 = t61 * qJ(5);
t140 = -qJ(5) + t247;
t212 = pkin(4) * t133 + t136;
t90 = pkin(2) + t212;
t222 = t132 * t90 + t134 * t140;
t131 = sin(t145);
t221 = t131 * t132;
t220 = t131 * t134;
t219 = t132 * t133;
t218 = t133 * t134;
t217 = t142 * t148;
t215 = t148 * t139;
t214 = t151 * t139;
t143 = t148 ^ 2;
t210 = -t151 ^ 2 + t143;
t207 = qJD(2) * t152;
t205 = qJD(3) * t151;
t203 = qJDD(1) * t149;
t200 = pkin(1) * t207;
t50 = t245 * t55;
t126 = -pkin(2) - t243;
t194 = t142 * t208;
t193 = t142 * t206;
t192 = t142 * t205;
t191 = pkin(4) * t41 + t129;
t187 = t147 * t54 - t50;
t186 = -t147 * t83 + t245 * t82;
t185 = qJD(3) * t233;
t184 = t245 * t110 - t111 * t147;
t183 = -t132 * t140 + t134 * t90;
t38 = pkin(3) * t193 - t139 * t234 - t211;
t182 = g(2) * t220 + g(3) * t221 + t38 * t86 - t65 * t40;
t181 = -t139 * t196 + t254 * t142 - t147 * t214;
t95 = -t201 - t241;
t180 = t256 * t148 + t95 * t205;
t179 = t142 * t202;
t105 = t126 - t136;
t174 = -t139 * t195 + t147 * t215;
t170 = -t147 * t51 - t50;
t18 = -t170 - t224;
t138 = qJDD(3) + qJDD(4);
t67 = pkin(7) * t139 + (qJD(1) * t207 + t203) * pkin(1);
t25 = -t94 * t205 + qJDD(3) * pkin(3) - t148 * t67 + (-t192 - t215) * pkin(8);
t26 = -t94 * t206 + t151 * t67 + (-t193 + t214) * pkin(8);
t160 = t170 * qJD(4) - t147 * t26 + t245 * t25;
t21 = t142 * t173 + t181;
t3 = t138 * pkin(4) + t21 * qJ(5) - t63 * qJD(5) + t160;
t158 = t147 * t25 + t51 * t189 - t55 * t204 + t245 * t26;
t22 = t41 * t142 + t174;
t4 = -t22 * qJ(5) - t61 * qJD(5) + t158;
t172 = t16 * t40 + t168 * t4 - t18 * t41 - t3 * t86 - t251;
t52 = t148 * t185 + t151 * t200;
t53 = -t148 * t200 + t151 * t185;
t169 = t147 * t53 + t82 * t189 - t83 * t204 + t245 * t52;
t166 = -t142 * t95 + t251 - t67;
t165 = -g(2) * t218 - g(3) * t219 - t168 * t38 + t65 * t41;
t164 = -t176 + t179;
t163 = pkin(1) * t194 + t124 * t154 + t126 * t139;
t162 = -pkin(7) * qJDD(3) + (t201 - t241) * qJD(3);
t161 = -qJDD(3) * t124 + (t126 * t142 - t200) * qJD(3);
t159 = -t227 * qJD(4) - t147 * t52 + t245 * t53;
t13 = pkin(4) * t22 + qJDD(5) + t38;
t156 = g(1) * t131 + g(2) * t219 - g(3) * t218 + t65 * t61 - t158;
t155 = g(2) * t221 - g(3) * t220 - t65 * t63 + t160 - t239;
t153 = cos(qJ(1));
t150 = sin(qJ(1));
t137 = t142 ^ 2;
t125 = t245 * pkin(3) + pkin(4);
t104 = qJDD(3) * t151 - t148 * t154;
t103 = qJDD(3) * t148 + t151 * t154;
t93 = t130 + t129;
t81 = t168 * qJ(5);
t75 = t95 * t206;
t68 = t139 * t143 + 0.2e1 * t148 * t192;
t60 = t61 ^ 2;
t42 = -0.2e1 * t210 * t142 * qJD(3) + 0.2e1 * t148 * t214;
t37 = t81 + t226;
t36 = t184 - t225;
t30 = t81 + t227;
t29 = t186 - t225;
t28 = t138 * t168 - t141 * t41;
t27 = t138 * t86 - t141 * t40;
t24 = -t60 + t248;
t20 = -t223 + t228;
t19 = t187 + t224;
t14 = -t181 + (-t199 + t61) * t141;
t8 = -t21 * t86 - t40 * t63;
t7 = t159 + t171;
t6 = t169 + t229;
t5 = -t168 * t21 - t22 * t86 + t40 * t61 - t41 * t63;
t1 = [qJDD(1), -g(2) * t153 - g(3) * t150, g(2) * t150 - g(3) * t153, t139, (t139 * t152 - t194) * pkin(1) - t176 + t211, ((-qJDD(1) - t139) * t149 + (-qJD(1) - t142) * t207) * pkin(1) + t251, t68, t42, t103, t104, 0, t75 + t161 * t148 + (-t163 - t256) * t151, t163 * t148 + t161 * t151 + t180, t8, t5, t27, t28, 0, t105 * t22 + t186 * t138 + t159 * t141 + t93 * t61 + t165, -t105 * t21 - t227 * t138 - t169 * t141 + t93 * t63 + t182, t21 * t29 - t22 * t30 - t6 * t61 - t63 * t7 + t172, t4 * t30 + t18 * t6 + t3 * t29 + t16 * t7 + t13 * (t105 - t246) + t35 * (t130 + t191) - g(2) * (pkin(1) * t153 + t183) - g(3) * (pkin(1) * t150 + t222); 0, 0, 0, t139, t164 + t211, (-t203 + (-qJD(2) + t142) * t209) * pkin(1) + t251, t68, t42, t103, t104, 0, t75 + t162 * t148 + (t164 - t66 - t249) * t151, t162 * t151 + (-t179 + t249) * t148 + t180, t8, t5, t27, t28, 0, t184 * t138 + t253 * t141 - t22 * t234 + t250 * t61 + t165, -t226 * t138 + t252 * t141 + t21 * t234 + t250 * t63 + t182, t21 * t36 - t22 * t37 - t231 * t63 - t232 * t61 + t172, t4 * t37 + t3 * t36 + t13 * (-t234 - t246) - g(2) * t183 - g(3) * t222 + (t191 - t202) * t35 + t232 * t18 + t231 * t16; 0, 0, 0, 0, 0, 0, -t148 * t137 * t151, t210 * t137, t215, t214, qJDD(3), -g(1) * t151 + t166 * t148, g(1) * t148 + t166 * t151, t235, t24, t14, -t174, t138, -t187 * t141 + (t245 * t138 - t141 * t204 - t61 * t217) * pkin(3) + t155, t228 * t141 + (-t147 * t138 - t141 * t189 - t63 * t217) * pkin(3) + t156, t125 * t21 + (t18 + t19) * t63 + (-t16 + t20) * t61 + (-t147 * t22 + (t147 * t63 - t245 * t61) * qJD(4)) * pkin(3), t3 * t125 - t18 * t20 - t16 * t19 - pkin(4) * t236 - g(1) * t212 - t251 * (-pkin(3) * t148 - pkin(4) * t131) + (-t35 * t217 + t4 * t147 + (-t147 * t16 + t245 * t18) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t24, t14, -t174, t138, -t170 * t141 + t155, t141 * t188 + t156, pkin(4) * t21 - t230 * t61, t230 * t18 + (t131 * t251 - t236 - t239 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 - t248, t16 * t63 + t18 * t61 + t13 + t176;];
tau_reg = t1;
