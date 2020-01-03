% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:09
% EndTime: 2019-12-31 19:04:18
% DurationCPUTime: 2.98s
% Computational Cost: add. (4087->331), mult. (9602->473), div. (0->0), fcn. (6165->8), ass. (0->172)
t139 = cos(qJ(4));
t136 = sin(qJ(4));
t189 = t136 * qJD(3);
t137 = sin(qJ(3));
t197 = qJD(1) * t137;
t101 = t139 * t197 + t189;
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t187 = t139 * qJD(3);
t99 = t136 * t197 - t187;
t155 = t101 * t138 - t135 * t99;
t47 = t101 * t135 + t138 * t99;
t230 = t47 * t155;
t232 = -pkin(8) - pkin(7);
t179 = qJD(4) * t232;
t140 = cos(qJ(3));
t186 = t140 * qJD(1);
t159 = pkin(3) * t137 - pkin(7) * t140;
t105 = t159 * qJD(1);
t123 = sin(pkin(9)) * pkin(1) + pkin(6);
t112 = t123 * qJD(1);
t238 = t140 * qJD(2) - t112 * t137;
t42 = t105 * t136 + t139 * t238;
t247 = -t42 + (pkin(8) * t186 + t179) * t136;
t202 = t139 * t140;
t151 = pkin(4) * t137 - pkin(8) * t202;
t41 = t105 * t139 - t136 * t238;
t246 = qJD(1) * t151 - t139 * t179 + t41;
t184 = qJD(1) * qJD(3);
t245 = qJD(3) * qJD(4) + t140 * t184;
t244 = t155 ^ 2 - t47 ^ 2;
t122 = -qJD(4) + t186;
t119 = -qJD(5) + t122;
t193 = qJD(4) * t137;
t172 = qJD(1) * t193;
t180 = t136 * t245 + t139 * t172;
t190 = qJD(5) * t138;
t191 = qJD(5) * t135;
t65 = t136 * t172 - t139 * t245;
t18 = t101 * t191 + t135 * t180 + t138 * t65 + t190 * t99;
t243 = -t119 * t47 - t18;
t125 = t137 * t184;
t188 = t137 * qJD(2);
t79 = t112 * t140 + t188;
t70 = qJD(3) * pkin(7) + t79;
t124 = -cos(pkin(9)) * pkin(1) - pkin(2);
t96 = -pkin(3) * t140 - pkin(7) * t137 + t124;
t73 = t96 * qJD(1);
t37 = t136 * t73 + t139 * t70;
t71 = t238 * qJD(3);
t108 = t159 * qJD(3);
t95 = qJD(1) * t108;
t16 = -qJD(4) * t37 - t136 * t71 + t139 * t95;
t10 = pkin(4) * t125 + t65 * pkin(8) + t16;
t192 = qJD(4) * t139;
t194 = qJD(4) * t136;
t15 = t136 * t95 + t139 * t71 + t192 * t73 - t194 * t70;
t11 = -pkin(8) * t180 + t15;
t36 = -t136 * t70 + t139 * t73;
t27 = -pkin(8) * t101 + t36;
t24 = -pkin(4) * t122 + t27;
t28 = -pkin(8) * t99 + t37;
t1 = (qJD(5) * t24 + t11) * t138 + t135 * t10 - t28 * t191;
t69 = -qJD(3) * pkin(3) - t238;
t45 = pkin(4) * t99 + t69;
t242 = t45 * t47 - t1;
t104 = t123 * t202;
t52 = t136 * t96 + t104;
t240 = t36 * t122 + t15;
t239 = t37 * t122 - t16;
t237 = t140 * t180;
t174 = t140 * t187;
t177 = t136 * t193;
t236 = t174 - t177;
t182 = qJD(4) + qJD(5);
t217 = t138 * t28;
t6 = t135 * t24 + t217;
t2 = -qJD(5) * t6 + t10 * t138 - t135 * t11;
t235 = -t155 * t45 + t2;
t19 = qJD(5) * t155 - t135 * t65 + t138 * t180;
t234 = -t119 * t155 - t19;
t131 = t137 ^ 2;
t154 = qJD(1) * t131 - t122 * t140;
t233 = -t122 * t177 - t154 * t187;
t231 = pkin(4) * t136;
t114 = t232 * t136;
t115 = t232 * t139;
t62 = t114 * t135 - t115 * t138;
t229 = t62 * qJD(5) + t135 * t247 + t138 * t246;
t61 = t114 * t138 + t115 * t135;
t228 = -t61 * qJD(5) + t135 * t246 - t138 * t247;
t175 = t140 * t189;
t103 = t135 * t139 + t136 * t138;
t55 = t182 * t103;
t31 = t135 * t175 + t137 * t55 - t138 * t174;
t206 = t135 * t136;
t102 = -t138 * t139 + t206;
t83 = t102 * t137;
t227 = t19 * t83 + t31 * t47;
t203 = t137 * t139;
t205 = t136 * t137;
t32 = -t191 * t205 + (t182 * t203 + t175) * t138 + t236 * t135;
t82 = t103 * t137;
t226 = t119 * t32 - t125 * t82;
t160 = t180 * t139;
t225 = -t137 * t160 - t174 * t99;
t224 = -t102 * t186 - t138 * t192 - t139 * t190 + t182 * t206;
t223 = -t103 * t186 + t55;
t222 = t123 * t137 * t189 + t108 * t139;
t221 = t101 * t99;
t220 = t135 * t28;
t219 = t136 * t69;
t218 = t137 * t99;
t216 = t139 * t69;
t72 = t79 * qJD(3);
t213 = t72 * t136;
t212 = t72 * t137;
t211 = t72 * t139;
t210 = t99 * t122;
t209 = t101 * t122;
t208 = t122 * t136;
t207 = t122 * t139;
t204 = t136 * t140;
t141 = qJD(3) ^ 2;
t201 = t141 * t137;
t200 = t141 * t140;
t199 = t99 * qJD(4);
t198 = -t140 ^ 2 + t131;
t113 = qJD(1) * t124;
t196 = qJD(3) * t137;
t195 = qJD(3) * t140;
t142 = qJD(1) ^ 2;
t181 = t137 * t142 * t140;
t178 = t101 * t195;
t176 = t137 * t192;
t169 = t140 * t18 + t155 * t196;
t168 = t101 * t196 + t65 * t140;
t167 = -t65 + t199;
t164 = t122 * t176;
t163 = t101 * t176;
t162 = t140 * t125;
t161 = pkin(4) * t194 - t188 - (qJD(1) * t231 + t112) * t140;
t158 = t155 * t32 - t18 * t82;
t85 = t139 * t96;
t40 = -pkin(8) * t203 + t85 + (-t123 * t136 - pkin(4)) * t140;
t44 = -pkin(8) * t205 + t52;
t20 = -t135 * t44 + t138 * t40;
t21 = t135 * t40 + t138 * t44;
t157 = -t136 * t37 - t139 * t36;
t156 = t136 * t36 - t139 * t37;
t152 = 0.2e1 * qJD(3) * t113;
t150 = t140 * t19 - t196 * t47;
t149 = t154 * t136;
t148 = -t119 * t31 + t125 * t83;
t147 = t175 + t176;
t29 = t96 * t192 + t136 * t108 + (-t137 * t187 - t140 * t194) * t123;
t144 = qJD(4) * t157 - t16 * t136 + t15 * t139;
t143 = t212 + t71 * t140 + (-t137 * t79 - t140 * t238) * qJD(3);
t129 = -pkin(4) * t139 - pkin(3);
t88 = (t123 + t231) * t137;
t56 = pkin(4) * t147 + t123 * t195;
t51 = -t123 * t204 + t85;
t38 = pkin(4) * t180 + t72;
t30 = -qJD(4) * t52 + t222;
t23 = -pkin(8) * t147 + t29;
t22 = t151 * qJD(3) + (-t104 + (pkin(8) * t137 - t96) * t136) * qJD(4) + t222;
t8 = t138 * t27 - t220;
t7 = -t135 * t27 - t217;
t5 = t138 * t24 - t220;
t4 = -qJD(5) * t21 - t135 * t23 + t138 * t22;
t3 = qJD(5) * t20 + t135 * t22 + t138 * t23;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t162, -0.2e1 * t198 * t184, t200, -0.2e1 * t162, -t201, 0, -t123 * t200 + t137 * t152, t123 * t201 + t140 * t152, t143, t143 * t123, t101 * t236 - t203 * t65, -t163 + (-t178 + (t65 + t199) * t137) * t136 + t225, t168 - t233, t99 * t176 + (t137 * t180 + t195 * t99) * t136, t164 + t237 + (-t149 - t218) * qJD(3), (-t122 - t186) * t196, -t30 * t122 + (-t16 + (t123 * t99 + t219) * qJD(3)) * t140 + (t123 * t180 + t213 + t69 * t192 + (qJD(1) * t51 + t36) * qJD(3)) * t137, t29 * t122 + (t15 + (t101 * t123 + t216) * qJD(3)) * t140 + (-t69 * t194 - t123 * t65 + t211 + (-qJD(1) * t52 - t37) * qJD(3)) * t137, -t29 * t99 - t52 * t180 - t30 * t101 + t51 * t65 + t157 * t195 + (qJD(4) * t156 - t15 * t136 - t16 * t139) * t137, t15 * t52 + t16 * t51 + t37 * t29 + t36 * t30 + (t195 * t69 + t212) * t123, -t155 * t31 + t18 * t83, -t158 + t227, -t148 + t169, t19 * t82 + t32 * t47, t150 + t226, (-t119 - t186) * t196, -t4 * t119 - t2 * t140 + t88 * t19 + t45 * t32 + t38 * t82 + t56 * t47 + (qJD(1) * t20 + t5) * t196, t1 * t140 + t3 * t119 - t88 * t18 - t45 * t31 - t38 * t83 + t56 * t155 + (-qJD(1) * t21 - t6) * t196, -t1 * t82 - t155 * t4 + t18 * t20 - t19 * t21 + t2 * t83 - t3 * t47 + t31 * t5 - t32 * t6, t1 * t21 + t2 * t20 + t3 * t6 + t38 * t88 + t4 * t5 + t45 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, -t200, 0, t71 * t137 - t72 * t140 + (-t137 * t238 + t140 * t79) * qJD(3), 0, 0, 0, 0, 0, 0, t164 - t237 + (-t149 + t218) * qJD(3), t168 + t233, t163 + (t137 * t167 + t178) * t136 + t225, (-qJD(3) * t156 - t72) * t140 + (qJD(3) * t69 + t144) * t137, 0, 0, 0, 0, 0, 0, -t150 + t226, t148 + t169, t158 + t227, -t1 * t83 - t140 * t38 + t196 * t45 - t2 * t82 - t31 * t6 - t32 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, t198 * t142, 0, t181, 0, 0, -t113 * t197, -t113 * t186, 0, 0, -t101 * t207 - t65 * t136, (-t65 + t210) * t139 + (-t180 + t209) * t136, -t122 * t192 + (t122 * t202 + (-t101 + t189) * t137) * qJD(1), -t208 * t99 - t160, t122 * t194 + (-t122 * t204 + (t99 + t187) * t137) * qJD(1), t122 * t197, -pkin(3) * t180 - t211 + t41 * t122 - t79 * t99 + (pkin(7) * t207 + t219) * qJD(4) + (-t36 * t137 + (-pkin(7) * t196 - t140 * t69) * t136) * qJD(1), pkin(3) * t65 - t79 * t101 - t42 * t122 + t213 + (-pkin(7) * t208 + t216) * qJD(4) + (-t69 * t202 + (-pkin(7) * t187 + t37) * t137) * qJD(1), t41 * t101 + t42 * t99 + ((qJD(4) * t101 - t180) * pkin(7) + t240) * t139 + (pkin(7) * t167 + t239) * t136, -t72 * pkin(3) + pkin(7) * t144 - t36 * t41 - t37 * t42 - t69 * t79, -t103 * t18 - t155 * t224, t18 * t102 - t103 * t19 - t155 * t223 + t224 * t47, t224 * t119 + (qJD(3) * t103 - t155) * t197, t19 * t102 + t223 * t47, t223 * t119 + (-qJD(3) * t102 + t47) * t197, t119 * t197, t38 * t102 + t129 * t19 + t161 * t47 + t223 * t45 + t229 * t119 + (qJD(3) * t61 - t5) * t197, t38 * t103 - t129 * t18 + t161 * t155 - t224 * t45 - t228 * t119 + (-qJD(3) * t62 + t6) * t197, -t1 * t102 - t2 * t103 + t155 * t229 + t61 * t18 - t62 * t19 - t223 * t6 + t224 * t5 + t228 * t47, t1 * t62 + t38 * t129 + t161 * t45 + t2 * t61 - t228 * t6 - t229 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t101 ^ 2 - t99 ^ 2, -t65 - t210, -t221, -t180 - t209, t125, -t69 * t101 - t239, t69 * t99 - t240, 0, 0, t230, t244, t243, -t230, t234, t125, t7 * t119 + (-t101 * t47 + t119 * t191 + t125 * t138) * pkin(4) + t235, -t8 * t119 + (-t101 * t155 + t119 * t190 - t125 * t135) * pkin(4) + t242, t6 * t155 + t8 * t47 - t5 * t47 + t7 * t155 + (-t135 * t19 + t138 * t18 + (t135 * t155 - t138 * t47) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t1 * t135 - t101 * t45 + t138 * t2 + (-t135 * t5 + t138 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t244, t243, -t230, t234, t125, -t6 * t119 + t235, -t5 * t119 + t242, 0, 0;];
tauc_reg = t9;
