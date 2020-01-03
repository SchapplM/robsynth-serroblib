% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:28
% EndTime: 2020-01-03 11:50:33
% DurationCPUTime: 1.68s
% Computational Cost: add. (3119->221), mult. (6722->327), div. (0->0), fcn. (6871->6), ass. (0->194)
t162 = cos(qJ(4));
t250 = t162 * pkin(3);
t153 = pkin(4) + t250;
t257 = -t153 / 0.2e1;
t158 = sin(pkin(8));
t163 = cos(qJ(3));
t221 = t163 * t158;
t160 = sin(qJ(4));
t161 = sin(qJ(3));
t226 = t160 * t161;
t125 = -t158 * t226 + t162 * t221;
t120 = t125 * qJ(5);
t159 = cos(pkin(8));
t139 = -t159 * pkin(2) - t158 * pkin(6) - pkin(1);
t114 = -t163 * t159 * qJ(2) - t161 * t139;
t229 = t158 * t161;
t103 = -pkin(7) * t229 - t114;
t227 = t160 * t103;
t133 = t163 * t139;
t194 = pkin(7) * t221;
t237 = qJ(2) * t161;
t92 = -t194 + t133 + (-pkin(3) - t237) * t159;
t79 = t162 * t92;
t49 = -t79 + t227;
t39 = -t120 - t49;
t264 = t125 ^ 2;
t263 = -t79 / 0.2e1;
t262 = -t92 / 0.2e1;
t223 = t162 * t161;
t225 = t160 * t163;
t138 = t223 + t225;
t122 = t138 * t158;
t233 = t122 * qJ(5);
t100 = t162 * t103;
t174 = t159 * t237 - t133;
t102 = -t174 - t194;
t228 = t160 * t102;
t55 = -t100 - t228;
t42 = t55 + t233;
t261 = t42 * pkin(4);
t183 = -t100 / 0.2e1;
t222 = t162 * t163;
t137 = t222 - t226;
t260 = t137 / 0.2e1;
t259 = -t138 / 0.2e1;
t258 = t138 / 0.2e1;
t256 = t158 / 0.2e1;
t255 = t159 / 0.2e1;
t254 = pkin(3) * t160;
t253 = pkin(4) * t122;
t252 = pkin(4) * t125;
t124 = t138 * t159;
t251 = t124 * pkin(4);
t50 = -t160 * t92 - t100;
t40 = -t50 - t233;
t242 = t40 * t137;
t38 = -t159 * pkin(4) + t39;
t248 = t38 * t259 + t242 / 0.2e1;
t247 = t38 - t39;
t246 = pkin(3) * qJD(3);
t245 = pkin(3) * qJD(4);
t244 = pkin(4) * qJD(4);
t127 = t137 * t159;
t167 = t127 * t254 / 0.2e1 + t124 * t257;
t224 = t162 * t102;
t56 = t224 - t227;
t43 = -t120 + t56;
t3 = (-t43 / 0.2e1 + t38 / 0.2e1) * t138 + (-t40 / 0.2e1 - t42 / 0.2e1) * t137 + t167;
t243 = t3 * qJD(1);
t5 = t247 * t122;
t241 = t5 * qJD(1);
t166 = t39 * t258 - t242 / 0.2e1 + t248;
t7 = t251 / 0.2e1 + t166;
t240 = t7 * qJD(1);
t8 = (-t40 - t42) * t125 + (t38 - t43) * t122;
t239 = t8 * qJD(1);
t206 = pkin(3) * t229 + t158 * qJ(2);
t99 = t206 + t253;
t9 = -t247 * t40 + t99 * t252;
t238 = t9 * qJD(1);
t195 = pkin(3) * t221;
t10 = t40 * t43 + t38 * t42 + t99 * (t195 + t252);
t236 = t10 * qJD(1);
t11 = -t38 * t124 + t40 * t127 + t99 * t158;
t235 = t11 * qJD(1);
t12 = -t40 * t122 - t38 * t125;
t234 = t12 * qJD(1);
t191 = pkin(3) * t255;
t176 = t191 + t102 / 0.2e1;
t13 = t183 + t100 / 0.2e1 + (t262 + t176) * t160;
t232 = t13 * qJD(1);
t15 = t176 * t162 + t263;
t231 = t15 * qJD(1);
t156 = t158 ^ 2;
t230 = t156 * t163;
t101 = t206 * t125;
t24 = -t122 * t195 + t55 * t159 - t101;
t220 = t24 * qJD(1);
t181 = t206 * t122;
t25 = t125 * t195 + t56 * t159 - t181;
t219 = t25 * qJD(1);
t27 = -t49 * t159 - t181;
t218 = t27 * qJD(1);
t28 = t50 * t159 - t101;
t217 = t28 * qJD(1);
t184 = -t160 * t122 / 0.2e1;
t45 = (t257 - pkin(4) / 0.2e1) * t125 + (t184 - t221 / 0.2e1) * pkin(3);
t216 = t45 * qJD(1);
t46 = -t127 * t122 + t124 * t125;
t215 = t46 * qJD(1);
t190 = t250 / 0.2e1;
t175 = t190 + t257;
t171 = pkin(4) / 0.2e1 + t175;
t48 = t171 * t122;
t214 = t48 * qJD(1);
t170 = -t122 * t259 + t125 * t260;
t51 = t256 + t170;
t213 = t51 * qJD(1);
t121 = t122 ^ 2;
t61 = t121 - t264;
t212 = t61 * qJD(1);
t67 = -t158 * t122 - t124 * t159;
t211 = t67 * qJD(1);
t68 = t158 * t125 + t127 * t159;
t210 = t68 * qJD(1);
t69 = t121 + t264;
t209 = t69 * qJD(1);
t169 = -t225 / 0.2e1 - t223 / 0.2e1;
t80 = (t259 + t169) * t159;
t72 = t80 * qJD(1);
t168 = -t222 / 0.2e1 + t226 / 0.2e1;
t81 = (-t137 / 0.2e1 + t168) * t159;
t74 = t81 * qJD(1);
t84 = t156 * t237 + t174 * t159;
t208 = t84 * qJD(1);
t85 = qJ(2) * t230 - t114 * t159;
t207 = t85 * qJD(1);
t149 = t159 ^ 2 + t156;
t205 = qJD(1) * t125;
t204 = qJD(1) * t159;
t203 = qJD(2) * t159;
t202 = qJD(3) * t159;
t134 = t149 * t161;
t201 = t134 * qJD(1);
t135 = t149 * t163;
t200 = t135 * qJD(1);
t136 = (t161 ^ 2 - t163 ^ 2) * t156;
t199 = t136 * qJD(1);
t140 = t149 * qJ(2);
t198 = t140 * qJD(1);
t197 = t149 * qJD(1);
t196 = qJD(3) + qJD(4);
t193 = pkin(4) * t205;
t192 = t160 * t245;
t189 = t161 * t230;
t188 = t161 * t202;
t187 = t163 * t202;
t186 = t161 * t204;
t185 = t163 * t204;
t182 = pkin(3) * t196;
t180 = t196 * t125;
t179 = t196 * t159;
t178 = -qJD(3) + t204;
t177 = qJD(1) * t189;
t173 = t178 * t161;
t172 = t178 * t163;
t131 = (-t153 + t250) * t254;
t164 = (t40 * t162 / 0.2e1 + (t39 / 0.2e1 - t38 / 0.2e1) * t160) * pkin(3) + t40 * t257;
t2 = -t261 / 0.2e1 + t164;
t65 = t171 * t138;
t165 = -t2 * qJD(1) - t65 * qJD(2) - t131 * qJD(3);
t113 = t125 * t204;
t112 = t122 * t204;
t83 = t138 * t255 + t169 * t159;
t82 = t137 * t255 + t168 * t159;
t76 = t83 * qJD(2);
t75 = t82 * qJD(2);
t73 = t81 * qJD(2);
t71 = t80 * qJD(2);
t70 = t122 * t205;
t64 = pkin(4) * t259 + t175 * t138;
t58 = t113 - t180;
t57 = -t196 * t122 + t112;
t54 = -t196 * t138 - t72;
t53 = -t196 * t137 - t74;
t52 = t256 - t170;
t47 = t253 / 0.2e1 - t175 * t122;
t44 = pkin(3) * t184 + t125 * t257 + t195 / 0.2e1 + t252 / 0.2e1;
t16 = t159 * t190 + t227 + t263 - t224 / 0.2e1;
t14 = 0.2e1 * t183 - t228 / 0.2e1 + (t191 + t262) * t160;
t6 = -t251 / 0.2e1 + t166;
t4 = t43 * t258 + t42 * t260 + t167 + t248;
t1 = t261 / 0.2e1 + t164;
t17 = [0, 0, 0, 0, 0, t149 * qJD(2), t140 * qJD(2), -qJD(3) * t189, t136 * qJD(3), t158 * t188, t158 * t187, 0, t134 * qJD(2) + t85 * qJD(3), t135 * qJD(2) - t84 * qJD(3), -t122 * t180, t196 * t61, t122 * t179, t125 * t179, 0, -t67 * qJD(2) - t24 * qJD(3) - t28 * qJD(4), t68 * qJD(2) + t25 * qJD(3) + t27 * qJD(4), t46 * qJD(2) + t8 * qJD(3) + t5 * qJD(4) + t69 * qJD(5), t11 * qJD(2) + t10 * qJD(3) + t9 * qJD(4) + t12 * qJD(5); 0, 0, 0, 0, 0, t197, t198, 0, 0, 0, 0, 0, t201, t200, 0, 0, 0, 0, 0, t196 * t83 - t211, t196 * t82 + t210, t215, t235 + (-t124 * t137 + t127 * t138) * qJD(2) + t4 * qJD(3) + t6 * qJD(4) + t52 * qJD(5); 0, 0, 0, 0, 0, 0, 0, -t177, t199, t158 * t173, t158 * t172, 0, t114 * qJD(3) + t207, qJD(3) * t174 - t208, -t70, t212, t57, t58, 0, t55 * qJD(3) + t14 * qJD(4) - t220 + t76, -t56 * qJD(3) + t16 * qJD(4) + t219 + t75, t239 + (t153 * t122 - t125 * t254) * qJD(3) + t47 * qJD(4), t236 + t4 * qJD(2) + (t42 * t153 + t254 * t43) * qJD(3) + t1 * qJD(4) + t44 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t212, t57, t58, 0, t14 * qJD(3) + t50 * qJD(4) - t217 + t76, t16 * qJD(3) + t49 * qJD(4) + t218 + t75, t47 * qJD(3) + t122 * t244 + t241, t6 * qJD(2) + t1 * qJD(3) - t244 * t40 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, t52 * qJD(2) + t44 * qJD(3) + t234; 0, 0, 0, 0, 0, -t197, -t198, 0, 0, 0, 0, 0, t188 - t201, t187 - t200, 0, 0, 0, 0, 0, -t196 * t80 + t211, -t196 * t81 - t210, -t215, -t3 * qJD(3) + t7 * qJD(4) - t51 * qJD(5) - t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t172, 0, 0, 0, 0, 0, t54, t53, 0, -t243 + (t137 * t254 - t138 * t153) * qJD(3) + t64 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, 0, t64 * qJD(3) - t138 * t244 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213; 0, 0, 0, 0, 0, 0, 0, t177, -t199, -t158 * t186, -t158 * t185, 0, -t161 * t203 - t207, -t163 * t203 + t208, t70, -t212, -t112, -t113, 0, t13 * qJD(4) + t220 + t71, t15 * qJD(4) - t219 + t73, -t48 * qJD(4) - t239, t3 * qJD(2) + t2 * qJD(4) + t45 * qJD(5) - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186, -t185, 0, 0, 0, 0, 0, t72, t74, 0, t65 * qJD(4) + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, -t162 * t245, 0, t131 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t182 + t232, -t162 * t182 + t231, -t214, -pkin(4) * t192 - t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t212, -t112, -t113, 0, -t13 * qJD(3) + t217 + t71, -t15 * qJD(3) - t218 + t73, t48 * qJD(3) - t241, -t7 * qJD(2) - t2 * qJD(3) - qJD(5) * t252 - t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t74, 0, -t65 * qJD(3) - t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160 * t246 - t232, t162 * t246 - t231, t214, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t51 * qJD(2) - t45 * qJD(3) + t125 * t244 - t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t17;
