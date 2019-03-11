% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:35
% EndTime: 2019-03-09 03:49:42
% DurationCPUTime: 2.68s
% Computational Cost: add. (3754->303), mult. (9743->385), div. (0->0), fcn. (7662->8), ass. (0->167)
t147 = qJD(3) - qJD(5);
t156 = cos(qJ(6));
t151 = cos(pkin(10));
t233 = cos(qJ(3));
t195 = t233 * t151;
t176 = qJD(1) * t195;
t150 = sin(pkin(10));
t155 = sin(qJ(3));
t210 = t155 * t150;
t194 = qJD(1) * t210;
t110 = -t176 + t194;
t124 = t233 * t150 + t155 * t151;
t112 = t124 * qJD(1);
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t201 = qJD(5) * t157;
t202 = qJD(5) * t154;
t135 = qJD(3) * t176;
t90 = qJD(3) * t194 - t135;
t115 = t124 * qJD(3);
t91 = qJD(1) * t115;
t185 = -t110 * t201 + t112 * t202 - t154 * t91 + t157 * t90;
t199 = qJD(6) * t156;
t153 = sin(qJ(6));
t200 = qJD(6) * t153;
t239 = t154 * t110 + t157 * t112;
t11 = -t147 * t199 - t156 * t185 - t200 * t239;
t172 = t153 * t147 - t156 * t239;
t12 = -qJD(6) * t172 - t153 * t185;
t241 = -t157 * t110 + t154 * t112;
t245 = qJD(6) + t241;
t229 = t172 * t245;
t53 = t156 * t147 + t153 * t239;
t230 = t53 * t245;
t262 = (t11 - t230) * t156 + (-t12 + t229) * t153;
t228 = t172 * t239;
t25 = qJD(5) * t239 - t154 * t90 - t157 * t91;
t220 = t153 * t25;
t255 = t156 * t245;
t252 = t245 * t255 + t220;
t259 = t228 + t252;
t222 = t11 * t153;
t258 = t172 * t255 - t222;
t216 = t239 * t147;
t256 = t25 + t216;
t225 = pkin(7) + qJ(2);
t131 = t225 * t151;
t126 = qJD(1) * t131;
t107 = t155 * t126;
t130 = t225 * t150;
t125 = qJD(1) * t130;
t76 = -t233 * t125 - t107;
t206 = qJD(4) - t76;
t186 = -t156 * t25 + t200 * t245;
t219 = t153 * t245;
t177 = t219 * t241 + t186;
t227 = t53 * t239;
t254 = t177 - t227;
t215 = t241 * t147;
t253 = t185 + t215;
t251 = t239 * pkin(5);
t142 = -t151 * pkin(2) - pkin(1);
t129 = t142 * qJD(1) + qJD(2);
t56 = t110 * pkin(3) - t112 * qJ(4) + t129;
t37 = -t110 * pkin(4) - t56;
t10 = pkin(5) * t241 - pkin(9) * t239 + t37;
t158 = -pkin(3) - pkin(4);
t246 = -t112 * pkin(8) + t206;
t44 = t158 * qJD(3) + t246;
t149 = qJD(3) * qJ(4);
t77 = -t155 * t125 + t233 * t126;
t52 = t110 * pkin(8) + t77;
t47 = t149 + t52;
t18 = t154 * t44 + t157 * t47;
t15 = -t147 * pkin(9) + t18;
t3 = t156 * t10 - t153 * t15;
t250 = t239 * t3;
t4 = t153 * t10 + t156 * t15;
t249 = t239 * t4;
t248 = t245 * t239;
t247 = t239 * t241;
t244 = t239 ^ 2 - t241 ^ 2;
t148 = qJD(3) * qJD(4);
t198 = qJD(1) * qJD(2);
t193 = qJD(2) * t233;
t174 = qJD(1) * t193;
t192 = qJD(3) * t233;
t214 = -t125 * t192 + t151 * t174;
t45 = t148 + (-qJD(3) * t126 - t150 * t198) * t155 + t214;
t31 = t91 * pkin(8) + t45;
t191 = t155 * t198;
t203 = qJD(3) * t155;
t46 = -t125 * t203 + t126 * t192 + t150 * t174 + t151 * t191;
t35 = t90 * pkin(8) + t46;
t1 = t154 * t35 + t157 * t31 + t44 * t201 - t202 * t47;
t243 = t241 * t37 - t1;
t2 = t154 * t31 - t157 * t35 + t47 * t201 + t44 * t202;
t242 = -t239 * t37 - t2;
t240 = -qJD(6) + t245;
t79 = -t155 * t130 + t233 * t131;
t171 = t157 * qJ(4) + t154 * t158;
t128 = -pkin(9) + t171;
t72 = t112 * pkin(3) + t110 * qJ(4);
t48 = -t112 * pkin(4) - t72;
t238 = (-pkin(9) * t241 + qJD(6) * t128 - t251 + t48) * t245 - t2;
t237 = (t245 * pkin(9) + t251) * t245 + t2;
t236 = qJD(3) * (t76 + t107) + t150 * t191 - t214;
t17 = -t154 * t47 + t157 * t44;
t14 = t147 * pkin(5) - t17;
t123 = -t195 + t210;
t168 = t157 * t123 - t154 * t124;
t73 = t123 * pkin(3) - t124 * qJ(4) + t142;
t50 = -t123 * pkin(4) - t73;
t75 = t154 * t123 + t157 * t124;
t16 = -pkin(5) * t168 - t75 * pkin(9) + t50;
t78 = t233 * t130 + t155 * t131;
t60 = -t124 * pkin(8) + t78;
t61 = t123 * pkin(8) + t79;
t22 = t154 * t60 + t157 * t61;
t114 = t150 * t203 - t151 * t192;
t33 = qJD(5) * t168 - t157 * t114 + t154 * t115;
t21 = t154 * t61 - t157 * t60;
t57 = -t130 * t192 + t151 * t193 + (-qJD(2) * t150 - qJD(3) * t131) * t155;
t38 = t115 * pkin(8) + t57;
t58 = t124 * qJD(2) + t79 * qJD(3);
t39 = t114 * pkin(8) + t58;
t7 = -qJD(5) * t21 + t154 * t39 + t157 * t38;
t235 = t14 * t33 - (qJD(6) * t16 + t7) * t245 + (qJD(6) * t10 + t1) * t168 + t2 * t75 - t22 * t25;
t234 = t112 ^ 2;
t232 = t14 * t75;
t231 = t16 * t25;
t226 = t75 * t25;
t170 = -t154 * qJ(4) + t157 * t158;
t224 = -qJD(5) * t170 + t154 * t52 - t246 * t157;
t223 = qJD(5) * t171 + t246 * t154 + t157 * t52;
t217 = t56 * t112;
t213 = t112 * t110;
t208 = t57 * qJD(3);
t207 = t58 * qJD(3);
t204 = t150 ^ 2 + t151 ^ 2;
t42 = t91 * pkin(3) + t90 * qJ(4) - t112 * qJD(4);
t49 = t115 * pkin(3) + t114 * qJ(4) - t124 * qJD(4);
t197 = t75 * t200;
t196 = t245 * t199;
t181 = t204 * qJD(1) ^ 2;
t179 = t147 * t245;
t178 = t147 ^ 2;
t173 = t245 * t33 + t226;
t30 = -t91 * pkin(4) - t42;
t36 = -t115 * pkin(4) - t49;
t167 = 0.2e1 * t204 * t198;
t166 = t77 * qJD(3) - t46;
t163 = -pkin(9) * t25 + (t14 + t17) * t245;
t162 = -t128 * t25 + (-t14 + t224) * t245;
t160 = 0.2e1 * t112 * qJD(3);
t127 = pkin(5) - t170;
t104 = t110 ^ 2;
t71 = t149 + t77;
t70 = t135 + (t110 - t194) * qJD(3);
t69 = -t135 + (t110 + t194) * qJD(3);
t68 = -qJD(3) * pkin(3) + t206;
t34 = qJD(5) * t75 - t154 * t114 - t157 * t115;
t9 = t34 * pkin(5) - t33 * pkin(9) + t36;
t8 = qJD(5) * t22 + t154 * t38 - t157 * t39;
t6 = t25 * pkin(5) + pkin(9) * t185 + t30;
t5 = t156 * t6;
t13 = [0, 0, 0, 0, 0, t167, qJ(2) * t167, -t112 * t114 - t90 * t124, t114 * t110 - t112 * t115 + t90 * t123 - t124 * t91, -t114 * qJD(3), -t115 * qJD(3), 0, t129 * t115 + t142 * t91 - t207, -t129 * t114 - t142 * t90 - t208, t49 * t110 + t56 * t115 + t42 * t123 + t73 * t91 - t207, -t57 * t110 + t58 * t112 - t68 * t114 - t71 * t115 - t45 * t123 + t46 * t124 - t78 * t90 - t79 * t91, -t49 * t112 + t56 * t114 - t42 * t124 + t73 * t90 + t208, t42 * t73 + t45 * t79 + t46 * t78 + t56 * t49 + t71 * t57 + t68 * t58, -t185 * t75 + t239 * t33, -t168 * t185 - t239 * t34 - t241 * t33 - t226, -t33 * t147, t34 * t147, 0, t8 * t147 - t168 * t30 + t241 * t36 + t50 * t25 + t37 * t34, t7 * t147 - t185 * t50 + t239 * t36 + t30 * t75 + t37 * t33, t172 * t197 + (t11 * t75 - t172 * t33) * t156 (t153 * t172 - t156 * t53) * t33 + (-t222 - t12 * t156 + (t153 * t53 + t156 * t172) * qJD(6)) * t75, -t11 * t168 + t156 * t173 - t172 * t34 - t197 * t245, t12 * t168 - t153 * t173 - t196 * t75 - t53 * t34, -t168 * t25 + t245 * t34, t21 * t12 + t3 * t34 - t5 * t168 + t8 * t53 + (t231 + t9 * t245 + (t15 * t168 - t22 * t245 + t232) * qJD(6)) * t156 + t235 * t153, t21 * t11 - t4 * t34 - t8 * t172 + (-(-qJD(6) * t22 + t9) * t245 - t231 + (-qJD(6) * t15 + t6) * t168 - qJD(6) * t232) * t153 + t235 * t156; 0, 0, 0, 0, 0, -t181, -qJ(2) * t181, 0, 0, 0, 0, 0, t160, -t69, t160, -t104 - t234, t69, t71 * t110 - t68 * t112 + t42, 0, 0, 0, 0, 0, -t25 + t216, t185 - t215, 0, 0, 0, 0, 0, t177 + t227, -t228 + t252; 0, 0, 0, 0, 0, 0, 0, t213, -t104 + t234, t70, 0, 0, -t129 * t112 + t166, t129 * t110 + t236, -t72 * t110 + t166 - t217, pkin(3) * t90 - t91 * qJ(4) + (t71 - t77) * t112 + (t68 - t206) * t110, -t56 * t110 + t72 * t112 + 0.2e1 * t148 - t236, -t46 * pkin(3) + t45 * qJ(4) + t206 * t71 - t56 * t72 - t68 * t77, -t247, -t244, t253, t256, 0, t147 * t223 - t241 * t48 - t242, -t147 * t224 - t239 * t48 - t243, t258, -t262, -t259, t254, t248, t127 * t12 + t162 * t153 - t238 * t156 + t223 * t53 + t250, t127 * t11 + t238 * t153 + t162 * t156 - t172 * t223 - t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t70, -qJD(3) ^ 2 - t234, -t71 * qJD(3) + t217 + t46, 0, 0, 0, 0, 0, -t112 * t241 - t154 * t178, -t112 * t239 - t157 * t178, 0, 0, 0, 0, 0, -t112 * t255 + (t153 * t179 - t12) * t157 + (-t147 * t53 - t196 - t220) * t154, t112 * t219 + (t156 * t179 - t11) * t157 + (t147 * t172 + t186) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, t244, -t253, -t256, 0, -t18 * t147 + t242, -t17 * t147 + t243, -t258, t262, t259, -t254, -t248, -pkin(5) * t12 + t163 * t153 - t237 * t156 - t18 * t53 - t250, -pkin(5) * t11 + t237 * t153 + t163 * t156 + t172 * t18 + t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172 * t53, t172 ^ 2 - t53 ^ 2, t11 + t230, -t12 - t229, t25, -t153 * t1 + t14 * t172 + t240 * t4 + t5, -t156 * t1 + t14 * t53 - t153 * t6 + t240 * t3;];
tauc_reg  = t13;
