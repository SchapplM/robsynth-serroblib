% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:42
% EndTime: 2019-03-09 03:09:47
% DurationCPUTime: 2.35s
% Computational Cost: add. (3842->324), mult. (9064->442), div. (0->0), fcn. (6115->8), ass. (0->164)
t157 = cos(qJ(3));
t201 = t157 * qJD(1);
t140 = -qJD(5) + t201;
t150 = sin(pkin(10));
t152 = cos(pkin(10));
t202 = t152 * qJD(3);
t155 = sin(qJ(3));
t208 = qJD(1) * t155;
t116 = -t150 * t208 + t202;
t203 = t150 * qJD(3);
t117 = t152 * t208 + t203;
t154 = sin(qJ(5));
t156 = cos(qJ(5));
t68 = -t156 * t116 + t154 * t117;
t223 = t68 * t140;
t199 = qJD(1) * qJD(3);
t192 = t157 * t199;
t184 = t150 * t192;
t185 = t156 * t157 * t202;
t204 = qJD(5) * t156;
t32 = t154 * (qJD(5) * t117 + t184) - qJD(1) * t185 - t116 * t204;
t249 = t32 - t223;
t248 = t68 ^ 2;
t119 = t154 * t150 - t156 * t152;
t100 = t119 * t155;
t191 = t155 * t199;
t120 = t156 * t150 + t154 * t152;
t110 = t120 * qJD(5);
t195 = t157 * t203;
t56 = t110 * t155 + t154 * t195 - t185;
t247 = -t100 * t191 + t56 * t140;
t205 = qJD(5) * t154;
t239 = -t150 * t205 + t152 * t204;
t246 = t119 * t201 + t239;
t166 = t120 * t157;
t220 = -qJD(1) * t166 + t110;
t174 = t154 * t116 + t156 * t117;
t207 = qJD(3) * t155;
t245 = t157 * t32 + t174 * t207;
t235 = t174 ^ 2;
t214 = t152 * t157;
t170 = pkin(4) * t155 - pkin(8) * t214;
t181 = pkin(3) * t155 - qJ(4) * t157;
t124 = t181 * qJD(1);
t143 = sin(pkin(9)) * pkin(1) + pkin(7);
t132 = t143 * qJD(1);
t241 = t157 * qJD(2) - t155 * t132;
t61 = t152 * t124 - t150 * t241;
t38 = qJD(1) * t170 + t61;
t197 = t150 * t201;
t62 = t150 * t124 + t152 * t241;
t48 = -pkin(8) * t197 + t62;
t232 = pkin(8) + qJ(4);
t129 = t232 * t150;
t130 = t232 * t152;
t78 = -t154 * t129 + t156 * t130;
t244 = qJD(4) * t120 + qJD(5) * t78 - t154 * t48 + t156 * t38;
t173 = -t156 * t129 - t154 * t130;
t243 = -qJD(4) * t119 + qJD(5) * t173 - t154 * t38 - t156 * t48;
t163 = qJD(3) * t166;
t238 = qJD(5) * t174;
t33 = qJD(1) * t163 + t238;
t189 = -t157 * t33 + t68 * t207;
t242 = t140 * t174;
t57 = t239 * t155 + t163;
t99 = t120 * t155;
t226 = t57 * t140 - t99 * t191;
t240 = t189 + t226;
t237 = t245 - t247;
t144 = -cos(pkin(9)) * pkin(1) - pkin(2);
t113 = -t157 * pkin(3) - t155 * qJ(4) + t144;
t102 = t152 * t113;
t215 = t152 * t155;
t54 = -pkin(8) * t215 + t102 + (-t143 * t150 - pkin(4)) * t157;
t217 = t150 * t155;
t123 = t143 * t214;
t74 = t150 * t113 + t123;
t60 = -pkin(8) * t217 + t74;
t175 = t154 * t54 + t156 * t60;
t165 = t170 * qJD(3);
t107 = qJD(3) * t181 - t155 * qJD(4);
t196 = t143 * t207;
t65 = t152 * t107 + t150 * t196;
t46 = t165 + t65;
t134 = t155 * t143;
t216 = t150 * t157;
t96 = t150 * t107;
t53 = t96 + (-pkin(8) * t216 - t152 * t134) * qJD(3);
t236 = -qJD(5) * t175 - t154 * t53 + t156 * t46;
t187 = -qJD(3) * pkin(3) + qJD(4);
t88 = -t241 + t187;
t63 = -t116 * pkin(4) + t88;
t13 = t68 * pkin(5) - qJ(6) * t174 + t63;
t234 = t13 * t174;
t233 = t174 * t68;
t231 = qJ(6) * t208 - t243;
t230 = pkin(5) * t208 + t244;
t229 = t100 * t33 + t56 * t68;
t147 = t155 * qJD(2);
t104 = t157 * t132 + t147;
t79 = pkin(4) * t197 + t104;
t228 = -t220 * pkin(5) + t246 * qJ(6) + t120 * qJD(6) + t79;
t85 = (qJD(4) + t241) * qJD(3);
t94 = t107 * qJD(1);
t37 = t150 * t94 + t152 * t85;
t90 = qJD(3) * qJ(4) + t104;
t93 = t113 * qJD(1);
t42 = t150 * t93 + t152 * t90;
t41 = -t150 * t90 + t152 * t93;
t25 = -pkin(4) * t201 - t117 * pkin(8) + t41;
t30 = t116 * pkin(8) + t42;
t8 = -t154 * t30 + t156 * t25;
t222 = qJD(6) - t8;
t219 = qJD(3) * t173;
t218 = qJD(3) * t78;
t158 = qJD(3) ^ 2;
t213 = t158 * t155;
t212 = t158 * t157;
t206 = qJD(3) * t157;
t95 = qJD(3) * t147 + t132 * t206;
t98 = pkin(4) * t195 + t143 * t206;
t106 = pkin(4) * t217 + t134;
t148 = t155 ^ 2;
t210 = -t157 ^ 2 + t148;
t133 = qJD(1) * t144;
t209 = qJD(1) * t148;
t198 = t143 * t216;
t75 = pkin(4) * t184 + t95;
t145 = -t152 * pkin(4) - pkin(3);
t36 = -t150 * t85 + t152 * t94;
t26 = qJD(1) * t165 + t36;
t31 = -pkin(8) * t184 + t37;
t188 = t154 * t31 - t156 * t26 + t30 * t204 + t25 * t205;
t186 = pkin(5) * t191;
t183 = qJ(6) * t191;
t182 = t174 * t57 - t99 * t32;
t180 = -t36 * t150 + t37 * t152;
t179 = -t150 * t41 + t152 * t42;
t9 = t154 * t25 + t156 * t30;
t176 = -t154 * t60 + t156 * t54;
t171 = 0.2e1 * qJD(3) * t133;
t169 = -t9 * t140 - t188;
t168 = -t154 * t26 - t156 * t31 - t25 * t204 + t205 * t30;
t167 = t154 * t46 + t156 * t53 + t54 * t204 - t205 * t60;
t2 = -t186 + t188;
t5 = t33 * pkin(5) + t32 * qJ(6) - qJD(6) * t174 + t75;
t162 = -qJ(4) * t207 + (t187 - t88) * t157;
t160 = t33 - t242;
t159 = qJD(1) ^ 2;
t73 = t102 - t198;
t66 = -t152 * t196 + t96;
t64 = t119 * pkin(5) - t120 * qJ(6) + t145;
t34 = t99 * pkin(5) + t100 * qJ(6) + t106;
t24 = pkin(5) * t174 + t68 * qJ(6);
t16 = t157 * pkin(5) - t176;
t15 = -t157 * qJ(6) + t175;
t14 = -t32 - t223;
t10 = t57 * pkin(5) + t56 * qJ(6) + t100 * qJD(6) + t98;
t7 = -t140 * qJ(6) + t9;
t6 = t140 * pkin(5) + t222;
t4 = -pkin(5) * t207 - t236;
t3 = qJ(6) * t207 - t157 * qJD(6) + t167;
t1 = -t140 * qJD(6) - t168 + t183;
t11 = [0, 0, 0, 0, 0.2e1 * t157 * t191, -0.2e1 * t210 * t199, t212, -t213, 0, -t143 * t212 + t155 * t171, t143 * t213 + t157 * t171, t95 * t217 + (-qJD(1) * t65 - t36) * t157 + ((-t116 * t143 + t150 * t88) * t157 + (t41 + (t73 + t198) * qJD(1)) * t155) * qJD(3), t95 * t215 + (qJD(1) * t66 + t37) * t157 + ((t117 * t143 + t152 * t88) * t157 + (-t42 + (-t74 + t123) * qJD(1)) * t155) * qJD(3), t66 * t116 - t65 * t117 + (-t150 * t37 - t152 * t36) * t155 + (-t150 * t42 - t152 * t41 + (-t150 * t74 - t152 * t73) * qJD(1)) * t206, t36 * t73 + t37 * t74 + t41 * t65 + t42 * t66 + (t155 * t95 + t206 * t88) * t143, t32 * t100 - t174 * t56, -t182 + t229, t245 + t247, t226 - t189 (-t140 - t201) * t207, -t236 * t140 + t188 * t157 + t98 * t68 + t106 * t33 + t75 * t99 + t63 * t57 + (qJD(1) * t176 + t8) * t207, t167 * t140 - t168 * t157 + t98 * t174 - t106 * t32 - t75 * t100 - t63 * t56 + (-t175 * qJD(1) - t9) * t207, t10 * t68 + t13 * t57 + t4 * t140 + t2 * t157 + t34 * t33 + t5 * t99 + (-qJD(1) * t16 - t6) * t207, -t1 * t99 - t2 * t100 - t15 * t33 - t16 * t32 + t174 * t4 - t3 * t68 - t6 * t56 - t7 * t57, -t1 * t157 - t10 * t174 + t5 * t100 + t13 * t56 - t3 * t140 + t34 * t32 + (qJD(1) * t15 + t7) * t207, t1 * t15 + t13 * t10 + t2 * t16 + t7 * t3 + t5 * t34 + t6 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, -t212 (-t116 * t155 - t150 * t209) * qJD(3) (t117 * t155 - t152 * t209) * qJD(3) (t116 * t152 + t117 * t150) * t206, -t95 * t157 + t180 * t155 + (t155 * t88 + t157 * t179) * qJD(3), 0, 0, 0, 0, 0, t240, t237, t240, t182 + t229, -t237, -t1 * t100 + t13 * t207 - t5 * t157 + t2 * t99 - t7 * t56 + t6 * t57; 0, 0, 0, 0, -t155 * t159 * t157, t210 * t159, 0, 0, 0, t104 * qJD(3) - t133 * t208 - t95, -t133 * t201, t104 * t116 - t95 * t152 + (t150 * t162 - t155 * t41 + t157 * t61) * qJD(1), -t104 * t117 + t95 * t150 + (t152 * t162 + t155 * t42 - t157 * t62) * qJD(1), -t62 * t116 + t61 * t117 + (qJD(4) * t116 + t201 * t41 + t37) * t152 + (qJD(4) * t117 + t201 * t42 - t36) * t150, -t95 * pkin(3) + qJ(4) * t180 + qJD(4) * t179 - t88 * t104 - t41 * t61 - t42 * t62, -t32 * t120 + t174 * t246, t32 * t119 - t120 * t33 - t174 * t220 - t246 * t68, -t246 * t140 + (qJD(3) * t120 - t174) * t208, t220 * t140 + (-qJD(3) * t119 + t68) * t208, t140 * t208, t75 * t119 + t145 * t33 - t79 * t68 + t220 * t63 + t244 * t140 + (-t8 + t219) * t208, t75 * t120 - t145 * t32 - t79 * t174 + t246 * t63 + t243 * t140 + (t9 - t218) * t208, t5 * t119 + t64 * t33 - t228 * t68 + t230 * t140 + t220 * t13 + (t6 + t219) * t208, -t1 * t119 + t2 * t120 + t173 * t32 + t174 * t230 - t220 * t7 + t231 * t68 + t246 * t6 - t78 * t33, -t5 * t120 + t64 * t32 + t228 * t174 + t231 * t140 - t246 * t13 + (-t7 + t218) * t208, t1 * t78 - t13 * t228 - t173 * t2 + t230 * t6 - t231 * t7 + t5 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t117 + t203) * t201 (-t116 + t202) * t201, -t116 ^ 2 - t117 ^ 2, -t42 * t116 + t41 * t117 + t95, 0, 0, 0, 0, 0, t160, -t249, t160, -t235 - t248, t249, -t174 * t6 + t68 * t7 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t235 - t248, t14, -t120 * t192 - t238 - t242, t191, -t174 * t63 + t169, -t8 * t140 + t63 * t68 + t168, -t24 * t68 + t169 + 0.2e1 * t186 - t234, pkin(5) * t32 - t33 * qJ(6) + (t7 - t9) * t174 + (t6 - t222) * t68, 0.2e1 * t183 - t13 * t68 + t24 * t174 + (-0.2e1 * qJD(6) + t8) * t140 - t168, -t2 * pkin(5) + t1 * qJ(6) - t13 * t24 + t222 * t7 - t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191 + t233, t14, -t140 ^ 2 - t235, t7 * t140 + t2 + t234;];
tauc_reg  = t11;
