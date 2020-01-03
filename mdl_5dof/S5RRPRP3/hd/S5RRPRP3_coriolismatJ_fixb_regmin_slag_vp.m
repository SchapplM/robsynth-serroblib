% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:16
% EndTime: 2019-12-31 19:51:19
% DurationCPUTime: 1.40s
% Computational Cost: add. (2468->205), mult. (4855->247), div. (0->0), fcn. (4872->6), ass. (0->162)
t197 = qJD(1) + qJD(2);
t156 = sin(pkin(8));
t157 = cos(pkin(8));
t158 = sin(qJ(4));
t239 = cos(qJ(4));
t129 = t239 * t156 + t158 * t157;
t164 = -t156 * t158 + t239 * t157;
t84 = pkin(4) * t129 - qJ(5) * t164;
t250 = t197 * t84;
t124 = t164 ^ 2;
t125 = t129 ^ 2;
t66 = t124 - t125;
t249 = t197 * t66;
t87 = t124 + t125;
t248 = t197 * t87;
t151 = -t157 * pkin(3) - pkin(2);
t160 = cos(qJ(2));
t238 = pkin(1) * t160;
t142 = t151 - t238;
t189 = t151 / 0.2e1 + t142 / 0.2e1;
t247 = t189 * t129;
t175 = -pkin(4) * t164 - t129 * qJ(5);
t74 = t151 + t175;
t67 = t74 - t238;
t193 = t74 / 0.2e1 + t67 / 0.2e1;
t246 = t193 * t129;
t245 = t197 * t164;
t244 = t197 * t129;
t154 = t156 ^ 2;
t155 = t157 ^ 2;
t146 = t154 + t155;
t243 = t197 * t146;
t159 = sin(qJ(2));
t237 = pkin(2) * t159;
t236 = t159 * pkin(1);
t5 = t67 * t84;
t9 = t84 * t74;
t100 = t129 * t238;
t101 = t164 * t238;
t28 = t100 * t129 + t101 * t164;
t85 = t87 * qJD(3);
t235 = t28 * qJD(2) + t85;
t234 = pkin(1) * qJD(1);
t233 = pkin(1) * qJD(2);
t232 = t5 * qJD(1);
t231 = t67 * t164;
t230 = t67 * t129;
t150 = qJ(3) + t236;
t153 = t157 * pkin(7);
t123 = t150 * t157 + t153;
t186 = (-pkin(7) - t150) * t156;
t71 = t123 * t158 - t239 * t186;
t229 = t71 * t129;
t72 = t239 * t123 + t158 * t186;
t228 = t72 * t164;
t227 = t74 * t164;
t226 = t74 * t129;
t143 = qJ(3) * t157 + t153;
t187 = (-pkin(7) - qJ(3)) * t156;
t90 = t143 * t158 - t239 * t187;
t225 = t90 * t129;
t91 = t239 * t143 + t158 * t187;
t224 = t91 * t164;
t223 = t175 * qJD(4) + qJD(5) * t164;
t119 = t129 * qJD(5);
t222 = qJD(4) * t84 - t119;
t18 = t228 + t229;
t221 = qJD(1) * t18;
t57 = t84 * t164;
t19 = -t57 + t230;
t220 = qJD(1) * t19;
t58 = t84 * t129;
t20 = -t58 - t231;
t219 = qJD(1) * t20;
t218 = qJD(1) * t28;
t217 = qJD(4) * t71;
t216 = qJD(4) * t90;
t15 = t100 * t71 + t101 * t72 + t67 * t236;
t215 = t15 * qJD(1);
t68 = t72 * qJD(4);
t102 = t146 * t150;
t73 = (-t237 + (t102 - t236) * t160) * pkin(1);
t214 = t73 * qJD(1);
t86 = t91 * qJD(4);
t194 = t159 * t234;
t103 = t164 * t194;
t195 = t159 * t233;
t105 = t164 * t195;
t213 = -t103 - t105;
t184 = t146 * t160;
t122 = pkin(1) * t184;
t144 = t146 * qJD(3);
t212 = t122 * qJD(2) + t144;
t118 = t125 * qJD(5);
t181 = t129 * t195;
t211 = t118 - t181;
t183 = t239 * t238;
t172 = t183 / 0.2e1;
t196 = t158 * t238;
t179 = -t196 / 0.2e1;
t210 = t156 * t179 + t157 * t172;
t178 = t196 / 0.2e1;
t209 = t156 * t172 + t157 * t178;
t173 = -t183 / 0.2e1;
t208 = t156 * t173 + t157 * t179;
t207 = t156 * t178 + t157 * t173;
t206 = qJD(1) * t102;
t205 = qJD(1) * t122;
t204 = qJD(1) * t142;
t203 = qJD(2) * t151;
t202 = qJD(4) * qJ(5);
t201 = t164 * qJD(3);
t200 = t164 * qJD(4);
t199 = t129 * qJD(3);
t198 = t129 * qJD(4);
t192 = qJD(1) * t230;
t191 = t164 * t204;
t190 = t129 * t204;
t185 = pkin(1) * t197;
t182 = t129 * t194;
t180 = t156 * t194;
t177 = t159 * t185;
t176 = t164 * t244;
t167 = t101 * qJ(5) / 0.2e1 - t100 * pkin(4) / 0.2e1;
t1 = -t193 * t84 + t167;
t174 = -t1 * qJD(1) + t9 * qJD(2);
t25 = t224 + t225;
t152 = t236 / 0.2e1;
t7 = t152 + (-t90 / 0.2e1 - t71 / 0.2e1) * t129 - (t91 / 0.2e1 + t72 / 0.2e1) * t164;
t171 = -qJD(1) * t7 + qJD(2) * t25;
t10 = t208 + t57 - t246;
t21 = -t57 + t226;
t170 = qJD(1) * t10 - qJD(2) * t21;
t11 = t164 * t193 + t210 + t58;
t22 = -t58 - t227;
t169 = qJD(1) * t11 - qJD(2) * t22;
t141 = t146 * qJ(3);
t161 = (qJ(3) + t150) * (t154 / 0.2e1 + t155 / 0.2e1);
t70 = -t236 / 0.2e1 + t161;
t168 = qJD(1) * t70 + qJD(2) * t141;
t166 = -t230 / 0.2e1 - t226 / 0.2e1;
t16 = t209 + t246;
t165 = qJD(1) * t16 + qJD(2) * t226;
t31 = -t164 * t189 + t207;
t163 = qJD(1) * t31 - t164 * t203;
t30 = t208 - t247;
t162 = qJD(1) * t30 - t129 * t203;
t145 = t156 * t195;
t89 = t164 * t119;
t88 = t164 * t198;
t82 = t197 * t125;
t69 = t152 + t161;
t59 = t66 * qJD(4);
t54 = t84 * qJD(3);
t33 = t208 + t247;
t32 = t207 + (t142 + t151) * t164 / 0.2e1;
t17 = t166 + t209;
t13 = -t58 - t227 / 0.2e1 - t231 / 0.2e1 + t210;
t12 = -t57 - t166 + t208;
t8 = t224 / 0.2e1 + t228 / 0.2e1 + t225 / 0.2e1 + t229 / 0.2e1 + t152;
t2 = t9 / 0.2e1 + t5 / 0.2e1 + t167;
t3 = [0, 0, 0, 0, -t195, -t160 * t233, -t157 * t195, t145, t212, qJD(2) * t73 + qJD(3) * t102, t88, t59, 0, 0, 0, t142 * t198 - t105, t142 * t200 + t181, qJD(4) * t19 - t105 + t89, t235, qJD(4) * t20 + t211, qJD(2) * t15 + qJD(3) * t18 + qJD(4) * t5 - t67 * t119; 0, 0, 0, 0, -t177, -t160 * t185, -t157 * t177, t145 + t180, t205 + t212, t214 + t69 * qJD(3) + (qJ(3) * t184 - t237) * t233, t88, t59, 0, 0, 0, qJD(4) * t33 + t213, qJD(4) * t32 + t181 + t182, qJD(4) * t12 + t213 + t89, t218 + t235, qJD(4) * t13 - t182 + t211, t215 + (t100 * t90 + t101 * t91 + t74 * t236) * qJD(2) + t8 * qJD(3) + t2 * qJD(4) + t17 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, t243, qJD(2) * t69 + t206, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, qJD(2) * t8 + t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t249, t200, -t198, 0, qJD(2) * t33 + t190 - t68, qJD(2) * t32 + t191 + t217, qJD(2) * t12 + t220 - t68, t223, qJD(2) * t13 - t217 + t219, t232 + t2 * qJD(2) + (-pkin(4) * t72 - qJ(5) * t71) * qJD(4) + t72 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t200, t82, qJD(2) * t17 - t192 + t68; 0, 0, 0, 0, t194, t160 * t234, t157 * t194, -t180, t144 - t205, qJD(3) * t70 - t214, t88, t59, 0, 0, 0, -qJD(4) * t30 + t103, -qJD(4) * t31 - t182, -qJD(4) * t10 + t103 + t89, t85 - t218, -qJD(4) * t11 + t118 + t182, -qJD(3) * t7 - qJD(4) * t1 - qJD(5) * t16 - t215; 0, 0, 0, 0, 0, 0, 0, 0, t144, t141 * qJD(3), t88, t59, 0, 0, 0, t151 * t198, t151 * t200, qJD(4) * t21 + t89, t85, qJD(4) * t22 + t118, qJD(3) * t25 + qJD(4) * t9 - t119 * t74; 0, 0, 0, 0, 0, 0, 0, 0, t243, t168, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t249, t200, -t198, 0, -t162 - t86, -t163 + t216, -t170 - t86, t223, -t169 - t216, (-pkin(4) * t91 - qJ(5) * t90) * qJD(4) + t91 * qJD(5) + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t200, t82, -t165 + t86; 0, 0, 0, 0, 0, 0, 0, 0, -t243, -qJD(2) * t70 - t206, 0, 0, 0, 0, 0, t198, t200, t198, -t248, -t200, qJD(2) * t7 - t221 + t222; 0, 0, 0, 0, 0, 0, 0, 0, -t243, -t168, 0, 0, 0, 0, 0, t198, t200, t198, -t248, -t200, -t171 + t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t245, t244, 0, -t245, t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, -t249, 0, 0, 0, qJD(2) * t30 - t190 - t199, qJD(2) * t31 - t191 - t201, qJD(2) * t10 - t199 - t220, 0, qJD(2) * t11 + t201 - t219, qJD(2) * t1 - t232 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, -t249, 0, 0, 0, -t199 + t162, -t201 + t163, -t199 + t170, 0, t201 + t169, -t174 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, -t245, -t244, 0, t245, -t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, 0, -t82, qJD(2) * t16 + t192 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, 0, -t82, t199 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
