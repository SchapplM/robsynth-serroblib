% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:16
% EndTime: 2019-03-08 21:49:23
% DurationCPUTime: 2.46s
% Computational Cost: add. (2424->340), mult. (5711->452), div. (0->0), fcn. (3731->8), ass. (0->184)
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t113 = cos(pkin(6));
t199 = qJD(1) * t113;
t116 = sin(qJ(2));
t112 = sin(pkin(6));
t200 = qJD(1) * t112;
t177 = t116 * t200;
t83 = qJD(2) * pkin(8) + t177;
t222 = -t115 * t83 + t118 * t199;
t243 = qJD(4) - t222;
t232 = pkin(4) + pkin(8);
t119 = cos(qJ(2));
t198 = qJD(2) * t112;
t167 = qJD(1) * t198;
t153 = t119 * t167;
t194 = qJD(3) * t113;
t244 = qJD(1) * t194 + t153;
t47 = -qJD(3) * pkin(3) + t243;
t195 = qJD(2) * t118;
t197 = qJD(2) * t115;
t102 = qJD(5) + t197;
t186 = qJD(2) * qJD(3);
t164 = t118 * t186;
t160 = pkin(5) * t164;
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t188 = qJD(5) * t117;
t189 = qJD(5) * t114;
t190 = qJD(3) * t118;
t29 = t244 * t115 + t83 * t190;
t21 = pkin(4) * t164 + t29;
t148 = pkin(9) * t115 - qJ(4) * t118;
t187 = t115 * qJD(4);
t126 = t148 * qJD(3) - t187;
t165 = t115 * t186;
t214 = pkin(3) * t165 + t116 * t167;
t30 = qJD(2) * t126 + t214;
t120 = -pkin(3) - pkin(9);
t203 = pkin(4) * t197 + t243;
t31 = t120 * qJD(3) + t203;
t170 = t119 * t200;
t163 = -t115 * qJ(4) - pkin(2);
t71 = t120 * t118 + t163;
t45 = qJD(2) * t71 - t170;
t162 = -t114 * t30 + t117 * t21 - t45 * t188 - t31 * t189;
t2 = -t160 - t162;
t11 = t114 * t31 + t117 * t45;
t6 = qJ(6) * t102 + t11;
t242 = t102 * t6 - t2;
t109 = qJD(3) * qJ(4);
t56 = t115 * t199 + t118 * t83;
t48 = -t109 - t56;
t208 = t114 * t119;
t192 = qJD(3) * t115;
t104 = pkin(3) * t192;
t58 = t104 + t126;
t80 = t232 * t190;
t89 = t232 * t115;
t241 = -t114 * t80 - t117 * t58 - t89 * t188 + t71 * t189 + (t115 * t208 + t116 * t117) * t200;
t193 = qJD(3) * t114;
t73 = t117 * t195 + t193;
t215 = t73 * t102;
t50 = qJD(5) * t73 - t114 * t165;
t240 = t50 - t215;
t169 = t114 * t195;
t191 = qJD(3) * t117;
t75 = -t169 + t191;
t220 = t102 * t75;
t51 = qJD(3) * t188 - qJD(5) * t169 - t117 * t165;
t239 = -t51 + t220;
t172 = t102 * t188;
t207 = t115 * t117;
t182 = t102 * t207;
t238 = qJD(2) * (t114 * t190 + t182) + qJD(3) * t75 + t172;
t121 = qJD(3) ^ 2;
t131 = -qJ(4) * t190 - t187;
t37 = qJD(2) * t131 + t214;
t63 = t104 + t131;
t237 = (-t63 + t177) * qJD(2) - pkin(8) * t121 - t37;
t175 = t119 * t198;
t122 = qJD(2) ^ 2;
t209 = t112 * t122;
t179 = t116 * t209;
t210 = t112 * t116;
t180 = t115 * t210;
t39 = -qJD(3) * t180 + (t175 + t194) * t118;
t236 = (t118 * t175 + t39) * qJD(3) - t115 * t179;
t155 = t115 * t175;
t67 = t113 * t115 + t118 * t210;
t40 = qJD(3) * t67 + t155;
t235 = (t40 + t155) * qJD(3) + t118 * t179;
t145 = t114 * t89 + t117 * t71;
t206 = t117 * t119;
t234 = -t145 * qJD(5) - t114 * t58 + t117 * t80 + (t114 * t116 - t115 * t206) * t200;
t233 = t75 ^ 2;
t231 = qJ(6) * t190 + qJD(6) * t115 - t241;
t230 = pkin(5) * t190 + t234;
t108 = qJD(3) * qJD(4);
t184 = t244 * t118 - t83 * t192;
t23 = -t108 - t184;
t16 = -pkin(4) * t165 - t23;
t3 = pkin(5) * t51 + qJ(6) * t50 - qJD(6) * t75 + t16;
t228 = t114 * t3;
t227 = t117 * t3;
t44 = pkin(4) * t195 + t56;
t35 = t109 + t44;
t12 = pkin(5) * t73 - qJ(6) * t75 + t35;
t226 = t12 * t75;
t225 = t75 * t73;
t150 = pkin(5) * t117 + qJ(6) * t114;
t136 = -pkin(4) - t150;
t224 = -t150 * qJD(5) + t117 * qJD(6) + t136 * t197 - t243;
t105 = pkin(3) * t197;
t61 = t148 * qJD(2) + t105;
t223 = t114 * t44 + t117 * t61;
t221 = qJD(2) * pkin(2);
t219 = t118 * t75;
t218 = t16 * t114;
t217 = t16 * t117;
t216 = t50 * t117;
t87 = -pkin(3) * t118 + t163;
t213 = qJD(2) * t87;
t212 = t102 * t115;
t211 = t102 * t120;
t205 = t121 * t115;
t204 = t121 * t118;
t10 = -t114 * t45 + t117 * t31;
t202 = qJD(6) - t10;
t90 = t232 * t118;
t110 = t115 ^ 2;
t111 = t118 ^ 2;
t201 = t110 - t111;
t196 = qJD(2) * t116;
t185 = -t114 * t21 - t117 * t30 - t31 * t188;
t183 = t114 * t211;
t181 = t117 * t211;
t178 = t115 * t122 * t118;
t176 = t112 * t196;
t174 = t120 * t190;
t173 = t102 * t189;
t171 = t118 * t188;
t159 = t73 * t170;
t158 = t75 * t170;
t154 = qJ(6) * t164;
t95 = t117 * t164;
t5 = -pkin(5) * t102 + t202;
t151 = t114 * t5 + t117 * t6;
t149 = -pkin(5) * t114 + qJ(6) * t117;
t147 = -t114 * t61 + t117 * t44;
t144 = -t114 * t71 + t117 * t89;
t142 = -qJD(3) * t222 + t184;
t141 = qJD(3) * t56 - t29;
t140 = -qJD(2) * t111 + t212;
t139 = t102 * t114;
t135 = t102 * t11 + t162;
t66 = -t113 * t118 + t180;
t42 = -t112 * t206 + t114 * t66;
t41 = t112 * t208 + t117 * t66;
t134 = t45 * t189 + t185;
t84 = -t170 - t221;
t129 = qJD(3) * (t170 + t84 - t221);
t8 = qJD(5) * t42 + t114 * t176 - t40 * t117;
t128 = -t102 * t8 + t41 * t164 + t39 * t73 + t67 * t51;
t57 = -t170 + t213;
t127 = qJD(3) * (-t170 - t57 - t213);
t125 = -qJD(3) * t73 - t102 * t139 + t95;
t9 = qJD(5) * t41 + t40 * t114 + t117 * t176;
t124 = t102 * t9 + t42 * t164 - t39 * t75 + t50 * t67;
t123 = t115 * t29 - t118 * t23 + (t115 * t48 + t118 * t47) * qJD(3);
t86 = qJ(4) - t149;
t85 = t120 * t95;
t79 = t232 * t192;
t78 = -qJ(4) * t195 + t105;
t54 = t150 * t118 + t90;
t46 = t57 * t197;
t38 = pkin(5) * t75 + qJ(6) * t73;
t33 = -pkin(5) * t115 - t144;
t32 = qJ(6) * t115 + t145;
t15 = (t149 * qJD(5) + qJD(6) * t114) * t118 + (-pkin(8) + t136) * t192;
t14 = -pkin(5) * t195 - t147;
t13 = qJ(6) * t195 + t223;
t1 = qJD(6) * t102 - t134 + t154;
t4 = [0, 0, -t179, -t119 * t209, 0, 0, 0, 0, 0, -t235, -t236 (t115 * t40 + t118 * t39 + (-t115 * t67 + t118 * t66) * qJD(3)) * qJD(2), t235, t236, -t23 * t67 + t29 * t66 - t39 * t48 + t40 * t47 + (-t119 * t37 + t57 * t196) * t112, 0, 0, 0, 0, 0, t128, -t124, t128, t41 * t50 - t42 * t51 - t73 * t9 + t75 * t8, t124, t1 * t42 + t12 * t39 - t2 * t41 + t3 * t67 + t5 * t8 + t6 * t9; 0, 0, 0, 0, 0.2e1 * t115 * t164, -0.2e1 * t201 * t186, t204, -t205, 0, -pkin(8) * t204 + t115 * t129, pkin(8) * t205 + t118 * t129 (-t110 - t111) * t153 + t123, t115 * t127 - t237 * t118, t237 * t115 + t118 * t127, t37 * t87 + t57 * t63 + (-t116 * t57 + (-t115 * t47 + t118 * t48) * t119) * t200 + t123 * pkin(8), -t75 * t171 + (t118 * t50 + t75 * t192) * t114 (-t114 * t73 + t117 * t75) * t192 + (t114 * t51 + t216 + (t114 * t75 + t117 * t73) * qJD(5)) * t118, -t102 * t171 - t50 * t115 + (t140 * t114 + t219) * qJD(3), t118 * t173 - t51 * t115 + (t140 * t117 - t118 * t73) * qJD(3) (t102 + t197) * t190, t90 * t51 - t79 * t73 + (-t35 * t191 + t162) * t115 + t234 * t102 + (-t159 - t35 * t189 + t217 + (t144 * qJD(2) + t10) * qJD(3)) * t118, -t90 * t50 - t79 * t75 + ((qJD(3) * t35 + qJD(5) * t45) * t114 + t185) * t115 + t241 * t102 + (-t158 - t35 * t188 - t218 + (-t145 * qJD(2) - t11) * qJD(3)) * t118, t15 * t73 + t51 * t54 + (-t12 * t191 - t2) * t115 + t230 * t102 + (-t159 - t12 * t189 + t227 + (-qJD(2) * t33 - t5) * qJD(3)) * t118, -t32 * t51 - t33 * t50 - t230 * t75 - t231 * t73 + t151 * t192 + (-t1 * t117 - t114 * t2 + (t114 * t6 - t117 * t5) * qJD(5)) * t118, -t15 * t75 + t50 * t54 + (-t12 * t193 + t1) * t115 + t231 * t102 + (t158 + t12 * t188 + t228 + (qJD(2) * t32 + t6) * qJD(3)) * t118, t1 * t32 + t2 * t33 + t3 * t54 + t231 * t6 - t230 * t5 + (-t118 * t170 + t15) * t12; 0, 0, 0, 0, -t178, t201 * t122, 0, 0, 0, -t84 * t197 + t141, -t84 * t195 - t142, 0, -t78 * t195 - t141 + t46, 0.2e1 * t108 + (t115 * t78 + t118 * t57) * qJD(2) + t142, -pkin(3) * t29 - qJ(4) * t23 - t243 * t48 - t47 * t56 - t57 * t78, -t75 * t139 - t216 (-t51 - t220) * t117 + (t50 + t215) * t114, -t173 + t95 + (-t114 * t212 - t219) * qJD(2), -t172 + (-t182 + (t73 - t193) * t118) * qJD(2), -t102 * t195, t85 + qJ(4) * t51 + t218 - t147 * t102 + t203 * t73 + (t117 * t35 - t183) * qJD(5) + (-t10 * t118 + t35 * t207) * qJD(2), -qJ(4) * t50 + t217 + t223 * t102 + t203 * t75 + (-t114 * t35 - t181) * qJD(5) + (t11 * t118 + (-t115 * t35 - t174) * t114) * qJD(2), t102 * t14 + t228 + t51 * t86 + t85 - t224 * t73 + (t117 * t12 - t183) * qJD(5) + (t118 * t5 + t12 * t207) * qJD(2), t13 * t73 - t14 * t75 + (-t6 * t197 + t120 * t50 + t2 + (-t120 * t73 - t6) * qJD(5)) * t117 + (-t5 * t197 - t120 * t51 - t1 + (t120 * t75 - t5) * qJD(5)) * t114, -t102 * t13 - t227 + t50 * t86 + t224 * t75 + (t114 * t12 + t181) * qJD(5) + (-t118 * t6 + (t115 * t12 + t174) * t114) * qJD(2), -t13 * t6 - t14 * t5 + t3 * t86 - t224 * t12 + (qJD(5) * t151 + t1 * t114 - t117 * t2) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, -t110 * t122 - t121, qJD(3) * t48 + t29 + t46, 0, 0, 0, 0, 0, t125, -t238, t125, t239 * t114 + t240 * t117, t238, -qJD(3) * t12 + t242 * t117 + (t102 * t5 + t1) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, -t73 ^ 2 + t233, -t240, t239, t164, -t35 * t75 + t135, t10 * t102 + t35 * t73 + t134, -t38 * t73 + t135 + 0.2e1 * t160 - t226, pkin(5) * t50 - qJ(6) * t51 + (-t11 + t6) * t75 + (t5 - t202) * t73, 0.2e1 * t154 - t12 * t73 + t38 * t75 + (0.2e1 * qJD(6) - t10) * t102 - t134, -pkin(5) * t2 + qJ(6) * t1 - t11 * t5 - t12 * t38 + t202 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164 + t225, -t240, -t102 ^ 2 - t233, t226 - t242;];
tauc_reg  = t4;
