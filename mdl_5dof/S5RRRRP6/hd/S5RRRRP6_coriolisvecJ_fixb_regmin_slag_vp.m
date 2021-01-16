% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:24
% EndTime: 2021-01-16 00:10:33
% DurationCPUTime: 2.20s
% Computational Cost: add. (4049->318), mult. (9897->431), div. (0->0), fcn. (6831->6), ass. (0->179)
t228 = cos(qJ(3));
t179 = t228 * qJD(3);
t168 = pkin(2) * t179;
t141 = cos(qJ(2));
t231 = pkin(7) + pkin(6);
t122 = t231 * t141;
t114 = qJD(1) * t122;
t138 = sin(qJ(3));
t102 = t138 * t114;
t139 = sin(qJ(2));
t121 = t231 * t139;
t112 = qJD(1) * t121;
t75 = -t228 * t112 - t102;
t240 = -t168 + t75;
t180 = qJD(1) * t228;
t200 = t138 * t141;
t101 = -qJD(1) * t200 - t139 * t180;
t140 = cos(qJ(4));
t153 = -t138 * t139 + t228 * t141;
t150 = t153 * qJD(3);
t78 = t153 * qJD(2) + t150;
t145 = t78 * qJD(1);
t189 = qJD(2) + qJD(3);
t169 = t140 * t189;
t137 = sin(qJ(4));
t192 = qJD(4) * t137;
t159 = qJD(4) * t169 + t101 * t192 + t140 * t145;
t82 = -t137 * t101 - t169;
t194 = qJD(1) * t139;
t99 = t138 * t194 - t141 * t180;
t96 = qJD(4) + t99;
t226 = t82 * t96;
t239 = t159 - t226;
t190 = qJD(1) * qJD(2);
t238 = -0.2e1 * t190;
t215 = qJD(2) * pkin(2);
t104 = -t112 + t215;
t183 = qJD(2) * t231;
t166 = qJD(1) * t183;
t105 = t139 * t166;
t106 = t141 * t166;
t193 = qJD(3) * t138;
t34 = t104 * t193 - t138 * t105 + t228 * t106 + t114 * t179;
t152 = t140 * t101 - t137 * t189;
t40 = -t152 * qJD(4) + t137 * t145;
t10 = t40 * pkin(4) + t34;
t191 = qJD(4) * t140;
t229 = t82 * pkin(4);
t177 = qJD(5) + t229;
t72 = t228 * t104 - t102;
t59 = -t189 * pkin(3) - t72;
t38 = t177 + t59;
t237 = t10 * t137 + t38 * t191;
t236 = t10 * t140 - t38 * t192;
t225 = t152 * t96;
t235 = t40 - t225;
t103 = t228 * t114;
t74 = -t138 * t112 + t103;
t165 = pkin(2) * t193 - t74;
t234 = -t228 * t121 - t138 * t122;
t201 = t137 * qJ(5);
t233 = -t140 * qJD(5) + t99 * t201;
t232 = t152 ^ 2;
t110 = t228 * t139 + t200;
t79 = t189 * t110;
t68 = t79 * qJD(1);
t230 = t68 * pkin(4);
t227 = t140 * pkin(4);
t224 = -qJ(5) - pkin(8);
t132 = -t141 * pkin(2) - pkin(1);
t120 = t132 * qJD(1);
t56 = t99 * pkin(3) + t101 * pkin(8) + t120;
t73 = t138 * t104 + t103;
t60 = t189 * pkin(8) + t73;
t23 = -t137 * t60 + t140 * t56;
t15 = qJ(5) * t152 + t23;
t14 = t96 * pkin(4) + t15;
t223 = t14 - t15;
t134 = t140 * qJ(5);
t163 = -t101 * pkin(4) + t99 * t134;
t129 = t138 * pkin(2) + pkin(8);
t196 = -qJ(5) - t129;
t170 = qJD(4) * t196;
t69 = -t101 * pkin(3) + t99 * pkin(8);
t58 = pkin(2) * t194 + t69;
t55 = t140 * t58;
t222 = -t140 * t170 + t163 + t55 + (qJD(5) - t240) * t137;
t174 = qJD(4) * t224;
t175 = -t137 * t72 + t140 * t69;
t221 = t137 * qJD(5) - t140 * t174 + t163 + t175;
t217 = t137 * t58 + t140 * t75;
t220 = -t137 * t170 - t140 * t168 + t217 + t233;
t218 = t137 * t69 + t140 * t72;
t219 = -t137 * t174 + t218 + t233;
t71 = -pkin(3) * t153 - t110 * pkin(8) + t132;
t87 = -t138 * t121 + t228 * t122;
t80 = t140 * t87;
t216 = t137 * t71 + t80;
t213 = t101 * t99;
t212 = t137 * t68;
t211 = t137 * t78;
t210 = t137 * t99;
t209 = t140 * t68;
t208 = t140 * t78;
t207 = t140 * t152;
t172 = t140 * t96;
t206 = t140 * t99;
t205 = t34 * t140;
t204 = t159 * t137;
t203 = t96 * t101;
t185 = pkin(4) * t192;
t88 = pkin(4) * t210;
t202 = t185 + t88 + t165;
t143 = qJD(1) ^ 2;
t199 = t141 * t143;
t142 = qJD(2) ^ 2;
t198 = t142 * t139;
t197 = t142 * t141;
t195 = t139 ^ 2 - t141 ^ 2;
t187 = t139 * t215;
t37 = t79 * pkin(3) - t78 * pkin(8) + t187;
t113 = t139 * t183;
t115 = t141 * t183;
t42 = qJD(3) * t234 - t228 * t113 - t138 * t115;
t188 = t137 * t37 + t140 * t42 + t71 * t191;
t51 = t59 * t192;
t181 = t110 * t191;
t178 = t139 * t190;
t176 = t23 * t101 + t51;
t28 = t68 * pkin(3) + (-pkin(8) * t150 + (t139 * pkin(2) - t153 * pkin(8)) * qJD(2)) * qJD(1);
t33 = t104 * t179 - t228 * t105 - t138 * t106 - t114 * t193;
t173 = -t137 * t28 - t140 * t33 - t56 * t191 + t60 * t192;
t171 = pkin(1) * t238;
t130 = -t228 * pkin(2) - pkin(3);
t24 = t137 * t56 + t140 * t60;
t167 = -t24 * t101 + t34 * t137 + t59 * t191;
t46 = -t88 + t73;
t164 = -t46 + t185;
t162 = -t129 * t68 + t59 * t99;
t16 = -t82 * qJ(5) + t24;
t161 = -t137 * t16 - t14 * t140;
t160 = -qJ(5) * t78 - qJD(5) * t110;
t158 = t14 * t101 - t236;
t157 = t40 * qJ(5) + t173;
t156 = -t96 * t191 - t212;
t155 = t96 * t192 - t209;
t154 = t120 * t101 - t34;
t151 = -t16 * t101 + t38 * t206 + t237;
t26 = t140 * t28;
t149 = -t24 * qJD(4) - t137 * t33 + t26;
t147 = -qJ(5) * t159 + t149;
t1 = qJD(5) * t152 + t147 + t230;
t3 = -t82 * qJD(5) - t157;
t148 = t161 * qJD(4) - t1 * t137 - t14 * t206 + t3 * t140 - t16 * t210;
t146 = t120 * t99 - t33;
t43 = t87 * qJD(3) - t138 * t113 + t228 * t115;
t131 = -pkin(3) - t227;
t119 = t140 * pkin(8) + t134;
t118 = t224 * t137;
t117 = t130 - t227;
t108 = t140 * t129 + t134;
t107 = t196 * t137;
t81 = t82 ^ 2;
t66 = t140 * t71;
t57 = t137 * t110 * pkin(4) - t234;
t48 = t101 ^ 2 - t99 ^ 2;
t45 = (-qJD(1) * t110 - t101) * t189;
t44 = t99 * t189 + t145;
t36 = t140 * t37;
t27 = -t110 * t201 + t216;
t21 = -pkin(4) * t153 - t110 * t134 - t137 * t87 + t66;
t19 = (t181 + t211) * pkin(4) + t43;
t12 = -t101 * t152 + t96 * t172 + t212;
t11 = -t96 ^ 2 * t137 - t82 * t101 + t209;
t8 = -t152 * t172 + t204;
t6 = -qJ(5) * t181 + (-qJD(4) * t87 + t160) * t137 + t188;
t5 = -t235 * t137 + t239 * t140;
t4 = t79 * pkin(4) - t137 * t42 + t36 + t160 * t140 + (-t80 + (qJ(5) * t110 - t71) * t137) * qJD(4);
t2 = [0, 0, 0, 0.2e1 * t141 * t178, t195 * t238, t197, -t198, 0, -pkin(6) * t197 + t139 * t171, pkin(6) * t198 + t141 * t171, -t101 * t78 + t110 * t145, t101 * t79 - t110 * t68 + t145 * t153 - t78 * t99, t78 * t189, -t79 * t189, 0, t132 * t68 + t120 * t79 - t43 * t189 + (-qJD(1) * t153 + t99) * t187, pkin(2) * t110 * t178 - t101 * t187 + t120 * t78 + t132 * t145 - t42 * t189, -t78 * t207 + (t140 * t159 + t152 * t192) * t110, (t137 * t152 - t140 * t82) * t78 + (-t204 - t140 * t40 + (t137 * t82 + t207) * qJD(4)) * t110, -t155 * t110 - t152 * t79 - t153 * t159 + t78 * t172, t110 * t156 + t153 * t40 - t96 * t211 - t82 * t79, -t153 * t68 + t96 * t79, (-t191 * t87 + t36) * t96 + t66 * t68 - (-t191 * t60 + t26) * t153 + t23 * t79 + t43 * t82 - t234 * t40 + t59 * t181 + ((-qJD(4) * t71 - t42) * t96 - t87 * t68 - (-qJD(4) * t56 - t33) * t153 + t34 * t110 + t59 * t78) * t137, -(-t192 * t87 + t188) * t96 - t216 * t68 - t173 * t153 - t24 * t79 - t43 * t152 - t234 * t159 + t59 * t208 + (-t51 + t205) * t110, -t1 * t153 + t110 * t237 + t14 * t79 + t19 * t82 + t21 * t68 + t38 * t211 + t4 * t96 + t57 * t40, t110 * t236 - t152 * t19 + t153 * t3 + t159 * t57 - t16 * t79 + t38 * t208 - t27 * t68 - t6 * t96, -t21 * t159 - t27 * t40 + t4 * t152 - t6 * t82 + t161 * t78 + (-t1 * t140 - t137 * t3 + (t137 * t14 - t140 * t16) * qJD(4)) * t110, t1 * t21 + t10 * t57 + t14 * t4 + t16 * t6 + t38 * t19 + t3 * t27; 0, 0, 0, -t139 * t199, t195 * t143, 0, 0, 0, t143 * pkin(1) * t139, pkin(1) * t199, -t213, t48, t44, t45, 0, t74 * t189 + (-t189 * t193 - t99 * t194) * pkin(2) + t154, t75 * t189 + (t101 * t194 - t189 * t179) * pkin(2) + t146, t8, t5, t12, t11, t203, t130 * t40 - t55 * t96 + t165 * t82 + (-qJD(4) * t129 * t96 - t34) * t140 + (t240 * t96 + t162) * t137 + t176, t130 * t159 + (t129 * t192 + t217) * t96 - t165 * t152 + (-t168 * t96 + t162) * t140 + t167, t107 * t68 + t117 * t40 + t202 * t82 + t210 * t38 - t222 * t96 + t158, -t108 * t68 + t117 * t159 - t152 * t202 + t220 * t96 + t151, -t107 * t159 - t108 * t40 - t152 * t222 + t220 * t82 + t148, t1 * t107 + t10 * t117 + t3 * t108 - t14 * t222 - t16 * t220 + t202 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t48, t44, t45, 0, t73 * t189 + t154, t72 * t189 + t146, t8, t5, t12, t11, t203, -pkin(3) * t40 + pkin(8) * t156 - t175 * t96 + t210 * t59 - t73 * t82 + t176 - t205, -pkin(3) * t159 + pkin(8) * t155 + t152 * t73 + t206 * t59 + t218 * t96 + t167, t118 * t68 + t131 * t40 - t46 * t82 - t221 * t96 + (qJD(4) * t229 + t38 * t99) * t137 + t158, -t119 * t68 + t131 * t159 - t152 * t164 + t219 * t96 + t151, -t118 * t159 - t119 * t40 - t152 * t221 + t219 * t82 + t148, t1 * t118 + t10 * t131 + t3 * t119 - t14 * t221 - t16 * t219 + t164 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 * t82, -t81 + t232, t159 + t226, -t40 - t225, t68, t152 * t59 + t24 * t96 + t149, t23 * t96 + t59 * t82 + t173, 0.2e1 * t230 + t16 * t96 - (-t177 - t38) * t152 + t147, -t232 * pkin(4) + t15 * t96 + (qJD(5) + t38) * t82 + t157, -pkin(4) * t159 - t223 * t82, t223 * t16 + (t152 * t38 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t239, -t81 - t232, -t14 * t152 + t16 * t82 + t10;];
tauc_reg = t2;
