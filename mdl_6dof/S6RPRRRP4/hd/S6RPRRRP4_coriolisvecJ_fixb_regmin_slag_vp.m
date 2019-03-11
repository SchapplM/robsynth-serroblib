% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:45
% EndTime: 2019-03-09 06:08:54
% DurationCPUTime: 3.15s
% Computational Cost: add. (7157->337), mult. (18881->438), div. (0->0), fcn. (14926->8), ass. (0->190)
t150 = sin(qJ(5));
t153 = cos(qJ(5));
t154 = cos(qJ(4));
t148 = sin(pkin(10));
t149 = cos(pkin(10));
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t129 = t155 * t148 + t152 * t149;
t124 = t129 * qJD(3);
t116 = qJD(1) * t124;
t208 = t155 * t149;
t209 = t152 * t148;
t170 = -t208 + t209;
t121 = t170 * qJD(1);
t237 = pkin(7) + qJ(2);
t133 = t237 * t148;
t130 = qJD(1) * t133;
t134 = t237 * t149;
t131 = qJD(1) * t134;
t251 = -t155 * t130 - t152 * t131;
t72 = -t116 * pkin(8) - qJD(2) * t121 + t251 * qJD(3);
t122 = t129 * qJD(1);
t86 = -t122 * pkin(8) + t251;
t84 = qJD(3) * pkin(3) + t86;
t185 = qJD(4) * t84 + t72;
t151 = sin(qJ(4));
t204 = qJD(4) * t151;
t200 = t155 * qJD(3);
t137 = t149 * qJD(1) * t200;
t189 = qJD(1) * t209;
t115 = -qJD(3) * t189 + t137;
t162 = t129 * qJD(2);
t161 = qJD(1) * t162;
t172 = t152 * t130 - t155 * t131;
t73 = -t115 * pkin(8) + t172 * qJD(3) - t161;
t87 = -t121 * pkin(8) - t172;
t187 = t151 * t73 - t87 * t204;
t14 = t185 * t154 + t187;
t197 = qJD(3) + qJD(4);
t83 = t154 * t87;
t46 = t151 * t84 + t83;
t40 = t197 * pkin(9) + t46;
t141 = -t149 * pkin(2) - pkin(1);
t132 = t141 * qJD(1) + qJD(2);
t106 = t121 * pkin(3) + t132;
t174 = -t151 * t121 + t154 * t122;
t98 = t154 * t121 + t151 * t122;
t47 = t98 * pkin(4) - pkin(9) * t174 + t106;
t20 = t150 * t47 + t153 * t40;
t203 = qJD(4) * t154;
t163 = -t154 * t115 + t151 * t116 + t121 * t203 + t122 * t204;
t182 = t151 * t115 + t154 * t116;
t249 = t174 * qJD(4);
t63 = t182 + t249;
t27 = t116 * pkin(3) + t63 * pkin(4) + t163 * pkin(9);
t26 = t153 * t27;
t158 = -t20 * qJD(5) - t150 * t14 + t26;
t179 = t153 * t197;
t202 = qJD(5) * t150;
t35 = -qJD(5) * t179 + t153 * t163 + t174 * t202;
t90 = t150 * t197 + t153 * t174;
t1 = t63 * pkin(5) + t35 * qJ(6) - t90 * qJD(6) + t158;
t88 = t150 * t174 - t179;
t10 = -t88 * qJ(6) + t20;
t260 = qJD(5) + t98;
t248 = t10 * t260 + t1;
t271 = t153 * t260;
t201 = qJD(5) * t153;
t164 = t153 * t14 + t150 * t27 + t47 * t201 - t40 * t202;
t36 = qJD(5) * t90 - t150 * t163;
t3 = -t36 * qJ(6) - t88 * qJD(6) + t164;
t19 = -t150 * t40 + t153 * t47;
t9 = -t90 * qJ(6) + t19;
t7 = pkin(5) * t260 + t9;
t273 = -t248 * t150 + t3 * t153 - t271 * t7;
t224 = t150 * t260;
t272 = pkin(5) * t224;
t252 = t98 * t197;
t269 = -t163 + t252;
t263 = t153 * t98;
t33 = t35 * t150;
t268 = -t33 + (t201 + t263) * t90;
t264 = t150 * t98;
t267 = -qJ(6) * t264 + t153 * qJD(6);
t60 = t150 * t63;
t229 = t201 * t260 + t60;
t239 = t90 * t174;
t266 = t260 * t263 + t229 - t239;
t265 = t106 * t98;
t262 = t174 * t98;
t194 = pkin(3) * t203;
t82 = t151 * t87;
t49 = t154 * t86 - t82;
t261 = t49 - t194;
t258 = t174 ^ 2 - t98 ^ 2;
t67 = pkin(4) * t174 + t98 * pkin(9);
t145 = t153 * qJ(6);
t257 = -t174 * pkin(5) - t98 * t145;
t238 = t174 * t88;
t253 = t260 * t174;
t62 = t153 * t63;
t165 = t202 * t260 - t62;
t94 = -t129 * pkin(8) - t155 * t133 - t152 * t134;
t171 = t152 * t133 - t155 * t134;
t95 = -pkin(8) * t170 - t171;
t57 = t151 * t95 - t154 * t94;
t207 = t174 * qJD(3);
t250 = t207 - t182;
t45 = t154 * t84 - t82;
t39 = -t197 * pkin(4) - t45;
t37 = t39 * t202;
t247 = -t19 * t174 + t37;
t15 = t151 * t72 - t154 * t73 + t87 * t203 + t84 * t204;
t246 = -t106 * t174 - t15;
t245 = t15 * t150 + t20 * t174 + t39 * t201;
t244 = t90 ^ 2;
t243 = t7 - t9;
t242 = t122 * pkin(3);
t241 = t154 * pkin(3);
t123 = t170 * qJD(3);
t173 = -t151 * t129 - t154 * t170;
t70 = t173 * qJD(4) - t154 * t123 - t151 * t124;
t240 = t39 * t70;
t236 = -qJ(6) - pkin(9);
t142 = t151 * pkin(3) + pkin(9);
t206 = -qJ(6) - t142;
t181 = qJD(5) * t206;
t52 = t242 + t67;
t51 = t153 * t52;
t235 = t153 * t181 - t51 + (-qJD(6) + t261) * t150 + t257;
t233 = -t150 * t36 - t88 * t201;
t232 = t150 * t67 + t153 * t45;
t231 = t150 * t52 + t153 * t49;
t58 = t151 * t94 + t154 * t95;
t54 = t153 * t58;
t103 = t154 * t129 - t151 * t170;
t110 = pkin(3) * t170 + t141;
t59 = -pkin(4) * t173 - t103 * pkin(9) + t110;
t230 = t150 * t59 + t54;
t227 = pkin(3) * qJD(4);
t226 = t15 * t153;
t225 = t150 * t90;
t34 = t153 * t35;
t222 = t153 * t88;
t221 = t153 * t90;
t217 = t150 * t181 + t153 * t194 - t231 + t267;
t186 = qJD(5) * t236;
t216 = t150 * t186 - t232 + t267;
t188 = -t150 * t45 + t153 * t67;
t215 = -t150 * qJD(6) + t153 * t186 - t188 + t257;
t211 = t103 * t150;
t205 = t148 ^ 2 + t149 ^ 2;
t198 = qJD(1) * qJD(2);
t159 = -t133 * t200 + qJD(2) * t208 + (-qJD(2) * t148 - qJD(3) * t134) * t152;
t74 = -t124 * pkin(8) + t159;
t157 = t171 * qJD(3) - t162;
t75 = t123 * pkin(8) + t157;
t23 = -t57 * qJD(4) + t151 * t75 + t154 * t74;
t71 = t103 * qJD(4) - t151 * t123 + t154 * t124;
t30 = t124 * pkin(3) + t71 * pkin(4) - t70 * pkin(9);
t196 = t150 * t30 + t153 * t23 + t59 * t201;
t195 = -t34 + t233;
t192 = t90 * t202;
t191 = -t153 * pkin(5) - pkin(4);
t190 = t103 * t201;
t183 = t205 * qJD(1) ^ 2;
t48 = t151 * t86 + t83;
t178 = pkin(3) * t204 - t48;
t177 = -t142 * t63 + t39 * t98;
t176 = t222 + t225;
t175 = -qJ(6) * t70 - qJD(6) * t103;
t169 = -t260 * t264 - t165;
t168 = 0.2e1 * t205 * t198;
t6 = t36 * pkin(5) + t15;
t167 = -t192 - t34;
t24 = t58 * qJD(4) + t151 * t74 - t154 * t75;
t143 = -pkin(4) - t241;
t136 = t153 * pkin(9) + t145;
t135 = t236 * t150;
t127 = t153 * t142 + t145;
t126 = t206 * t150;
t85 = t88 ^ 2;
t56 = t153 * t59;
t31 = t88 * pkin(5) + qJD(6) + t39;
t29 = t153 * t30;
t22 = -qJ(6) * t211 + t230;
t16 = -pkin(5) * t173 - t103 * t145 - t150 * t58 + t56;
t5 = -qJ(6) * t190 + (-qJD(5) * t58 + t175) * t150 + t196;
t4 = t71 * pkin(5) - t150 * t23 + t29 + t175 * t153 + (-t54 + (qJ(6) * t103 - t59) * t150) * qJD(5);
t2 = [0, 0, 0, 0, 0, t168, qJ(2) * t168, t115 * t129 - t122 * t123, -t115 * t170 - t116 * t129 + t121 * t123 - t124 * t122, -qJD(3) * t123, -t124 * qJD(3), 0, t157 * qJD(3) + t141 * t116 + t132 * t124, -t159 * qJD(3) + t141 * t115 - t132 * t123, -t163 * t103 + t174 * t70, -t103 * t63 - t163 * t173 - t174 * t71 - t70 * t98, t70 * t197, -t71 * t197, 0, t110 * t63 + t106 * t71 - t24 * t197 + (-t116 * t173 + t124 * t98) * pkin(3), -t110 * t163 + t106 * t70 - t23 * t197 + (t116 * t103 + t124 * t174) * pkin(3), t167 * t103 + t70 * t221, -t176 * t70 + (t33 - t153 * t36 + (t150 * t88 - t221) * qJD(5)) * t103, -t103 * t165 + t173 * t35 + t271 * t70 + t90 * t71, -t103 * t229 + t173 * t36 - t70 * t224 - t88 * t71, -t173 * t63 + t260 * t71 (-t201 * t58 + t29) * t260 + t56 * t63 - (-t201 * t40 + t26) * t173 + t19 * t71 + t24 * t88 + t57 * t36 + t39 * t190 + ((-qJD(5) * t59 - t23) * t260 - t58 * t63 - (-qJD(5) * t47 - t14) * t173 + t15 * t103 + t240) * t150 -(-t202 * t58 + t196) * t260 - t230 * t63 + t164 * t173 - t20 * t71 + t24 * t90 - t57 * t35 + t153 * t240 + (-t37 + t226) * t103, t16 * t35 - t22 * t36 - t4 * t90 - t5 * t88 + (-t10 * t150 - t153 * t7) * t70 + (-t1 * t153 - t3 * t150 + (-t10 * t153 + t150 * t7) * qJD(5)) * t103, t3 * t22 + t10 * t5 + t1 * t16 + t7 * t4 + t6 * (pkin(5) * t211 + t57) + t31 * ((t150 * t70 + t190) * pkin(5) + t24); 0, 0, 0, 0, 0, -t183, -qJ(2) * t183, 0, 0, 0, 0, 0, 0.2e1 * t122 * qJD(3), t137 + (-t121 - t189) * qJD(3), 0, 0, 0, 0, 0, t182 + t207 + 0.2e1 * t249, -t163 - t252, 0, 0, 0, 0, 0, t169 - t238, -t260 * t271 - t239 - t60 -(t222 - t225) * t98 - t167 + t233, -t31 * t174 + t248 * t153 + (-t260 * t7 + t3) * t150; 0, 0, 0, 0, 0, 0, 0, t121 * t122, -t121 ^ 2 + t122 ^ 2, t137 + (t121 - t189) * qJD(3), 0, 0, -t132 * t122 - t161, t132 * t121 + t170 * t198, t262, t258, t269, t250, 0, t48 * t197 + (-t122 * t98 - t197 * t204) * pkin(3) + t246, -t174 * t242 + t265 + t49 * t197 + (-t197 * t227 - t185) * t154 - t187, t268, -t176 * t98 - t192 + t195, t266, t169 + t238, -t253, t143 * t36 - t51 * t260 + t178 * t88 + (-qJD(5) * t142 * t260 - t15) * t153 + (t260 * t261 + t177) * t150 + t247, -t143 * t35 + (t142 * t202 + t231) * t260 + t178 * t90 + (-t194 * t260 + t177) * t153 + t245, t126 * t35 - t127 * t36 - t217 * t88 - t235 * t90 + t273, t3 * t127 + t1 * t126 + t6 * (t191 - t241) + t235 * t7 + (-t83 + (-t86 + t227) * t151 + t272) * t31 + t217 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, t258, t269, t250, 0, t46 * t197 + t246, t45 * t197 - t14 + t265, t268, -t224 * t90 - t263 * t88 + t195, t266, -t224 * t260 + t238 + t62, -t253, -pkin(4) * t36 - pkin(9) * t229 - t188 * t260 + t264 * t39 - t46 * t88 - t226 + t247, pkin(4) * t35 + t165 * pkin(9) + t232 * t260 + t263 * t39 - t46 * t90 + t245, t135 * t35 - t136 * t36 - t215 * t90 - t216 * t88 + t273, t3 * t136 + t1 * t135 + t6 * t191 + t215 * t7 + (-t46 + t272) * t31 + t216 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 * t88, -t85 + t244, t260 * t88 - t35, t260 * t90 - t36, t63, t20 * t260 - t39 * t90 + t158, t19 * t260 + t39 * t88 - t164, pkin(5) * t35 - t243 * t88, t243 * t10 + (-t31 * t90 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 - t244, t10 * t88 + t7 * t90 + t6;];
tauc_reg  = t2;
