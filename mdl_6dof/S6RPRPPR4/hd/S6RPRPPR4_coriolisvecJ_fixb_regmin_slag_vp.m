% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:23
% EndTime: 2019-03-09 02:48:29
% DurationCPUTime: 2.55s
% Computational Cost: add. (3578->338), mult. (9609->454), div. (0->0), fcn. (7370->8), ass. (0->164)
t158 = cos(pkin(9));
t231 = cos(qJ(3));
t200 = qJD(1) * t231;
t146 = t158 * t200;
t156 = sin(pkin(9));
t160 = sin(qJ(3));
t216 = t160 * t156;
t201 = qJD(1) * t216;
t120 = -t146 + t201;
t206 = qJD(6) - t120;
t250 = qJD(6) - t206;
t155 = sin(pkin(10));
t157 = cos(pkin(10));
t215 = t160 * t158;
t134 = t231 * t156 + t215;
t238 = t134 * qJD(1);
t90 = -t157 * qJD(3) + t155 * t238;
t249 = t120 * t90;
t159 = sin(qJ(6));
t161 = cos(qJ(6));
t180 = t155 * qJD(3) + t157 * t238;
t52 = t159 * t180 - t161 * t90;
t248 = t206 * t52;
t54 = t159 * t90 + t161 * t180;
t247 = t206 * t54;
t133 = t161 * t155 - t159 * t157;
t207 = qJD(6) * t161;
t208 = qJD(6) * t159;
t223 = -t133 * t120 + t155 * t207 - t157 * t208;
t229 = pkin(7) + qJ(2);
t140 = t229 * t156;
t142 = t229 * t158;
t89 = -t160 * t140 + t231 * t142;
t64 = qJD(2) * t134 + qJD(3) * t89;
t246 = t157 * t134 * qJD(5) - t64;
t70 = t133 * t134;
t135 = qJD(1) * t140;
t136 = qJD(1) * t142;
t81 = -t160 * t135 + t231 * t136;
t245 = t81 * qJD(3);
t131 = t159 * t155 + t161 * t157;
t244 = t238 * qJD(3);
t243 = t180 ^ 2;
t126 = t134 * qJD(3);
t107 = qJD(1) * t126;
t174 = t231 * t158 - t216;
t167 = t174 * qJD(2);
t175 = t231 * t135 + t160 * t136;
t55 = qJD(1) * t167 + (qJD(4) - t175) * qJD(3);
t49 = t155 * t55;
t145 = qJD(3) * t146;
t106 = -qJD(3) * t201 + t145;
t51 = t107 * pkin(3) - t106 * qJ(4) - qJD(4) * t238;
t17 = t157 * t51 - t49;
t11 = -t107 * pkin(4) - t17;
t149 = -t158 * pkin(2) - pkin(1);
t137 = t149 * qJD(1) + qJD(2);
t62 = t120 * pkin(3) - qJ(4) * t238 + t137;
t76 = qJD(3) * qJ(4) + t81;
t36 = t155 * t62 + t157 * t76;
t26 = t120 * qJ(5) + t36;
t242 = -t26 * t120 + t11;
t239 = -t231 * t140 - t160 * t142;
t24 = qJD(6) * t54 - t133 * t106;
t224 = t206 * t131;
t237 = t133 * t107 + t206 * t224;
t115 = t120 ^ 2;
t236 = t157 * t107 - t155 * t115 - t238 * t90;
t235 = t155 * t107 + t157 * t115 + t180 * t238;
t199 = t155 * qJ(5) + pkin(3);
t138 = -t157 * pkin(4) - t199;
t222 = qJ(4) * t107;
t178 = qJD(3) * pkin(3) - qJD(4) - t175;
t166 = qJ(5) * t180 + t178;
t34 = t90 * pkin(4) - t166;
t234 = -t106 * t138 + (qJD(4) - t34) * t120 + t222;
t221 = qJ(5) * t157;
t232 = pkin(4) + pkin(5);
t233 = t232 * t155 - t221;
t230 = pkin(8) * t155;
t228 = -pkin(8) + qJ(4);
t18 = t155 * t51 + t157 * t55;
t125 = t174 * qJD(3);
t58 = t126 * pkin(3) - t125 * qJ(4) - t134 * qJD(4);
t63 = t239 * qJD(3) + t167;
t29 = t155 * t58 + t157 * t63;
t77 = pkin(3) * t238 + t120 * qJ(4);
t43 = t155 * t77 - t157 * t175;
t78 = -pkin(3) * t174 - t134 * qJ(4) + t149;
t46 = t155 * t78 + t157 * t89;
t227 = t238 * t52;
t225 = t54 * t238;
t220 = qJD(4) * t180;
t99 = t155 * t106;
t100 = t157 * t106;
t212 = t180 * qJD(5);
t210 = t156 ^ 2 + t158 ^ 2;
t209 = qJD(4) * t157;
t205 = qJD(1) * qJD(2);
t31 = qJ(5) * t238 + t43;
t39 = -qJ(5) * t174 + t46;
t4 = t49 + (-pkin(8) * t106 - t51) * t157 - t232 * t107;
t6 = t107 * qJ(5) + t120 * qJD(5) + t18;
t5 = pkin(8) * t99 + t6;
t204 = -t159 * t5 + t161 * t4;
t198 = t120 * t180 + t99;
t60 = t155 * t63;
t28 = t157 * t58 - t60;
t35 = -t155 * t76 + t157 * t62;
t74 = t155 * t175;
t42 = t157 * t77 + t74;
t83 = t155 * t89;
t45 = t157 * t78 - t83;
t197 = t100 - t249;
t196 = t210 * qJD(1) ^ 2;
t195 = t155 * qJD(5) - t233 * t120 + t81;
t194 = t206 ^ 2;
t12 = t126 * qJ(5) - qJD(5) * t174 + t29;
t57 = t156 * qJD(2) * t200 + t205 * t215 + t245;
t193 = t131 * t107 - t223 * t206;
t191 = qJD(5) - t35;
t190 = t159 * t4 + t161 * t5;
t189 = -t90 ^ 2 - t243;
t13 = t90 * pkin(8) + t26;
t9 = -pkin(8) * t180 - t232 * t120 + t191;
t2 = t161 * t13 + t159 * t9;
t188 = t159 * t13 - t161 * t9;
t187 = pkin(4) * t155 - t221;
t186 = -t155 * t35 + t157 * t36;
t20 = t83 + (-pkin(8) * t134 - t78) * t157 + t232 * t174;
t30 = t134 * t230 + t39;
t185 = -t159 * t30 + t161 * t20;
t184 = t159 * t20 + t161 * t30;
t179 = 0.2e1 * t210 * t205;
t176 = -t131 * t106 + t180 * t208 - t90 * t207;
t139 = t228 * t155;
t173 = qJD(6) * t139 + t120 * t230 + t209 - t31;
t141 = t228 * t157;
t172 = qJD(4) * t155 - qJD(6) * t141 + t74 - (pkin(8) * t120 - t77) * t157 + t232 * t238;
t71 = t131 * t134;
t169 = -pkin(4) * t99 + t212 - t57;
t15 = -qJ(5) * t100 - t169;
t48 = t134 * t187 - t239;
t171 = t106 * t48 + t125 * t34 + t134 * t15;
t170 = -t106 * t239 - t125 * t178 + t134 * t57;
t165 = -pkin(3) * t106 - t222 + (-qJD(4) - t178) * t120;
t164 = (t155 * t180 - t157 * t90) * t120 + (-t155 ^ 2 - t157 ^ 2) * t106;
t127 = t232 * t157 + t199;
t85 = t90 * t209;
t44 = -t120 * t187 + t81;
t41 = -t134 * t233 + t239;
t40 = pkin(4) * t174 - t45;
t38 = qJD(6) * t71 - t133 * t125;
t37 = qJD(6) * t70 + t125 * t131;
t33 = -pkin(4) * t238 - t42;
t27 = t125 * t187 - t246;
t25 = -t120 * pkin(4) + t191;
t21 = -t232 * t90 + t166;
t19 = -t126 * pkin(4) - t28;
t16 = -t125 * t233 + t246;
t10 = (-pkin(5) * t155 + t221) * t106 + t169;
t8 = t125 * t230 + t12;
t7 = t60 + (-pkin(8) * t125 - t58) * t157 - t232 * t126;
t1 = [0, 0, 0, 0, 0, t179, qJ(2) * t179, t106 * t134 + t125 * t238, t106 * t174 - t134 * t107 - t125 * t120 - t126 * t238, t125 * qJD(3), -t126 * qJD(3), 0, -t64 * qJD(3) + t149 * t107 + t137 * t126, -t63 * qJD(3) + t149 * t106 + t137 * t125, t45 * t107 + t28 * t120 + t35 * t126 + t155 * t170 - t17 * t174 + t64 * t90, -t46 * t107 - t29 * t120 - t36 * t126 + t157 * t170 + t174 * t18 + t180 * t64, -t28 * t180 - t29 * t90 + (-t106 * t45 - t125 * t35 - t134 * t17) * t157 + (-t106 * t46 - t125 * t36 - t134 * t18) * t155, t17 * t45 - t178 * t64 + t18 * t46 - t239 * t57 + t35 * t28 + t36 * t29, -t40 * t107 + t11 * t174 - t19 * t120 - t25 * t126 + t155 * t171 + t27 * t90, -t12 * t90 + t19 * t180 + (t106 * t40 + t11 * t134 + t125 * t25) * t157 + (-t106 * t39 - t125 * t26 - t134 * t6) * t155, t39 * t107 + t12 * t120 + t26 * t126 - t157 * t171 - t174 * t6 - t180 * t27, t11 * t40 + t26 * t12 + t15 * t48 + t25 * t19 + t34 * t27 + t6 * t39, -t176 * t71 + t54 * t37, -t176 * t70 - t71 * t24 - t37 * t52 - t54 * t38, -t71 * t107 - t54 * t126 - t174 * t176 + t206 * t37, -t70 * t107 + t52 * t126 - t174 * t24 - t206 * t38, -t107 * t174 - t126 * t206 (-t159 * t8 + t161 * t7) * t206 - t185 * t107 + t204 * t174 + t188 * t126 + t16 * t52 + t41 * t24 - t10 * t70 + t21 * t38 + (-t174 * t2 - t184 * t206) * qJD(6) -(t159 * t7 + t161 * t8) * t206 + t184 * t107 - t190 * t174 + t2 * t126 + t16 * t54 - t41 * t176 + t10 * t71 + t21 * t37 + (t174 * t188 - t185 * t206) * qJD(6); 0, 0, 0, 0, 0, -t196, -qJ(2) * t196, 0, 0, 0, 0, 0, 0.2e1 * t244, t145 + (-t120 - t201) * qJD(3), t236, -t235, t164, t120 * t186 + t18 * t155 + t17 * t157 + t178 * t238, t236, t164, t235, -t11 * t157 - t34 * t238 + t6 * t155 + (t155 * t25 + t157 * t26) * t120, 0, 0, 0, 0, 0, t193 + t227, t225 + t237; 0, 0, 0, 0, 0, 0, 0, t238 * t120, t238 ^ 2 - t115, t145 + (t120 - t201) * qJD(3), 0, 0, -t137 * t238 + t245 - t57, t137 * t120 - t174 * t205, -t42 * t120 + t155 * t165 - t57 * t157 - t238 * t35 - t81 * t90, t43 * t120 + t57 * t155 + t157 * t165 - t180 * t81 + t238 * t36, t42 * t180 + t43 * t90 - t85 + (-t120 * t35 + t18) * t157 + (-t120 * t36 - t17 + t220) * t155, -t57 * pkin(3) - t35 * t42 - t36 * t43 + t178 * t81 + t186 * qJD(4) + (-t17 * t155 + t18 * t157) * qJ(4), t33 * t120 + t25 * t238 - t15 * t157 - t44 * t90 + (-qJD(5) * t90 - t234) * t155, t31 * t90 - t33 * t180 - t85 + (t120 * t25 + t6) * t157 + (t220 + t242) * t155, -t31 * t120 - t26 * t238 + t44 * t180 + (-t15 + t212) * t155 + t234 * t157, t15 * t138 - t25 * t33 - t26 * t31 - t34 * t44 + (qJ(4) * t6 + qJD(4) * t26) * t157 + (qJ(4) * t11 + qJD(4) * t25 - qJD(5) * t34) * t155, -t133 * t176 - t224 * t54, t131 * t176 - t133 * t24 - t223 * t54 + t224 * t52, t225 - t237, t193 - t227, t206 * t238 -(t161 * t139 - t159 * t141) * t107 + t127 * t24 + t10 * t131 - t188 * t238 + t195 * t52 + t223 * t21 - (t159 * t173 - t161 * t172) * t206 (t159 * t139 + t161 * t141) * t107 - t127 * t176 + t10 * t133 - t2 * t238 + t195 * t54 - t224 * t21 - (t159 * t172 + t161 * t173) * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t197, t189, t180 * t35 + t36 * t90 + t57, t198, t189, -t197, -t180 * t25 + t26 * t90 + t15, 0, 0, 0, 0, 0, -t24 - t247, t176 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180 * t90 - t244, t100 + t249, -t115 - t243, t180 * t34 + t242, 0, 0, 0, 0, 0, -t161 * t107 - t159 * t194 - t180 * t52, t159 * t107 - t161 * t194 - t180 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t52 ^ 2 + t54 ^ 2, -t176 + t248, -t24 + t247, -t107, -t250 * t2 - t21 * t54 + t204, t250 * t188 + t21 * t52 - t190;];
tauc_reg  = t1;
