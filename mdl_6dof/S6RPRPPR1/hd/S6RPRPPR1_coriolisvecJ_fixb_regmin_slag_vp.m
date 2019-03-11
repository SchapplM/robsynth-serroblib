% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:34
% EndTime: 2019-03-09 02:39:41
% DurationCPUTime: 2.16s
% Computational Cost: add. (3275->280), mult. (8061->406), div. (0->0), fcn. (5855->10), ass. (0->156)
t147 = cos(qJ(3));
t200 = cos(pkin(10));
t176 = t200 * t147;
t127 = qJD(1) * t176;
t140 = sin(pkin(10));
t145 = sin(qJ(3));
t184 = t145 * qJD(1);
t107 = t140 * t184 - t127;
t104 = qJD(6) + t107;
t216 = qJD(6) - t104;
t146 = cos(qJ(6));
t121 = t140 * t147 + t145 * t200;
t110 = t121 * qJD(1);
t139 = sin(pkin(11));
t142 = cos(pkin(11));
t94 = t142 * qJD(3) - t139 * t110;
t215 = t146 * t94;
t144 = sin(qJ(6));
t122 = t146 * t139 + t144 * t142;
t201 = t104 * t122;
t93 = t139 * qJD(3) + t142 * t110;
t162 = -t144 * t94 - t146 * t93;
t214 = t104 * t162;
t132 = sin(pkin(9)) * pkin(1) + pkin(7);
t188 = qJ(4) + t132;
t174 = t188 * qJD(1);
t96 = t147 * qJD(2) - t174 * t145;
t109 = t121 * qJD(3);
t101 = qJD(1) * t109;
t191 = t146 * t142;
t192 = t144 * t139;
t120 = -t191 + t192;
t202 = t104 * t120;
t213 = -t122 * t101 + t104 * t202;
t97 = t145 * qJD(2) + t147 * t174;
t182 = qJD(1) * qJD(3);
t177 = t145 * t182;
t102 = qJD(3) * t127 - t140 * t177;
t21 = -qJD(6) * t162 + t102 * t122;
t105 = t107 ^ 2;
t212 = t142 * pkin(8);
t181 = qJD(1) * qJD(4);
t150 = -t97 * qJD(3) - t145 * t181;
t80 = t96 * qJD(3) + t147 * t181;
t35 = t140 * t80 - t200 * t150;
t117 = t188 * t145;
t118 = t188 * t147;
t78 = t200 * t117 + t140 * t118;
t211 = t35 * t78;
t129 = t140 * pkin(3) + qJ(5);
t210 = pkin(8) + t129;
t155 = -t140 * t145 + t176;
t194 = t139 * t102;
t161 = -qJD(6) * t93 - t194;
t185 = qJD(6) * t146;
t207 = t102 * t191 + t185 * t94;
t20 = t144 * t161 + t207;
t209 = -t109 * t162 - t155 * t20;
t112 = t155 * qJD(3);
t186 = qJD(6) * t121;
t196 = t121 * t142;
t29 = t112 * t122 + t185 * t196 - t186 * t192;
t73 = t122 * t121;
t208 = -t73 * t101 - t29 * t104;
t36 = t140 * t150 + t200 * t80;
t33 = qJD(3) * qJD(5) + t36;
t128 = pkin(3) * t177;
t41 = t101 * pkin(4) - t102 * qJ(5) - t110 * qJD(5) + t128;
t9 = t139 * t41 + t142 * t33;
t178 = t200 * t97;
t206 = qJD(3) * pkin(3);
t89 = t96 + t206;
t49 = t140 * t89 + t178;
t44 = qJD(3) * qJ(5) + t49;
t134 = -cos(pkin(9)) * pkin(1) - pkin(2);
t159 = -t147 * pkin(3) + t134;
t154 = t159 * qJD(1);
t106 = qJD(4) + t154;
t63 = t107 * pkin(4) - t110 * qJ(5) + t106;
t16 = t139 * t63 + t142 * t44;
t86 = t140 * t97;
t56 = t200 * t96 - t86;
t70 = pkin(3) * t184 + t110 * pkin(4) + t107 * qJ(5);
t23 = t139 * t70 + t142 * t56;
t180 = t145 * t206;
t54 = t109 * pkin(4) - t112 * qJ(5) - t121 * qJD(5) + t180;
t175 = qJD(3) * t188;
t98 = t147 * qJD(4) - t145 * t175;
t99 = -t145 * qJD(4) - t147 * t175;
t62 = t140 * t99 + t200 * t98;
t19 = t139 * t54 + t142 * t62;
t72 = -pkin(4) * t155 - t121 * qJ(5) + t159;
t79 = -t140 * t117 + t118 * t200;
t31 = t139 * t72 + t142 * t79;
t51 = t144 * t93 - t215;
t205 = t110 * t51;
t204 = t35 * t155;
t203 = t162 * t110;
t199 = t107 * t139;
t198 = t112 * t139;
t197 = t121 * t139;
t193 = t142 * t102;
t148 = qJD(3) ^ 2;
t190 = t148 * t145;
t189 = t148 * t147;
t187 = t145 ^ 2 - t147 ^ 2;
t125 = qJD(1) * t134;
t8 = -t139 * t33 + t142 * t41;
t4 = t101 * pkin(5) - pkin(8) * t193 + t8;
t5 = -pkin(8) * t194 + t9;
t179 = -t144 * t5 + t146 * t4;
t15 = -t139 * t44 + t142 * t63;
t18 = -t139 * t62 + t142 * t54;
t22 = -t139 * t56 + t142 * t70;
t30 = -t139 * t79 + t142 * t72;
t55 = t140 * t96 + t178;
t61 = t140 * t98 - t200 * t99;
t173 = -t120 * t101 - t104 * t201;
t133 = -pkin(3) * t200 - pkin(4);
t172 = -t8 * t139 + t9 * t142;
t171 = t144 * t4 + t146 * t5;
t48 = t200 * t89 - t86;
t10 = pkin(8) * t94 + t16;
t6 = t107 * pkin(5) - t93 * pkin(8) + t15;
t2 = t146 * t10 + t144 * t6;
t170 = t144 * t10 - t146 * t6;
t28 = -t112 * t120 - t122 * t186;
t74 = t120 * t121;
t169 = t101 * t74 - t104 * t28;
t168 = t78 * t102 + t35 * t121;
t167 = -t109 * t51 + t155 * t21;
t166 = -t139 * t15 + t142 * t16;
t165 = t139 * t93 + t142 * t94;
t17 = -pkin(5) * t155 - pkin(8) * t196 + t30;
t25 = -pkin(8) * t197 + t31;
t164 = -t144 * t25 + t146 * t17;
t163 = t144 * t17 + t146 * t25;
t158 = 0.2e1 * qJD(3) * t125;
t116 = t210 * t142;
t157 = t110 * pkin(5) + qJD(5) * t139 + qJD(6) * t116 + t107 * t212 + t22;
t115 = t210 * t139;
t156 = pkin(8) * t199 - qJD(5) * t142 + qJD(6) * t115 + t23;
t43 = -qJD(3) * pkin(4) + qJD(5) - t48;
t153 = t112 * t43 + t168;
t152 = -t121 * t101 - t102 * t155 - t112 * t107;
t151 = -t101 * t129 + t102 * t133 + (-qJD(5) + t43) * t107;
t149 = qJD(1) ^ 2;
t123 = -t142 * pkin(5) + t133;
t60 = pkin(5) * t197 + t78;
t38 = pkin(5) * t198 + t61;
t34 = -pkin(5) * t199 + t55;
t27 = -pkin(5) * t94 + t43;
t26 = pkin(5) * t194 + t35;
t12 = -pkin(8) * t198 + t19;
t7 = t109 * pkin(5) - t112 * t212 + t18;
t1 = [0, 0, 0, 0, 0.2e1 * t147 * t177, -0.2e1 * t187 * t182, t189, -t190, 0, -t132 * t189 + t145 * t158, t132 * t190 + t147 * t158, -t79 * t101 - t62 * t107 - t49 * t109 + t61 * t110 - t48 * t112 + t155 * t36 + t168, t211 + t36 * t79 - t48 * t61 + t49 * t62 + (t106 + t154) * t180, t30 * t101 + t18 * t107 + t15 * t109 + t139 * t153 - t155 * t8 - t61 * t94, -t31 * t101 - t19 * t107 - t16 * t109 + t142 * t153 + t155 * t9 + t61 * t93, -t18 * t93 + t19 * t94 + (-t102 * t30 - t112 * t15 - t121 * t8) * t142 + (-t102 * t31 - t112 * t16 - t121 * t9) * t139, t15 * t18 + t16 * t19 + t30 * t8 + t31 * t9 + t43 * t61 + t211, -t162 * t28 - t20 * t74, t162 * t29 - t20 * t73 + t21 * t74 - t28 * t51, -t169 + t209, t167 + t208, -t101 * t155 + t104 * t109 (-t144 * t12 + t146 * t7) * t104 + t164 * t101 - t179 * t155 - t170 * t109 + t38 * t51 + t60 * t21 + t26 * t73 + t27 * t29 + (-t104 * t163 + t155 * t2) * qJD(6) -(t146 * t12 + t144 * t7) * t104 - t163 * t101 + t171 * t155 - t2 * t109 - t38 * t162 + t60 * t20 - t26 * t74 + t27 * t28 + (-t104 * t164 - t155 * t170) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t189, t109 * t110 + t152, -t48 * t109 + t49 * t112 + t36 * t121 - t204, -t109 * t94 + t139 * t152, t109 * t93 + t142 * t152, t165 * t112, t43 * t109 + t112 * t166 + t121 * t172 - t204, 0, 0, 0, 0, 0, -t167 + t208, t169 + t209; 0, 0, 0, 0, -t145 * t149 * t147, t187 * t149, 0, 0, 0, -t125 * t184, -t125 * t147 * qJD(1) (t49 - t55) * t110 + (-t48 + t56) * t107 + (-t101 * t140 - t102 * t200) * pkin(3), t48 * t55 - t49 * t56 + (-t106 * t184 + t140 * t36 - t200 * t35) * pkin(3), -t22 * t107 - t15 * t110 + t139 * t151 - t35 * t142 + t55 * t94, t23 * t107 + t16 * t110 + t35 * t139 + t142 * t151 - t55 * t93, t22 * t93 - t23 * t94 + (qJD(5) * t94 - t107 * t15 + t9) * t142 + (qJD(5) * t93 - t107 * t16 - t8) * t139, qJD(5) * t166 + t129 * t172 + t35 * t133 - t15 * t22 - t16 * t23 - t43 * t55, t20 * t122 + t162 * t202, -t20 * t120 - t122 * t21 + t162 * t201 + t202 * t51, t203 - t213, t173 + t205, -t104 * t110 (-t146 * t115 - t144 * t116) * t101 + t123 * t21 + t26 * t120 + t170 * t110 - t34 * t51 + t201 * t27 + (t144 * t156 - t146 * t157) * t104 -(-t144 * t115 + t146 * t116) * t101 + t123 * t20 + t26 * t122 + t2 * t110 + t34 * t162 - t202 * t27 + (t144 * t157 + t146 * t156) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 ^ 2 - t105, t49 * t107 + t48 * t110 + t128, t142 * t101 - t139 * t105 + t110 * t94, -t139 * t101 - t142 * t105 - t110 * t93, t165 * t107 + (-t139 ^ 2 - t142 ^ 2) * t102, t107 * t166 - t43 * t110 + t9 * t139 + t8 * t142, 0, 0, 0, 0, 0, t173 - t205, t203 + t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t93 + t194, t94 * t107 + t193, -t93 ^ 2 - t94 ^ 2, t15 * t93 - t16 * t94 + t35, 0, 0, 0, 0, 0, t21 - t214, t104 * t215 + (-t104 * t93 + t161) * t144 + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162 * t51, t162 ^ 2 - t51 ^ 2, t51 * t104 + t20, -t21 - t214, t101, t27 * t162 - t2 * t216 + t179, t170 * t216 + t27 * t51 - t171;];
tauc_reg  = t1;
