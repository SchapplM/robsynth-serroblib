% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:15
% EndTime: 2019-12-31 18:30:25
% DurationCPUTime: 3.54s
% Computational Cost: add. (4680->310), mult. (12448->425), div. (0->0), fcn. (9555->8), ass. (0->158)
t149 = sin(qJ(5));
t151 = cos(qJ(5));
t147 = cos(pkin(9));
t146 = sin(pkin(8));
t148 = cos(pkin(8));
t150 = sin(qJ(3));
t214 = cos(qJ(3));
t122 = t214 * t146 + t150 * t148;
t112 = t122 * qJD(1);
t145 = sin(pkin(9));
t189 = t145 * t112;
t230 = t147 * qJD(3) - t189;
t161 = t151 * t230;
t94 = t145 * qJD(3) + t147 * t112;
t49 = -t149 * t94 + t161;
t234 = t49 ^ 2;
t175 = qJD(1) * t214;
t134 = t148 * t175;
t185 = t150 * t146;
t176 = qJD(1) * t185;
t110 = -t134 + t176;
t107 = qJD(5) + t110;
t233 = t49 * t107;
t48 = t149 * t230 + t151 * t94;
t232 = t48 ^ 2;
t133 = qJD(3) * t134;
t102 = qJD(3) * t176 - t133;
t188 = t147 * t102;
t231 = t110 * t230 - t188;
t184 = t151 * t147;
t221 = -t149 * t145 + t184;
t222 = t221 * qJD(5);
t201 = -t110 * t221 - t222;
t121 = t151 * t145 + t149 * t147;
t115 = t121 * qJD(5);
t200 = t121 * t110 + t115;
t191 = t145 * t102;
t227 = t110 * t94 - t191;
t117 = t122 * qJD(3);
t103 = qJD(1) * t117;
t159 = t214 * t148 - t185;
t84 = t103 * t159;
t226 = t110 * t117 - t84;
t225 = t112 * t230;
t224 = t121 * t102;
t210 = pkin(6) + qJ(2);
t127 = t210 * t146;
t129 = t210 * t148;
t223 = -t214 * t127 - t150 * t129;
t22 = (qJD(5) * t94 - t191) * t149 - qJD(5) * t161 + t102 * t184;
t220 = -t200 * t48 - t22 * t221;
t219 = t121 * t103 - t201 * t107;
t108 = t110 ^ 2;
t218 = -t145 * t103 - t147 * t108;
t44 = t103 * pkin(3) + t102 * qJ(4) - t112 * qJD(4);
t124 = qJD(1) * t129;
t109 = t150 * t124;
t179 = qJD(1) * qJD(2);
t173 = t150 * t179;
t123 = qJD(1) * t127;
t168 = qJD(2) * t175;
t174 = qJD(3) * t214;
t183 = -t123 * t174 + t148 * t168;
t158 = -t146 * t173 + t183;
t50 = (qJD(4) - t109) * qJD(3) + t158;
t17 = -t145 * t50 + t147 * t44;
t12 = t103 * pkin(4) + pkin(7) * t188 + t17;
t18 = t145 * t44 + t147 * t50;
t13 = pkin(7) * t191 + t18;
t140 = -t148 * pkin(2) - pkin(1);
t125 = t140 * qJD(1) + qJD(2);
t60 = t110 * pkin(3) - t112 * qJ(4) + t125;
t83 = -t150 * t123 + t214 * t124;
t77 = qJD(3) * qJ(4) + t83;
t28 = -t145 * t77 + t147 * t60;
t15 = t110 * pkin(4) - t94 * pkin(7) + t28;
t29 = t145 * t60 + t147 * t77;
t20 = pkin(7) * t230 + t29;
t164 = t149 * t20 - t151 * t15;
t1 = -t164 * qJD(5) + t149 * t12 + t151 * t13;
t217 = t112 ^ 2;
t213 = t147 * pkin(7);
t78 = t112 * pkin(3) + t110 * qJ(4);
t82 = -t214 * t123 - t109;
t34 = -t145 * t82 + t147 * t78;
t21 = t112 * pkin(4) + t110 * t213 + t34;
t196 = t110 * t145;
t35 = t145 * t78 + t147 * t82;
t27 = pkin(7) * t196 + t35;
t209 = pkin(7) + qJ(4);
t126 = t209 * t145;
t128 = t209 * t147;
t89 = -t151 * t126 - t149 * t128;
t216 = qJD(4) * t221 + t89 * qJD(5) - t149 * t21 - t151 * t27;
t91 = -t149 * t126 + t151 * t128;
t215 = -t121 * qJD(4) - t91 * qJD(5) + t149 * t27 - t151 * t21;
t212 = t48 * t49;
t181 = qJD(3) * t150;
t54 = -t123 * t181 + t124 * t174 + t146 * t168 + t148 * t173;
t211 = t54 * t223;
t116 = t146 * t181 - t148 * t174;
t56 = t117 * pkin(3) + t116 * qJ(4) - t122 * qJD(4);
t63 = t159 * qJD(2) + t223 * qJD(3);
t25 = t145 * t56 + t147 * t63;
t79 = -pkin(3) * t159 - t122 * qJ(4) + t140;
t92 = -t150 * t127 + t214 * t129;
t38 = t145 * t79 + t147 * t92;
t207 = t147 * t103 - t145 * t108;
t206 = t112 * t49;
t204 = t48 * t112;
t203 = t94 * t112;
t202 = t94 * t145;
t199 = t102 * t122;
t198 = t110 * t112;
t195 = t116 * t145;
t194 = t116 * t147;
t192 = t122 * t145;
t182 = t146 ^ 2 + t148 ^ 2;
t24 = -t145 * t63 + t147 * t56;
t37 = -t145 * t92 + t147 * t79;
t172 = t182 * qJD(1) ^ 2;
t23 = t48 * qJD(5) - t224;
t170 = -t121 * t23 - t201 * t49;
t169 = t103 * t221 - t200 * t107;
t33 = -pkin(4) * t191 + t54;
t167 = t102 * t223 + t54 * t122;
t165 = -t145 * t28 + t147 * t29;
t6 = t149 * t15 + t151 * t20;
t26 = -pkin(4) * t159 - t122 * t213 + t37;
t30 = -pkin(7) * t192 + t38;
t9 = -t149 * t30 + t151 * t26;
t10 = t149 * t26 + t151 * t30;
t162 = t147 * t230;
t160 = 0.2e1 * t182 * t179;
t75 = -qJD(3) * pkin(3) + qJD(4) - t82;
t157 = -t75 * t116 + t167;
t156 = -t102 * t159 - t122 * t103 + t116 * t110;
t155 = t162 - t202;
t153 = pkin(3) * t102 - qJ(4) * t103 + (-qJD(4) + t75) * t110;
t2 = -t6 * qJD(5) + t151 * t12 - t149 * t13;
t64 = t122 * qJD(2) + qJD(3) * t92;
t143 = t147 ^ 2;
t141 = t145 ^ 2;
t139 = -t147 * pkin(4) - pkin(3);
t72 = t221 * t122;
t71 = t121 * t122;
t65 = pkin(4) * t192 - t223;
t55 = -pkin(4) * t196 + t83;
t53 = (-qJD(3) * t124 - t146 * t179) * t150 + t183;
t40 = -pkin(4) * t230 + t75;
t39 = -pkin(4) * t195 + t64;
t32 = -t121 * t116 + t222 * t122;
t31 = t115 * t122 + t221 * t116;
t16 = pkin(7) * t195 + t25;
t14 = t117 * pkin(4) + pkin(7) * t194 + t24;
t4 = -t10 * qJD(5) + t151 * t14 - t149 * t16;
t3 = t9 * qJD(5) + t149 * t14 + t151 * t16;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, qJ(2) * t160, -t112 * t116 - t199, -t112 * t117 + t156, -t116 * qJD(3), t226, -t117 * qJD(3), 0, -t64 * qJD(3) + t140 * t103 + t125 * t117, -t63 * qJD(3) - t140 * t102 - t125 * t116, -t92 * t103 - t63 * t110 + t64 * t112 + t82 * t116 - t83 * t117 + t159 * t53 + t167, t53 * t92 + t83 * t63 - t82 * t64 - t211, -t143 * t199 - t94 * t194, -t116 * t155 + 0.2e1 * t188 * t192, t94 * t117 - t147 * t156, -t141 * t199 + t195 * t230, t117 * t230 + t145 * t156, t226, t37 * t103 + t24 * t110 + t28 * t117 + t145 * t157 - t159 * t17 - t230 * t64, -t38 * t103 - t25 * t110 - t29 * t117 + t147 * t157 + t159 * t18 + t64 * t94, -t25 * t189 - t24 * t94 + (t25 * qJD(3) + t37 * t102 + t28 * t116 - t17 * t122) * t147 + (t38 * t102 + t29 * t116 - t18 * t122) * t145, t17 * t37 + t18 * t38 + t28 * t24 + t29 * t25 + t75 * t64 - t211, -t22 * t72 - t48 * t31, t22 * t71 - t72 * t23 - t31 * t49 - t48 * t32, t72 * t103 - t31 * t107 + t48 * t117 + t159 * t22, t23 * t71 - t32 * t49, -t71 * t103 - t32 * t107 + t117 * t49 + t159 * t23, t107 * t117 - t84, t9 * t103 + t4 * t107 - t117 * t164 - t159 * t2 + t65 * t23 + t40 * t32 + t33 * t71 - t39 * t49, t1 * t159 - t10 * t103 - t3 * t107 - t6 * t117 - t65 * t22 - t40 * t31 + t33 * t72 + t39 * t48, -t1 * t71 - t10 * t23 - t164 * t31 - t2 * t72 + t9 * t22 + t3 * t49 - t6 * t32 - t4 * t48, t1 * t10 - t164 * t4 + t2 * t9 + t6 * t3 + t33 * t65 + t40 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -qJ(2) * t172, 0, 0, 0, 0, 0, 0, 0.2e1 * t112 * qJD(3), t133 + (-t110 - t176) * qJD(3), -t108 - t217, t83 * t110 + t82 * t112, 0, 0, 0, 0, 0, 0, t207 + t225, -t203 + t218, (t162 + t202) * t110 + (t141 + t143) * t102, t165 * t110 - t75 * t112 + t18 * t145 + t17 * t147, 0, 0, 0, 0, 0, 0, t169 + t206, -t204 - t219, t170 - t220, t1 * t121 - t40 * t112 + t164 * t200 + t2 * t221 - t201 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, -t108 + t217, t133 + (t110 - t176) * qJD(3), -t198, 0, 0, t83 * qJD(3) - t125 * t112 - t54, t125 * t110 + (t82 + t109) * qJD(3) - t158, 0, 0, t227 * t147, t155 * t110 + (t141 - t143) * t102, -t203 - t218, -t231 * t145, t207 - t225, -t198, -t34 * t110 - t28 * t112 + t145 * t153 - t54 * t147 + t230 * t83, t35 * t110 + t29 * t112 + t54 * t145 + t147 * t153 - t83 * t94, t35 * t189 + t34 * t94 + (-qJD(4) * t189 - t28 * t110 + t18 + (t147 * qJD(4) - t35) * qJD(3)) * t147 + (qJD(4) * t94 - t29 * t110 - t17) * t145, -t54 * pkin(3) - t28 * t34 - t29 * t35 - t75 * t83 + t165 * qJD(4) + (-t17 * t145 + t18 * t147) * qJ(4), -t22 * t121 - t201 * t48, t170 + t220, -t204 + t219, -t200 * t49 - t221 * t23, t169 - t206, -t107 * t112, t89 * t103 + t215 * t107 + t112 * t164 + t139 * t23 + t200 * t40 - t221 * t33 + t49 * t55, -t91 * t103 - t216 * t107 + t6 * t112 + t33 * t121 - t139 * t22 - t201 * t40 - t55 * t48, t1 * t221 - t2 * t121 - t164 * t201 - t200 * t6 - t215 * t48 + t216 * t49 + t89 * t22 - t91 * t23, t1 * t91 + t33 * t139 - t164 * t215 + t2 * t89 + t216 * t6 - t40 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t231, -t230 ^ 2 - t94 ^ 2, -t230 * t29 + t28 * t94 + t54, 0, 0, 0, 0, 0, 0, t48 * t107 + t23, -t22 + t233, -t232 - t234, -t164 * t48 - t6 * t49 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, t232 - t234, -t22 - t233, t212, t224 + (-qJD(5) + t107) * t48, t103, t6 * t107 - t40 * t48 + t2, -t107 * t164 - t40 * t49 - t1, 0, 0;];
tauc_reg = t5;
