% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:23:00
% EndTime: 2020-01-03 11:23:06
% DurationCPUTime: 1.26s
% Computational Cost: add. (1040->173), mult. (2810->295), div. (0->0), fcn. (3066->8), ass. (0->169)
t119 = cos(pkin(9));
t117 = sin(pkin(8));
t122 = sin(qJ(5));
t197 = t117 * t122;
t151 = t119 * t197;
t120 = cos(pkin(8));
t123 = cos(qJ(5));
t188 = t123 * t120;
t134 = t151 + t188;
t116 = sin(pkin(9));
t121 = cos(pkin(7));
t190 = t121 * t119;
t118 = sin(pkin(7));
t193 = t120 * t118;
t81 = t116 * t193 + t190;
t217 = t134 * t81;
t79 = t81 ^ 2;
t216 = -t116 / 0.2e1;
t215 = -t119 / 0.2e1;
t196 = t117 * t123;
t191 = t121 * t116;
t83 = t119 * t193 - t191;
t55 = -t118 * t196 + t83 * t122;
t57 = t118 * t197 + t83 * t123;
t214 = t55 * t57;
t192 = t120 * t121;
t92 = -t121 * pkin(2) - t118 * qJ(3) - pkin(1);
t67 = qJ(2) * t192 + t117 * t92;
t58 = -t121 * qJ(4) + t67;
t76 = (pkin(3) * t117 - qJ(4) * t120 + qJ(2)) * t118;
t33 = t116 * t76 + t119 * t58;
t199 = t117 * t118;
t32 = -t116 * t58 + t119 * t76;
t29 = -pkin(4) * t199 - t32;
t198 = t117 * t121;
t66 = -qJ(2) * t198 + t120 * t92;
t139 = t121 * pkin(3) - t66;
t28 = t81 * pkin(4) - t83 * pkin(6) + t139;
t30 = pkin(6) * t199 + t33;
t6 = t122 * t30 - t123 * t28;
t1 = -t29 * t55 + t6 * t81;
t213 = t1 * qJD(1);
t85 = t116 * t118 + t120 * t190;
t212 = t122 * t85;
t211 = t123 * t85;
t7 = -t122 * t28 - t123 * t30;
t2 = t29 * t57 + t7 * t81;
t210 = t2 * qJD(1);
t195 = t119 * t118;
t80 = t120 * t191 - t195;
t3 = t139 * t198 - t32 * t80 + t33 * t85;
t209 = t3 * qJD(1);
t153 = t117 * t195;
t154 = t116 * t199;
t4 = -t139 * t193 + t33 * t153 - t32 * t154;
t208 = t4 * qJD(1);
t5 = -t32 * t83 - t33 * t81;
t207 = t5 * qJD(1);
t189 = t122 * t120;
t135 = -t119 * t196 + t189;
t129 = t135 * t81;
t200 = t116 * t117;
t156 = t57 * t200;
t125 = -t129 / 0.2e1 - t156 / 0.2e1;
t145 = t196 / 0.2e1;
t131 = -t212 / 0.2e1 + t121 * t145;
t8 = t125 + t131;
t206 = t8 * qJD(1);
t157 = t55 * t200;
t124 = -t217 / 0.2e1 + t157 / 0.2e1;
t152 = t121 * t197;
t130 = -t211 / 0.2e1 - t152 / 0.2e1;
t9 = t124 + t130;
t205 = t9 * qJD(1);
t204 = qJD(1) * t81;
t203 = qJD(5) * t81;
t112 = t117 ^ 2;
t113 = t118 ^ 2;
t202 = t112 * t113;
t201 = t112 * t118;
t12 = (t121 * t196 - t212) * t81 + t80 * t55;
t194 = t12 * qJD(1);
t13 = -(t152 + t211) * t81 + t80 * t57;
t187 = t13 * qJD(1);
t14 = t55 ^ 2 - t57 ^ 2;
t186 = t14 * qJD(1);
t15 = -t122 * t79 - t83 * t55;
t185 = t15 * qJD(1);
t127 = (t151 / 0.2e1 + t188 / 0.2e1) * t118;
t150 = t81 * t216;
t133 = t123 * t150 + t57 * t215;
t16 = t127 - t133;
t184 = t16 * qJD(1);
t126 = (t119 * t145 - t189 / 0.2e1) * t118;
t132 = t122 * t150 + t55 * t215;
t17 = t126 + t132;
t183 = t17 * qJD(1);
t20 = (-t157 + t217) * t118;
t182 = t20 * qJD(1);
t21 = (-t129 - t156) * t118;
t181 = t21 * qJD(1);
t114 = t120 ^ 2;
t104 = -t114 * t118 / 0.2e1;
t128 = t104 + (-t119 ^ 2 / 0.2e1 - t116 ^ 2 / 0.2e1) * t201;
t136 = t80 * t119 / 0.2e1 + t85 * t216;
t23 = t128 + t136;
t180 = t23 * qJD(1);
t24 = t80 * t83 - t85 * t81;
t179 = t24 * qJD(1);
t25 = t123 * t79 + t83 * t57;
t178 = t25 * qJD(1);
t27 = t113 * qJ(2) + (-t117 * t66 + t120 * t67) * t121;
t177 = t27 * qJD(1);
t31 = (t117 * t67 + t120 * t66) * t118;
t176 = t31 * qJD(1);
t137 = t81 * t215 + t116 * t83 / 0.2e1;
t35 = (-t121 / 0.2e1 + t137) * t117;
t175 = t35 * qJD(1);
t36 = -t81 * t153 + t83 * t154;
t174 = t36 * qJD(1);
t138 = t83 * t215 + t150;
t38 = -t193 / 0.2e1 + t138;
t173 = t38 * qJD(1);
t39 = (-t118 * t80 + t121 * t81) * t117;
t172 = t39 * qJD(1);
t40 = -t83 * t198 + t85 * t199;
t171 = t40 * qJD(1);
t42 = t122 * t81;
t170 = t42 * qJD(1);
t46 = t83 ^ 2 + t79;
t169 = t46 * qJD(1);
t52 = t116 * t202 + t81 * t193;
t168 = t52 * qJD(1);
t53 = t119 * t202 + t83 * t193;
t167 = t53 * qJD(1);
t77 = (0.1e1 / 0.2e1 + t112 / 0.2e1 + t114 / 0.2e1) * t118;
t166 = t77 * qJD(1);
t87 = (t112 + t114) * t113;
t165 = t87 * qJD(1);
t101 = t121 ^ 2 + t113;
t88 = t101 * t117;
t164 = t88 * qJD(1);
t89 = t101 * t120;
t163 = t89 * qJD(1);
t99 = t101 * qJ(2);
t162 = t99 * qJD(1);
t161 = qJD(1) * t118;
t160 = qJD(5) * t122;
t159 = qJD(5) * t123;
t158 = t101 * qJD(1);
t155 = qJD(1) * t214;
t149 = t117 * t161;
t148 = t121 * t161;
t147 = qJD(3) * t118 * t121;
t146 = qJD(4) * t199;
t144 = -qJD(5) - t204;
t143 = t81 * t149;
t142 = t83 * t149;
t141 = t117 * t148;
t140 = t120 * t148;
t78 = -t201 / 0.2e1 + t104 + t118 / 0.2e1;
t37 = t193 / 0.2e1 + t138;
t34 = t198 / 0.2e1 + t137 * t117;
t22 = t128 - t136;
t19 = t127 + t133;
t18 = t126 - t132;
t11 = -t125 + t131;
t10 = -t124 + t130;
t26 = [0, 0, 0, 0, 0, t101 * qJD(2), t99 * qJD(2), t88 * qJD(2) + t120 * t147, t89 * qJD(2) - t117 * t147, t87 * qJD(3), t27 * qJD(2) - t31 * qJD(3), t39 * qJD(2) + t52 * qJD(3) - t83 * t146, -t40 * qJD(2) + t53 * qJD(3) + t81 * t146, t24 * qJD(2) - t36 * qJD(3) + t46 * qJD(4), t3 * qJD(2) - t4 * qJD(3) + t5 * qJD(4), -qJD(5) * t214, t14 * qJD(5), -t55 * t203, -t57 * t203, 0, t12 * qJD(2) + t20 * qJD(3) - t15 * qJD(4) + t2 * qJD(5), t13 * qJD(2) + t21 * qJD(3) + t25 * qJD(4) + t1 * qJD(5); 0, 0, 0, 0, 0, t158, t162, t164, t163, 0, t78 * qJD(3) + t177, t172, -t171, t179, t209 + t22 * qJD(3) + t34 * qJD(4) + (t116 * t80 + t119 * t85 - t192) * qJD(2) * t117, 0, 0, 0, 0, 0, t11 * qJD(5) + t194, t10 * qJD(5) + t187; 0, 0, 0, 0, 0, 0, 0, t140, -t141, t165, t78 * qJD(2) - t176, t168, t167, -t174, t22 * qJD(2) + t37 * qJD(4) - t208, 0, 0, 0, 0, 0, t19 * qJD(5) + t182, t18 * qJD(5) + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t143, t169, t34 * qJD(2) + t37 * qJD(3) + t207, 0, 0, 0, 0, 0, -t185, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t186, t144 * t55, t144 * t57, 0, t11 * qJD(2) + t19 * qJD(3) + t7 * qJD(5) + t210, t10 * qJD(2) + t18 * qJD(3) + t6 * qJD(5) + t213; 0, 0, 0, 0, 0, -t158, -t162, -t164, -t163, 0, -t77 * qJD(3) - t177, -t172, t171, -t179, t23 * qJD(3) + t35 * qJD(4) - t209, 0, 0, 0, 0, 0, -t8 * qJD(5) - t194, -t9 * qJD(5) - t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, 0, 0, 0, t180, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t135 - t206, qJD(5) * t134 - t205; 0, 0, 0, 0, 0, 0, 0, -t140, t141, -t165, t77 * qJD(2) + t176, -t168, -t167, t174, -t23 * qJD(2) + t38 * qJD(4) + t208, 0, 0, 0, 0, 0, -t16 * qJD(5) - t182, -t17 * qJD(5) - t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, 0, 0, -t180, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 * t159 - t184, t116 * t160 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, -t143, -t169, -t35 * qJD(2) - t38 * qJD(3) - t207, 0, 0, 0, 0, 0, -t42 * qJD(5) + t185, -t159 * t81 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 - t170, t144 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t186, t55 * t204, t57 * t204, 0, t8 * qJD(2) + t16 * qJD(3) + t42 * qJD(4) - t210, t123 * t81 * qJD(4) + t9 * qJD(2) + t17 * qJD(3) - t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t123 * t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t26;
