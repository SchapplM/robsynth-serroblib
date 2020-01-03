% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:28
% EndTime: 2019-12-31 20:58:33
% DurationCPUTime: 1.54s
% Computational Cost: add. (2173->258), mult. (5677->286), div. (0->0), fcn. (3727->4), ass. (0->150)
t122 = -pkin(3) - pkin(4);
t120 = sin(qJ(2));
t119 = sin(qJ(3));
t121 = cos(qJ(2));
t170 = t119 * t121;
t194 = cos(qJ(3));
t89 = t194 * t120 + t170;
t172 = qJD(1) * t89;
t154 = t121 * pkin(2) + pkin(1);
t97 = t154 * qJD(1);
t205 = -t172 * qJ(4) - t97;
t152 = t194 * t121;
t143 = qJD(1) * t152;
t162 = qJD(1) * t120;
t151 = t119 * t162;
t73 = -t143 + t151;
t17 = t122 * t73 + qJD(5) - t205;
t196 = -pkin(7) - pkin(6);
t153 = qJD(2) * t196;
t145 = qJD(1) * t153;
t85 = t120 * t145;
t86 = t121 * t145;
t148 = t119 * t85 - t194 * t86;
t99 = t196 * t121;
t95 = qJD(1) * t99;
t82 = t194 * t95;
t185 = qJD(2) * pkin(2);
t98 = t196 * t120;
t93 = qJD(1) * t98;
t84 = t93 + t185;
t52 = t119 * t84 - t82;
t15 = t52 * qJD(3) + t148;
t116 = qJD(2) + qJD(3);
t171 = t119 * t120;
t139 = t116 * t171;
t164 = t116 * t143;
t46 = qJD(1) * t139 - t164;
t4 = t46 * qJ(5) - qJD(5) * t172 + t15;
t127 = t17 * t172 - t4;
t36 = t73 * pkin(3) + t205;
t190 = t36 * t172;
t206 = t15 + t190;
t204 = t116 * t172;
t160 = qJD(1) * qJD(2);
t203 = -0.2e1 * t160;
t150 = qJD(3) * t194;
t103 = pkin(2) * t150 + qJD(4);
t113 = t116 * qJD(4);
t201 = t103 * t116 + t113;
t181 = t119 * t95;
t51 = t194 * t84 + t181;
t200 = qJD(4) - t51;
t197 = t172 ^ 2;
t199 = -t116 ^ 2 - t197;
t198 = t73 ^ 2;
t35 = t197 - t198;
t195 = pkin(4) * t172;
t57 = -t119 * t99 - t194 * t98;
t193 = t15 * t57;
t191 = t36 * t73;
t189 = t172 * t73;
t188 = t97 * t73;
t187 = t97 * t172;
t48 = pkin(3) * t172 + t73 * qJ(4);
t54 = t194 * t93 + t181;
t58 = t119 * t98 - t194 * t99;
t56 = t116 * t89;
t47 = t56 * qJD(1);
t186 = qJ(4) * t47;
t109 = pkin(2) * t119 + qJ(4);
t184 = t109 * t47;
t183 = t116 * t56;
t182 = t116 * t73;
t161 = qJD(3) * t119;
t94 = t120 * t153;
t19 = t98 * t150 + t153 * t170 + t99 * t161 + t194 * t94;
t180 = t19 * t116;
t144 = qJD(2) * t152;
t20 = t58 * qJD(3) + t119 * t94 - t196 * t144;
t179 = t20 * t116;
t53 = t119 * t93 - t82;
t178 = t53 * t116;
t55 = -t121 * t150 + t139 - t144;
t50 = t55 * t116;
t177 = t73 * qJ(5);
t68 = t172 * qJ(5);
t33 = t68 + t54;
t174 = t103 - t33;
t173 = t103 - t54;
t124 = qJD(1) ^ 2;
t169 = t121 * t124;
t123 = qJD(2) ^ 2;
t168 = t123 * t120;
t167 = t123 * t121;
t29 = t68 + t51;
t166 = qJD(4) - t29;
t163 = t120 ^ 2 - t121 ^ 2;
t159 = pkin(2) * t162;
t158 = pkin(2) * t161;
t157 = t120 * t185;
t155 = t120 * t169;
t149 = t120 * t160;
t147 = -t119 * t86 - t84 * t150 - t95 * t161 - t194 * t85;
t146 = pkin(1) * t203;
t111 = -t194 * pkin(2) - pkin(3);
t142 = t121 * t149;
t32 = t53 + t177;
t141 = -t32 + t158;
t140 = -t53 + t158;
t88 = -t152 + t171;
t10 = t47 * t88 + t56 * t73;
t138 = t47 * qJ(5) + t73 * qJD(5) - t147;
t38 = t159 + t48;
t137 = t89 * qJ(4) + t154;
t136 = t116 * t51 + t147;
t135 = t54 * t116 + t147;
t8 = pkin(2) * t149 + t47 * pkin(3) + t46 * qJ(4) - qJD(4) * t172;
t134 = t17 * t73 + t138;
t30 = t52 + t177;
t133 = -t55 * qJ(4) + t89 * qJD(4) - t157;
t132 = -t116 * t151 + t164;
t131 = t15 * t89 + t172 * t20 - t19 * t73 - t46 * t57 - t58 * t47;
t130 = -t172 * t56 + t46 * t88 - t47 * t89 + t55 * t73;
t2 = -pkin(4) * t47 - t8;
t23 = t46 - t182;
t128 = t52 * t116 - t15;
t114 = t116 * qJ(4);
t108 = 0.2e1 * t113;
t107 = -pkin(4) + t111;
t100 = t116 * t158;
t49 = pkin(3) * t88 - t137;
t41 = t114 + t52;
t40 = qJ(5) * t88 + t58;
t39 = -qJ(5) * t89 + t57;
t37 = -pkin(3) * t116 + t200;
t31 = t122 * t88 + t137;
t28 = -t48 - t195;
t25 = t47 - t204;
t24 = t132 + t182;
t22 = t114 + t30;
t21 = -t38 - t195;
t16 = t122 * t116 + t200 - t68;
t12 = t113 - t147;
t11 = pkin(3) * t56 - t133;
t9 = -t172 * t55 - t46 * t89;
t7 = t55 * qJ(5) - t89 * qJD(5) + t20;
t6 = qJ(5) * t56 + qJD(5) * t88 + t19;
t5 = t122 * t56 + t133;
t3 = t113 + t138;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t142, t163 * t203, t167, -0.2e1 * t142, -t168, 0, -pkin(6) * t167 + t120 * t146, pkin(6) * t168 + t121 * t146, 0, 0, t9, t130, -t50, t10, -t183, 0, -t154 * t47 - t179 - t97 * t56 + (qJD(1) * t88 + t73) * t157, t154 * t46 + 0.2e1 * t172 * t157 + t97 * t55 - t180, t147 * t88 + t51 * t55 - t52 * t56 + t131, -t147 * t58 - 0.2e1 * t97 * t157 + t52 * t19 - t51 * t20 + t193, t9, -t50, -t130, 0, t183, t10, t11 * t73 + t36 * t56 + t47 * t49 + t8 * t88 - t179, -t12 * t88 - t37 * t55 - t41 * t56 + t131, -t11 * t172 + t36 * t55 + t46 * t49 - t8 * t89 + t180, t11 * t36 + t12 * t58 + t19 * t41 + t20 * t37 + t49 * t8 + t193, t9, -t130, t50, t10, -t183, 0, -t116 * t7 - t17 * t56 - t2 * t88 - t31 * t47 - t5 * t73, t116 * t6 - t17 * t55 + t172 * t5 + t2 * t89 - t31 * t46, t16 * t55 - t172 * t7 + t22 * t56 + t3 * t88 + t39 * t46 - t4 * t89 + t40 * t47 + t6 * t73, t16 * t7 + t17 * t5 + t2 * t31 + t22 * t6 + t3 * t40 + t39 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t163 * t124, 0, t155, 0, 0, t124 * pkin(1) * t120, pkin(1) * t169, 0, 0, t189, t35, t24, -t189, 0, 0, -t73 * t159 + t178 + t187 + (t82 + (-pkin(2) * t116 - t84) * t119) * qJD(3) - t148, -t188 + (-t116 * t150 - t162 * t172) * pkin(2) + t135, (t52 - t53) * t172 + (-t51 + t54) * t73 + (t194 * t46 - t119 * t47 + (t119 * t172 - t194 * t73) * qJD(3)) * pkin(2), t51 * t53 - t52 * t54 + (t97 * t162 - t194 * t15 - t119 * t147 + (-t119 * t51 + t194 * t52) * qJD(3)) * pkin(2), t189, t24, -t35, 0, t25, -t189, -t38 * t73 - t100 + t178 - t206, -t184 - t111 * t46 + (t140 + t41) * t172 + (t37 - t173) * t73, t172 * t38 - t135 - t191 + t201, t12 * t109 + t15 * t111 + t140 * t37 + t173 * t41 - t36 * t38, t189, -t35, t23, -t189, 0, 0, t32 * t116 + t21 * t73 - t100 + t127, -t116 * t33 - t172 * t21 + t134 + t201, t107 * t46 + t184 + (-t141 - t22) * t172 + (-t16 + t174) * t73, t4 * t107 + t3 * t109 + t141 * t16 - t17 * t21 + t174 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t35, t24, -t189, 0, 0, t128 + t187, t136 - t188, 0, 0, t189, t24, -t35, 0, t25, -t189, -t48 * t73 + t128 - t190, pkin(3) * t46 - t186 + (t41 - t52) * t172 + (t37 - t200) * t73, t172 * t48 + t108 - t136 - t191, -t15 * pkin(3) + t12 * qJ(4) + t200 * t41 - t36 * t48 - t37 * t52, t189, -t35, t23, -t189, 0, 0, t30 * t116 + t28 * t73 + t127, -t116 * t29 - t172 * t28 + t108 + t134, t186 + t122 * t46 + (-t22 + t30) * t172 + (-t16 + t166) * t73, t3 * qJ(4) + t4 * t122 - t16 * t30 + t166 * t22 - t17 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t24, t199, -t41 * t116 + t206, 0, 0, 0, 0, 0, 0, t189, t199, t23, -t22 * t116 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t204, t132 - t182, -t197 - t198, t16 * t172 - t22 * t73 + t2;];
tauc_reg = t1;
