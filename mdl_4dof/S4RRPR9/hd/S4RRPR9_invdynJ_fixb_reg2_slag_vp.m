% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:07
% EndTime: 2019-12-31 17:10:12
% DurationCPUTime: 2.91s
% Computational Cost: add. (2387->335), mult. (5652->472), div. (0->0), fcn. (3847->10), ass. (0->171)
t129 = cos(qJ(2));
t119 = g(3) * t129;
t127 = sin(qJ(2));
t128 = sin(qJ(1));
t130 = cos(qJ(1));
t159 = g(1) * t130 + g(2) * t128;
t138 = t159 * t127 - t119;
t178 = t127 * qJDD(1);
t110 = pkin(5) * t178;
t181 = qJD(1) * qJD(2);
t168 = t129 * t181;
t69 = -qJDD(2) * pkin(2) + pkin(5) * t168 + qJDD(3) + t110;
t134 = t69 - t138;
t126 = sin(qJ(4));
t216 = cos(qJ(4));
t123 = sin(pkin(7));
t188 = qJD(1) * t127;
t171 = t123 * t188;
t124 = cos(pkin(7));
t185 = qJD(2) * t124;
t83 = -t171 + t185;
t173 = t124 * t188;
t186 = qJD(2) * t123;
t84 = t173 + t186;
t32 = t126 * t84 - t216 * t83;
t224 = t32 ^ 2;
t35 = t126 * t83 + t216 * t84;
t223 = t35 ^ 2;
t187 = qJD(1) * t129;
t106 = -qJD(4) + t187;
t222 = t106 * t32;
t170 = qJD(4) * t216;
t182 = qJD(4) * t126;
t221 = -t123 * t182 + t124 * t170;
t115 = t129 * qJDD(1);
t142 = t127 * t181 - t115;
t145 = -t126 * t123 + t216 * t124;
t121 = t127 ^ 2;
t122 = t129 ^ 2;
t190 = t121 - t122;
t162 = qJD(1) * t190;
t156 = pkin(2) * t127 - qJ(3) * t129;
t71 = qJD(2) * t156 - t127 * qJD(3);
t157 = pkin(2) * t129 + qJ(3) * t127;
t93 = -pkin(1) - t157;
t30 = qJD(1) * t71 + qJDD(1) * t93;
t60 = -pkin(5) * t142 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t17 = t123 * t30 + t124 * t60;
t143 = t168 + t178;
t179 = t124 * qJDD(2);
t52 = t123 * t143 - t179;
t11 = -pkin(6) * t52 + t17;
t77 = t93 * qJD(1);
t112 = pkin(5) * t187;
t98 = qJD(2) * qJ(3) + t112;
t38 = -t123 * t98 + t124 * t77;
t18 = -pkin(3) * t187 - pkin(6) * t84 + t38;
t39 = t123 * t77 + t124 * t98;
t19 = pkin(6) * t83 + t39;
t148 = t126 * t19 - t18 * t216;
t16 = -t123 * t60 + t124 * t30;
t180 = t123 * qJDD(2);
t53 = t124 * t143 + t180;
t8 = pkin(3) * t142 - t53 * pkin(6) + t16;
t1 = -qJD(4) * t148 + t11 * t216 + t126 * t8;
t220 = -0.2e1 * pkin(1);
t219 = pkin(3) * t52;
t218 = pkin(5) * t83;
t217 = pkin(5) * t84;
t215 = pkin(3) * t123;
t214 = g(1) * t128;
t211 = g(2) * t130;
t210 = g(3) * t127;
t209 = t35 * t32;
t208 = pkin(6) + qJ(3);
t197 = t124 * t129;
t150 = pkin(3) * t127 - pkin(6) * t197;
t90 = t156 * qJD(1);
t49 = pkin(5) * t171 + t124 * t90;
t28 = qJD(1) * t150 + t49;
t198 = t124 * t127;
t199 = t123 * t129;
t144 = -pkin(5) * t198 - pkin(6) * t199;
t75 = t123 * t90;
t36 = qJD(1) * t144 + t75;
t96 = t208 * t123;
t97 = t208 * t124;
t43 = -t126 * t97 - t216 * t96;
t207 = qJD(3) * t145 + qJD(4) * t43 - t126 * t28 - t216 * t36;
t44 = -t126 * t96 + t216 * t97;
t88 = t123 * t216 + t126 * t124;
t206 = -qJD(3) * t88 - qJD(4) * t44 + t126 * t36 - t216 * t28;
t139 = t129 * t88;
t74 = t88 * qJD(4);
t205 = qJD(1) * t139 - t74;
t204 = t145 * t187 - t221;
t203 = t123 * t53;
t202 = t124 * t52;
t184 = qJD(2) * t127;
t177 = pkin(5) * t184;
t40 = t123 * t177 + t124 * t71;
t55 = pkin(5) * t197 + t123 * t93;
t201 = pkin(5) * qJDD(1);
t132 = qJD(1) ^ 2;
t200 = t122 * t132;
t195 = t127 * t130;
t194 = t128 * t129;
t193 = t129 * t130;
t92 = -qJD(2) * pkin(2) + pkin(5) * t188 + qJD(3);
t192 = qJD(3) - t92;
t191 = t130 * pkin(1) + t128 * pkin(5);
t189 = t121 + t122;
t183 = qJD(2) * t129;
t176 = pkin(5) + t215;
t174 = t123 * t187;
t167 = t123 * t178;
t166 = t123 * t115;
t165 = t124 * t178;
t164 = t124 * t115;
t163 = t126 * t53 + t216 * t52;
t160 = t127 * t168;
t158 = -t211 + t214;
t155 = t123 * t84 - t124 * t83;
t109 = pkin(3) * t124 + pkin(2);
t154 = t109 * t129 + t127 * t208;
t152 = qJD(1) * (-t84 + t186);
t151 = qJD(1) * (-t83 + t185);
t82 = t124 * t93;
t37 = -pkin(6) * t198 + t82 + (-pkin(5) * t123 - pkin(3)) * t129;
t42 = -pkin(6) * t123 * t127 + t55;
t14 = -t126 * t42 + t216 * t37;
t6 = t126 * t18 + t19 * t216;
t15 = t126 * t37 + t216 * t42;
t146 = -pkin(5) * qJDD(2) + t181 * t220;
t9 = t126 * t52 - t83 * t170 + t182 * t84 - t216 * t53;
t141 = pkin(1) * t132 + t159;
t131 = qJD(2) ^ 2;
t140 = pkin(5) * t131 + qJDD(1) * t220 + t211;
t136 = -t129 * t159 - t210;
t2 = -qJD(4) * t6 - t126 * t11 + t216 * t8;
t10 = qJD(4) * t35 + t163;
t120 = pkin(7) + qJ(4);
t117 = t130 * pkin(5);
t114 = cos(t120);
t113 = sin(t120);
t107 = t127 * t214;
t103 = t127 * t132 * t129;
t91 = t176 * t127;
t86 = qJDD(4) + t142;
t80 = qJDD(1) * t122 - 0.2e1 * t160;
t79 = t176 * t183;
t78 = pkin(3) * t174 + t112;
t68 = t145 * t127;
t67 = t88 * t127;
t64 = t113 * t128 + t114 * t193;
t63 = -t113 * t193 + t114 * t128;
t62 = t113 * t130 - t114 * t194;
t61 = t113 * t194 + t114 * t130;
t58 = t123 * t71;
t54 = -pkin(5) * t199 + t82;
t50 = -pkin(5) * t173 + t75;
t48 = -pkin(3) * t83 + t92;
t41 = -t124 * t177 + t58;
t29 = qJD(2) * t144 + t58;
t25 = qJD(2) * t139 + t221 * t127;
t24 = t127 * t74 - t145 * t183;
t23 = t69 + t219;
t22 = qJD(2) * t150 + t40;
t4 = -qJD(4) * t15 - t126 * t29 + t216 * t22;
t3 = qJD(4) * t14 + t126 * t22 + t216 * t29;
t5 = [0, 0, 0, 0, 0, qJDD(1), t158, t159, 0, 0, qJDD(1) * t121 + 0.2e1 * t160, -0.2e1 * qJD(2) * t162 + 0.2e1 * t115 * t127, qJDD(2) * t127 + t129 * t131, t80, qJDD(2) * t129 - t127 * t131, 0, t146 * t127 + (-t140 + t214) * t129, t127 * t140 + t129 * t146 - t107, 0.2e1 * t189 * t201 - t159, -g(1) * (-pkin(1) * t128 + t117) - g(2) * t191 + (t189 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), (t127 * t53 + t183 * t84) * t124, (-t202 - t203) * t127 - t155 * t183, (-t53 - t165) * t129 + (t124 * t162 + t127 * t84) * qJD(2), (t127 * t52 - t183 * t83) * t123, (t52 + t167) * t129 + (-t123 * t162 + t127 * t83) * qJD(2), t80, -t159 * t123 + (pkin(5) * t52 + t69 * t123 + (qJD(1) * t54 + t38) * qJD(2)) * t127 + (-t40 * qJD(1) - t54 * qJDD(1) - t16 + t158 * t124 + (t123 * t92 - t218) * qJD(2)) * t129, -t159 * t124 + (pkin(5) * t53 + t69 * t124 + (-qJD(1) * t55 - t39) * qJD(2)) * t127 + (t41 * qJD(1) + t55 * qJDD(1) + t17 - t158 * t123 + (t124 * t92 + t217) * qJD(2)) * t129, -t40 * t84 + t41 * t83 - t52 * t55 - t53 * t54 + t107 + (-t123 * t39 - t124 * t38) * t183 + (-t123 * t17 - t124 * t16 - t211) * t127, t17 * t55 + t39 * t41 + t16 * t54 + t38 * t40 - g(1) * t117 - g(2) * (t130 * t157 + t191) - t93 * t214 + (t127 * t69 + t183 * t92) * pkin(5), -t24 * t35 - t68 * t9, -t10 * t68 + t24 * t32 - t25 * t35 + t67 * t9, t106 * t24 + t129 * t9 + t184 * t35 + t68 * t86, t10 * t67 + t25 * t32, t10 * t129 + t106 * t25 - t184 * t32 - t67 * t86, -t106 * t184 - t129 * t86, -g(1) * t62 - g(2) * t64 + t10 * t91 - t106 * t4 - t129 * t2 + t14 * t86 - t148 * t184 + t23 * t67 + t25 * t48 + t32 * t79, -g(1) * t61 - g(2) * t63 + t1 * t129 + t106 * t3 - t15 * t86 - t184 * t6 + t23 * t68 - t24 * t48 + t35 * t79 - t9 * t91, -g(2) * t195 - t1 * t67 - t10 * t15 + t14 * t9 - t148 * t24 - t2 * t68 - t25 * t6 - t3 * t32 - t35 * t4 + t107, t1 * t15 + t6 * t3 + t2 * t14 - t148 * t4 + t23 * t91 + t48 * t79 - g(1) * (t130 * t215 + t117) - g(2) * (t109 * t193 + t195 * t208 + t191) + (-g(1) * (-pkin(1) - t154) - g(2) * t215) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t190 * t132, t178, t103, t115, qJDD(2), t127 * t141 - t110 - t119, t210 + (t141 - t201) * t129, 0, 0, -t124 * t84 * t187 + t203, -t123 * t52 + t53 * t124 + t155 * t187, t124 * t200 + t127 * t152 - t166, t174 * t83 - t202, -t123 * t200 + t127 * t151 - t164, t103, qJ(3) * t166 - pkin(2) * t52 - t134 * t124 + ((-qJ(3) * t186 - t38) * t127 + (t123 * t192 + t218 + t49) * t129) * qJD(1), qJ(3) * t164 - pkin(2) * t53 + t134 * t123 + ((-qJ(3) * t185 + t39) * t127 + (t124 * t192 - t217 - t50) * t129) * qJD(1), t49 * t84 - t50 * t83 + (-qJ(3) * t52 + qJD(3) * t83 + t187 * t38 + t17) * t124 + (qJ(3) * t53 + qJD(3) * t84 + t187 * t39 - t16) * t123 + t136, -t92 * t112 - t38 * t49 - t39 * t50 + (-t123 * t38 + t124 * t39) * qJD(3) - t134 * pkin(2) + (-t16 * t123 + t17 * t124 + t136) * qJ(3), -t204 * t35 - t9 * t88, -t88 * t10 - t145 * t9 + t204 * t32 + t205 * t35, t106 * t204 - t188 * t35 + t88 * t86, -t10 * t145 - t205 * t32, -t106 * t205 + t145 * t86 + t188 * t32, t106 * t188, -t10 * t109 - t106 * t206 + t114 * t138 - t145 * t23 + t148 * t188 - t205 * t48 - t32 * t78 + t43 * t86, t106 * t207 + t109 * t9 - t113 * t138 + t188 * t6 - t204 * t48 + t23 * t88 - t35 * t78 - t44 * t86, t1 * t145 - t10 * t44 - t148 * t204 - t2 * t88 + t205 * t6 - t206 * t35 - t207 * t32 + t43 * t9 + t136, -g(3) * t154 + t1 * t44 - t23 * t109 + t2 * t43 - t206 * t148 + t207 * t6 - t48 * t78 + t159 * (t109 * t127 - t129 * t208); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t152 + t167 - t179, t129 * t151 + t165 + t180, -t83 ^ 2 - t84 ^ 2, t38 * t84 - t39 * t83 + t134, 0, 0, 0, 0, 0, 0, -t35 * t106 + t10, -t9 + t222, -t223 - t224, -t148 * t35 + t32 * t6 + t134 + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, t223 - t224, -t9 - t222, -t209, -t163 + (-qJD(4) - t106) * t35, t86, -g(1) * t63 + g(2) * t61 - t6 * t106 + t113 * t210 - t48 * t35 + t2, g(1) * t64 - g(2) * t62 + t106 * t148 + t114 * t210 + t48 * t32 - t1, 0, 0;];
tau_reg = t5;
