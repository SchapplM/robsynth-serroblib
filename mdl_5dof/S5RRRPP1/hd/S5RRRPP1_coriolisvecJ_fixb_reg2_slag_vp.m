% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:58
% EndTime: 2019-12-31 20:50:02
% DurationCPUTime: 1.21s
% Computational Cost: add. (2289->227), mult. (4098->287), div. (0->0), fcn. (2571->6), ass. (0->154)
t122 = sin(pkin(8));
t123 = sin(qJ(3));
t125 = cos(qJ(3));
t193 = -qJ(4) - pkin(7);
t154 = qJD(3) * t193;
t133 = -t123 * qJD(4) + t125 * t154;
t184 = cos(pkin(8));
t152 = t184 * t125;
t181 = t122 * t123;
t138 = t152 - t181;
t126 = cos(qJ(2));
t188 = pkin(1) * qJD(1);
t163 = t126 * t188;
t115 = t125 * qJD(4);
t82 = t123 * t154 + t115;
t191 = t122 * t133 - t138 * t163 + t184 * t82;
t120 = t123 ^ 2;
t121 = t125 ^ 2;
t172 = t120 + t121;
t119 = qJD(1) + qJD(2);
t153 = t184 * t123;
t90 = t122 * t125 + t153;
t206 = t90 * t119;
t203 = t206 ^ 2;
t144 = t119 * t152;
t73 = t119 * t181 - t144;
t68 = t73 ^ 2;
t205 = -t68 - t203;
t204 = -t68 + t203;
t84 = t90 * qJD(3);
t66 = t119 * t84;
t143 = qJD(3) * t152;
t171 = qJD(3) * t123;
t157 = t122 * t171;
t97 = t119 * t157;
t67 = t119 * t143 - t97;
t124 = sin(qJ(2));
t187 = pkin(1) * qJD(2);
t162 = qJD(1) * t187;
t106 = t124 * t162;
t160 = t119 * t171;
t83 = pkin(3) * t160 + t106;
t137 = t66 * pkin(4) - t67 * qJ(5) + t83;
t13 = -qJD(5) * t206 + t137;
t112 = -t125 * pkin(3) - pkin(2);
t72 = t112 * t119 + qJD(4) - t163;
t25 = t73 * pkin(4) - qJ(5) * t206 + t72;
t202 = -t13 * t138 + t25 * t84;
t85 = t143 - t157;
t201 = -t13 * t90 - t25 * t85;
t200 = pkin(3) * t123;
t147 = t126 * t162;
t135 = qJD(4) * t119 + t147;
t164 = t124 * t188;
t95 = t119 * pkin(7) + t164;
t151 = qJ(4) * t119 + t95;
t141 = qJD(3) * t151;
t128 = -t135 * t123 - t125 * t141;
t43 = -t123 * t141 + t135 * t125;
t10 = t122 * t43 - t184 * t128;
t110 = t124 * pkin(1) + pkin(7);
t176 = -qJ(4) - t110;
t150 = t176 * t123;
t117 = t125 * qJ(4);
t88 = t125 * t110 + t117;
t49 = t122 * t88 - t184 * t150;
t199 = t10 * t49;
t103 = t125 * pkin(7) + t117;
t59 = t122 * t103 - t193 * t153;
t198 = t10 * t59;
t197 = t10 * t90;
t196 = t126 * pkin(1);
t195 = t25 * t206;
t194 = t206 * t73;
t11 = t122 * t128 + t184 * t43;
t192 = t122 * t82 - t184 * t133 - t90 * t163;
t190 = -t138 * t83 + t72 * t84;
t189 = t72 * t85 + t83 * t90;
t63 = t151 * t125;
t57 = t184 * t63;
t62 = t151 * t123;
t61 = qJD(3) * pkin(3) - t62;
t32 = t122 * t61 + t57;
t186 = t122 * t63;
t170 = qJD(3) * t125;
t96 = -t119 * pkin(2) - t163;
t185 = t123 * t106 + t96 * t170;
t183 = qJD(3) * t84;
t182 = t119 * t123;
t180 = t124 * t125;
t127 = qJD(3) ^ 2;
t179 = t127 * t123;
t116 = t127 * t125;
t149 = qJD(3) * t176;
t165 = t126 * t187;
t129 = (-qJD(4) - t165) * t123 + t125 * t149;
t55 = t123 * t149 + t125 * t165 + t115;
t23 = t122 * t55 - t184 * t129;
t178 = t23 * qJD(3);
t24 = t122 * t129 + t184 * t55;
t177 = t24 * qJD(3);
t35 = -t184 * t62 - t186;
t175 = qJD(5) - t35;
t174 = t172 * t147;
t173 = t120 - t121;
t169 = qJD(3) * t126;
t168 = -qJD(1) - t119;
t167 = -qJD(2) + t119;
t166 = pkin(3) * t182;
t114 = t124 * t187;
t113 = pkin(3) * t171;
t118 = t119 ^ 2;
t161 = t123 * t118 * t125;
t159 = t119 * t170;
t158 = t123 * t169;
t31 = t184 * t61 - t186;
t28 = -qJD(3) * pkin(4) + qJD(5) - t31;
t29 = qJD(3) * qJ(5) + t32;
t9 = qJD(3) * qJD(5) + t11;
t156 = t138 * t9 + t28 * t85 - t29 * t84 + t197;
t155 = t11 * t138 - t31 * t85 - t32 * t84 + t197;
t148 = t172 * qJD(2);
t146 = t123 * t159;
t33 = t84 * pkin(4) - t85 * qJ(5) - t90 * qJD(5) + t113;
t145 = -t33 + t164;
t142 = -t138 * t66 + t73 * t84;
t34 = -t122 * t62 + t57;
t140 = t34 * qJD(3) - t10;
t50 = t122 * t150 + t184 * t88;
t139 = t206 * t23 - t24 * t73 + t49 * t67 - t50 * t66;
t136 = -t119 * t96 - t147;
t134 = -t124 * t182 + t125 * t169;
t51 = -pkin(4) * t138 - t90 * qJ(5) + t112;
t1 = -t138 * t67 + t206 * t84 + t90 * t66 + t85 * t73;
t60 = t184 * t103 + t193 * t181;
t131 = -t191 * t73 + t192 * t206 + t59 * t67 - t60 * t66;
t130 = 0.2e1 * t206 * qJD(3);
t111 = -pkin(2) - t196;
t109 = -t184 * pkin(3) - pkin(4);
t107 = t122 * pkin(3) + qJ(5);
t98 = t112 - t196;
t94 = -0.2e1 * t146;
t93 = 0.2e1 * t146;
t92 = t114 + t113;
t86 = t96 * t171;
t79 = t85 * qJD(3);
t77 = -0.2e1 * t173 * t119 * qJD(3);
t46 = t51 - t196;
t42 = -t97 + (t144 + t73) * qJD(3);
t41 = -t97 + (t144 - t73) * qJD(3);
t36 = pkin(4) * t206 + t73 * qJ(5) + t166;
t26 = t114 + t33;
t20 = t206 * t85 + t67 * t90;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 * t114 - t106, t168 * t165, 0, 0, t93, t77, t116, t94, -t179, 0, t111 * t160 - t110 * t116 + t86 + (t168 * t180 - t158) * t187, t110 * t179 + t111 * t159 - t134 * t187 + t185, t119 * t148 * t196 + t174, ((qJD(1) * t111 + t96) * t124 + (qJD(1) * t110 + t95) * t126 * t172) * t187, t20, -t1, t79, t142, -t183, 0, t98 * t66 + t92 * t73 - t178 + t190, t206 * t92 + t98 * t67 - t177 + t189, t139 + t155, t11 * t50 - t31 * t23 + t32 * t24 + t72 * t92 + t83 * t98 + t199, t20, t79, t1, 0, t183, t142, t26 * t73 + t46 * t66 - t178 + t202, t139 + t156, -t206 * t26 - t46 * t67 + t177 + t201, t13 * t46 + t28 * t23 + t29 * t24 + t25 * t26 + t9 * t50 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t164 - t106, t167 * t163, 0, 0, t93, t77, t116, t94, -t179, 0, -pkin(2) * t160 - pkin(7) * t116 + t86 + (t167 * t180 + t158) * t188, -pkin(2) * t159 + pkin(7) * t179 + t134 * t188 + t185, -t172 * t119 * t163 + t174, ((-pkin(2) * qJD(2) - t96) * t124 + (pkin(7) * t148 - t172 * t95) * t126) * t188, t20, -t1, t79, t142, -t183, 0, -t73 * t164 + t112 * t66 + (t73 * t200 - t192) * qJD(3) + t190, -t206 * t164 + t112 * t67 + (t200 * t206 - t191) * qJD(3) + t189, t131 + t155, t198 + t11 * t60 + t83 * t112 + (-t164 + t113) * t72 + t191 * t32 - t192 * t31, t20, t79, t1, 0, t183, t142, -t192 * qJD(3) - t145 * t73 + t51 * t66 + t202, t131 + t156, t191 * qJD(3) + t145 * t206 - t51 * t67 + t201, t13 * t51 - t145 * t25 + t191 * t29 + t192 * t28 + t9 * t60 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t173 * t118, 0, t161, 0, 0, t136 * t123, t136 * t125, 0, 0, t194, t204, t42, -t194, 0, 0, -t73 * t166 - t206 * t72 + t140, t35 * qJD(3) - t166 * t206 + t72 * t73 - t11, (t32 - t34) * t206 + (-t31 + t35) * t73 + (-t122 * t66 - t184 * t67) * pkin(3), t31 * t34 - t32 * t35 + (-t184 * t10 + t11 * t122 - t72 * t182) * pkin(3), t194, t42, -t204, 0, 0, -t194, -t36 * t73 + t140 - t195, -t107 * t66 + t109 * t67 + (t29 - t34) * t206 + (t28 - t175) * t73, -t25 * t73 + t36 * t206 + (0.2e1 * qJD(5) - t35) * qJD(3) + t11, t10 * t109 + t9 * t107 + t175 * t29 - t25 * t36 - t28 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t41, t205, t206 * t31 + t32 * t73 + t83, 0, 0, 0, 0, 0, 0, t130, t205, -t41, t29 * t73 + (-qJD(5) - t28) * t206 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t42, -t203 - t127, -t29 * qJD(3) + t10 + t195;];
tauc_reg = t2;
