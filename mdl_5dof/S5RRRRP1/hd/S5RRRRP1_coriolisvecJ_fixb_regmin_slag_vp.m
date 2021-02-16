% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP1
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
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:44
% EndTime: 2021-01-15 23:52:51
% DurationCPUTime: 1.93s
% Computational Cost: add. (3747->251), mult. (9889->330), div. (0->0), fcn. (7141->6), ass. (0->169)
t128 = sin(qJ(4));
t125 = qJD(2) + qJD(3);
t130 = sin(qJ(2));
t213 = pkin(6) + pkin(7);
t111 = t213 * t130;
t105 = qJD(1) * t111;
t197 = qJD(2) * pkin(2);
t100 = -t105 + t197;
t209 = cos(qJ(3));
t131 = cos(qJ(2));
t112 = t213 * t131;
t107 = qJD(1) * t112;
t129 = sin(qJ(3));
t94 = t129 * t107;
t165 = t209 * t100 - t94;
t177 = qJD(1) * t209;
t190 = t129 * t131;
t93 = -qJD(1) * t190 - t130 * t177;
t88 = t93 * pkin(8);
t46 = t165 + t88;
t35 = t125 * pkin(3) + t46;
t208 = cos(qJ(4));
t98 = t209 * t107;
t154 = -t129 * t100 - t98;
t184 = qJD(1) * t130;
t92 = -t129 * t184 + t131 * t177;
t211 = t92 * pkin(8);
t47 = -t154 + t211;
t41 = t208 * t47;
t160 = -t128 * t35 - t41;
t104 = t209 * t130 + t190;
t138 = t125 * t104;
t136 = t138 * qJD(1);
t180 = qJD(2) * t213;
t163 = qJD(1) * t180;
t101 = t130 * t163;
t102 = t131 * t163;
t176 = t209 * qJD(3);
t183 = qJD(3) * t129;
t144 = t100 * t176 - t209 * t101 - t129 * t102 - t107 * t183;
t22 = -pkin(8) * t136 + t144;
t167 = t129 * t101 - t209 * t102;
t146 = t154 * qJD(3) + t167;
t71 = t92 * t125;
t23 = -t71 * pkin(8) + t146;
t172 = -t128 * t22 + t208 * t23;
t150 = t160 * qJD(4) + t172;
t175 = t208 * qJD(4);
t182 = qJD(4) * t128;
t156 = -t128 * t136 + t92 * t175 + t93 * t182 + t208 * t71;
t195 = t156 * qJ(5);
t142 = t150 - t195;
t158 = -t128 * t92 + t208 * t93;
t186 = t158 * qJD(5);
t2 = t142 + t186;
t65 = t128 * t93 + t208 * t92;
t173 = -pkin(4) * t65 + qJD(5);
t121 = -t131 * pkin(2) - pkin(1);
t110 = t121 * qJD(1);
t78 = -t92 * pkin(3) + t110;
t31 = t173 + t78;
t206 = t31 * t158;
t217 = t2 + t206;
t214 = t158 ^ 2;
t61 = t65 ^ 2;
t11 = -t61 + t214;
t134 = t208 * t138;
t149 = -qJD(1) * t134 + t158 * qJD(4) - t128 * t71;
t124 = qJD(4) + t125;
t192 = t158 * t124;
t7 = t149 - t192;
t204 = t78 * t158;
t216 = t150 + t204;
t194 = t65 * t124;
t6 = t156 - t194;
t205 = t158 * t65;
t193 = t65 * qJ(5);
t55 = t158 * qJ(5);
t147 = t128 * t23 + t35 * t175 - t47 * t182 + t208 * t22;
t143 = -t78 * t65 - t147;
t181 = qJD(2) * qJD(1);
t215 = -0.2e1 * t181;
t39 = t128 * t47;
t171 = t208 * t35 - t39;
t9 = t171 + t55;
t8 = t124 * pkin(4) + t9;
t212 = t8 - t9;
t210 = t93 * pkin(3);
t207 = pkin(3) * t124;
t203 = t93 * t92;
t166 = t129 * t105 - t98;
t49 = t166 - t211;
t198 = -t209 * t105 - t94;
t50 = t88 + t198;
t169 = -t128 * t50 + t208 * t49;
t14 = t169 - t193;
t120 = t209 * pkin(2) + pkin(3);
t179 = t208 * t129;
t70 = -t120 * t182 + (-t129 * t175 + (-t209 * t128 - t179) * qJD(3)) * pkin(2);
t202 = t14 - t70;
t199 = t128 * t49 + t208 * t50;
t15 = t55 + t199;
t191 = t128 * t129;
t69 = t120 * t175 + (-t129 * t182 + (t209 * t208 - t191) * qJD(3)) * pkin(2);
t201 = t15 - t69;
t200 = t208 * t46 - t39;
t196 = t110 * t93;
t133 = qJD(1) ^ 2;
t189 = t131 * t133;
t132 = qJD(2) ^ 2;
t188 = t132 * t130;
t187 = t132 * t131;
t185 = t130 ^ 2 - t131 ^ 2;
t123 = t130 * t197;
t122 = pkin(2) * t184;
t174 = t130 * t181;
t170 = -t128 * t46 - t41;
t164 = pkin(1) * t215;
t118 = pkin(2) * t174;
t34 = -pkin(4) * t158 - t210;
t10 = -t160 + t193;
t161 = -t10 * t158 + t8 * t65;
t56 = -t104 * pkin(8) - t209 * t111 - t129 * t112;
t153 = t129 * t111 - t209 * t112;
t155 = t129 * t130 - t209 * t131;
t57 = -t155 * pkin(8) - t153;
t159 = -t128 * t56 - t208 * t57;
t106 = t130 * t180;
t108 = t131 * t180;
t152 = -t209 * t106 - t129 * t108 - t111 * t176 - t112 * t183;
t29 = -t138 * pkin(8) + t152;
t145 = t153 * qJD(3) + t129 * t106 - t209 * t108;
t77 = t125 * t155;
t30 = t77 * pkin(8) + t145;
t157 = t128 * t30 + t56 * t175 - t57 * t182 + t208 * t29;
t151 = t208 * t155;
t81 = t155 * pkin(3) + t121;
t148 = t159 * qJD(4) - t128 * t29 + t208 * t30;
t76 = t208 * t104 - t128 * t155;
t141 = qJ(5) * t149 + t147;
t140 = (-t41 + (-t35 - t207) * t128) * qJD(4) + t172;
t139 = -t110 * t92 - t144;
t1 = qJD(5) * t65 + t141;
t135 = -t31 * t65 - t1;
t68 = t138 * pkin(3) + t123;
t51 = pkin(3) * t136 + t118;
t5 = -pkin(4) * t149 + t51;
t119 = t208 * pkin(3) + pkin(4);
t109 = t175 * t207;
t89 = pkin(2) * t179 + t128 * t120;
t86 = -pkin(2) * t191 + t208 * t120 + pkin(4);
t79 = t122 - t210;
t75 = t128 * t104 + t151;
t54 = t70 * t124;
t53 = t69 * t124;
t48 = -t92 ^ 2 + t93 ^ 2;
t42 = t75 * pkin(4) + t81;
t38 = -t93 * t125 - t136;
t32 = t122 + t34;
t25 = t76 * qJD(4) - t128 * t77 + t134;
t24 = qJD(4) * t151 + t104 * t182 + t128 * t138 + t208 * t77;
t18 = -t75 * qJ(5) - t159;
t17 = -t76 * qJ(5) - t128 * t57 + t208 * t56;
t16 = t25 * pkin(4) + t68;
t13 = t55 + t200;
t12 = t170 - t193;
t4 = t24 * qJ(5) - t76 * qJD(5) + t148;
t3 = -t25 * qJ(5) - t75 * qJD(5) + t157;
t19 = [0, 0, 0, 0.2e1 * t131 * t174, t185 * t215, t187, -t188, 0, -pkin(6) * t187 + t130 * t164, pkin(6) * t188 + t131 * t164, t71 * t104 + t93 * t77, -t104 * t136 + t93 * t138 - t71 * t155 - t77 * t92, -t77 * t125, -t138 * t125, 0, t155 * t118 - t92 * t123 + (0.2e1 * t104 * t110 + t145) * t125, t121 * t71 - t110 * t77 - t152 * t125 + (qJD(1) * t104 - t93) * t123, t156 * t76 + t158 * t24, t149 * t76 - t156 * t75 + t158 * t25 - t24 * t65, -t24 * t124, -t25 * t124, 0, t148 * t124 - t149 * t81 + t78 * t25 + t51 * t75 - t65 * t68, -t157 * t124 + t156 * t81 - t158 * t68 - t78 * t24 + t51 * t76, t4 * t124 - t149 * t42 - t16 * t65 + t31 * t25 + t5 * t75, -t3 * t124 + t156 * t42 - t158 * t16 - t31 * t24 + t5 * t76, -t1 * t75 - t10 * t25 + t149 * t18 - t156 * t17 + t158 * t4 - t2 * t76 + t8 * t24 + t3 * t65, t1 * t18 + t10 * t3 + t31 * t16 + t2 * t17 + t8 * t4 + t5 * t42; 0, 0, 0, -t130 * t189, t185 * t133, 0, 0, 0, t133 * pkin(1) * t130, pkin(1) * t189, t203, t48, 0, t38, 0, t92 * t122 + t196 - t166 * t125 + (-t98 + (-pkin(2) * t125 - t100) * t129) * qJD(3) + t167, t198 * t125 + (-t125 * t176 + t93 * t184) * pkin(2) + t139, t205, t11, t6, t7, 0, -t169 * t124 + t65 * t79 + t216 + t54, t199 * t124 + t158 * t79 + t143 - t53, -t14 * t124 + t32 * t65 + t217 + t54, t15 * t124 + t158 * t32 + t135 - t53, t149 * t89 - t156 * t86 - t158 * t202 - t201 * t65 + t161, t1 * t89 - t10 * t201 + t2 * t86 - t202 * t8 - t31 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t48, 0, t38, 0, -t154 * t125 + t146 + t196, t165 * t125 + t139, t205, t11, t6, t7, 0, -t170 * t124 - t210 * t65 + t140 + t204, t200 * t124 - t158 * t210 - t109 + t143, -t12 * t124 + t34 * t65 + t140 + t186 - t195 + t206, t13 * t124 + t158 * t34 - t109 + t135, -t119 * t156 - t12 * t158 - t13 * t65 + (t128 * t149 + (-t128 * t158 + t208 * t65) * qJD(4)) * pkin(3) + t161, -t10 * t13 + t2 * t119 - t8 * t12 - t31 * t34 + (t1 * t128 + (t10 * t208 - t128 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, t11, t6, t7, 0, -t160 * t124 + t216, t171 * t124 + t143, t10 * t124 - (-t173 - t31) * t158 + t142, -t214 * pkin(4) + t9 * t124 - (qJD(5) + t31) * t65 - t141, -pkin(4) * t156 + t212 * t65, t217 * pkin(4) + t212 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149 - t192, t156 + t194, -t61 - t214, -t10 * t65 - t158 * t8 + t5;];
tauc_reg = t19;
