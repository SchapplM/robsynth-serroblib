% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:32:02
% EndTime: 2019-03-09 08:32:07
% DurationCPUTime: 1.95s
% Computational Cost: add. (3261->296), mult. (8132->390), div. (0->0), fcn. (5667->6), ass. (0->167)
t126 = sin(pkin(9));
t127 = cos(pkin(9));
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t105 = t126 * t131 + t127 * t129;
t95 = t105 * qJD(1);
t208 = qJD(5) + t95;
t130 = cos(qJ(5));
t128 = sin(qJ(5));
t168 = t128 * qJD(2);
t167 = t129 * qJD(1);
t183 = t127 * t131;
t92 = -qJD(1) * t183 + t126 * t167;
t68 = -t130 * t92 + t168;
t159 = t208 * t68;
t94 = t105 * qJD(2);
t138 = qJD(1) * t94;
t169 = qJD(5) * t130;
t47 = qJD(5) * t168 - t128 * t138 - t92 * t169;
t213 = t47 - t159;
t165 = qJD(1) * qJD(2);
t163 = t129 * t165;
t116 = pkin(2) * t163;
t162 = t131 * t165;
t85 = -t126 * t163 + t127 * t162;
t158 = -t85 * qJ(4) + t116;
t178 = t95 * qJD(4);
t202 = pkin(3) + pkin(8);
t16 = t202 * t165 * t105 + t158 - t178;
t170 = qJD(5) * t128;
t193 = -qJ(3) - pkin(7);
t160 = qJD(2) * t193;
t89 = t131 * qJD(3) + t129 * t160;
t80 = t89 * qJD(1);
t90 = -t129 * qJD(3) + t131 * t160;
t81 = t90 * qJD(1);
t43 = t126 * t80 - t127 * t81;
t23 = t85 * pkin(4) + t43;
t164 = -t131 * pkin(2) - pkin(1);
t148 = t164 * qJD(1);
t109 = qJD(3) + t148;
t137 = -t95 * qJ(4) + t109;
t26 = t202 * t92 + t137;
t110 = t193 * t129;
t107 = qJD(1) * t110;
t103 = qJD(2) * pkin(2) + t107;
t111 = t193 * t131;
t108 = qJD(1) * t111;
t98 = t126 * t108;
t61 = t127 * t103 + t98;
t149 = qJD(4) - t61;
t199 = t95 * pkin(4);
t31 = -t202 * qJD(2) + t149 + t199;
t142 = t128 * t23 + t130 * t16 + t31 * t169 - t26 * t170;
t70 = t130 * qJD(2) + t128 * t92;
t78 = t130 * t138;
t48 = qJD(5) * t70 - t78;
t2 = -t48 * qJ(6) - t68 * qJD(6) + t142;
t11 = -t128 * t26 + t130 * t31;
t6 = -t70 * qJ(6) + t11;
t5 = pkin(5) * t208 + t6;
t206 = t208 * t5 - t2;
t12 = t128 * t31 + t130 * t26;
t161 = -t128 * t16 + t130 * t23;
t135 = -qJD(5) * t12 + t161;
t1 = t85 * pkin(5) + t47 * qJ(6) - t70 * qJD(6) + t135;
t7 = -t68 * qJ(6) + t12;
t207 = t208 * t7 + t1;
t215 = -t206 * t128 + t207 * t130;
t209 = t130 * t208;
t214 = t70 * t209;
t155 = t128 * t208;
t79 = t130 * t85;
t145 = -t208 * t155 + t79;
t212 = -0.2e1 * t165;
t91 = t95 ^ 2;
t210 = -t92 ^ 2 - t91;
t139 = t95 * qJD(2);
t64 = t127 * t107 + t98;
t176 = -qJD(4) + t64;
t119 = -t127 * pkin(2) - pkin(3);
t115 = -pkin(8) + t119;
t200 = t92 * pkin(4);
t184 = t127 * t108;
t62 = t126 * t103 - t184;
t59 = -qJD(2) * qJ(4) - t62;
t36 = -t59 - t200;
t205 = t115 * t85 + t208 * t36;
t204 = t70 ^ 2;
t201 = t5 - t6;
t166 = t130 * qJD(6);
t177 = qJ(6) - t115;
t157 = pkin(2) * t167 + t92 * qJ(4);
t32 = t202 * t95 + t157;
t63 = t126 * t107 - t184;
t45 = t63 - t200;
t40 = t130 * t45;
t198 = t177 * t170 - t166 + t92 * pkin(5) - t40 - (-qJ(6) * t95 - t32) * t128;
t102 = t177 * t130;
t182 = t130 * qJ(6);
t192 = t128 * t45 + t130 * t32;
t197 = -qJD(5) * t102 - t128 * qJD(6) - t95 * t182 - t192;
t65 = -t127 * t110 - t126 * t111;
t196 = t43 * t65;
t195 = t70 * t92;
t194 = t92 * t68;
t104 = t126 * t129 - t183;
t147 = -t105 * qJ(4) + t164;
t42 = t202 * t104 + t147;
t56 = t105 * pkin(4) + t65;
t191 = t128 * t56 + t130 * t42;
t44 = t126 * t81 + t127 * t80;
t189 = t128 * t85;
t188 = t128 * t94;
t187 = t130 * t47;
t186 = t130 * t94;
t185 = t36 * t128;
t133 = qJD(1) ^ 2;
t181 = t131 * t133;
t132 = qJD(2) ^ 2;
t180 = t132 * t129;
t179 = t132 * t131;
t175 = t199 - t176;
t173 = t129 ^ 2 - t131 ^ 2;
t172 = qJD(2) * t129;
t171 = qJD(5) * t115;
t122 = qJD(2) * qJD(4);
t38 = -t122 - t44;
t121 = pkin(2) * t172;
t117 = t126 * pkin(2) + qJ(4);
t54 = t126 * t89 - t127 * t90;
t156 = -qJ(6) * t104 - t42;
t152 = pkin(1) * t212;
t55 = t126 * t90 + t127 * t89;
t66 = t126 * t110 - t127 * t111;
t50 = t92 * pkin(3) + t137;
t146 = t50 * t95 + t43;
t144 = t104 * t170 - t186;
t97 = qJD(2) * t183 - t126 * t172;
t143 = -t97 * qJ(4) - t105 * qJD(4) + t121;
t19 = t202 * t94 + t143;
t33 = t97 * pkin(4) + t54;
t141 = t128 * t33 + t130 * t19 + t56 * t169 - t42 * t170;
t34 = -t94 * pkin(4) + t55;
t140 = -t208 * t209 - t189;
t136 = t43 * t105 + t54 * t95 - t55 * t92 + t65 * t85;
t20 = -pkin(4) * t138 - t38;
t134 = pkin(3) * t139 + t158;
t13 = t48 * pkin(5) + t20;
t101 = t177 * t128;
t87 = qJD(2) * t92;
t67 = t68 ^ 2;
t60 = t104 * pkin(3) + t147;
t58 = -qJD(2) * pkin(3) + t149;
t57 = -t104 * pkin(4) + t66;
t53 = t95 * pkin(3) + t157;
t52 = t130 * t56;
t41 = t130 * t48;
t35 = t94 * pkin(3) + t143;
t30 = t130 * t33;
t24 = t134 - t178;
t17 = pkin(5) * t68 + qJD(6) + t36;
t14 = t104 * t182 + t191;
t10 = t105 * pkin(5) + t156 * t128 + t52;
t4 = -qJ(6) * t144 + t104 * t166 + t141;
t3 = t97 * pkin(5) + t30 + t156 * t169 + (-qJ(6) * t94 - qJD(5) * t56 - qJD(6) * t104 - t19) * t128;
t8 = [0, 0, 0, 0.2e1 * t129 * t162, t173 * t212, t179, -t180, 0, -pkin(7) * t179 + t129 * t152, pkin(7) * t180 + t131 * t152, -t44 * t104 - t139 * t66 - t61 * t97 - t62 * t94 + t136, t196 + t44 * t66 - t61 * t54 + t62 * t55 + (t109 + t148) * t121, t38 * t104 - t138 * t66 + t58 * t97 + t59 * t94 + t136, -t24 * t104 - t35 * t92 - t50 * t94 + (-t60 * t95 + t54) * qJD(2), qJD(2) * t55 - t105 * t24 - t35 * t95 - t50 * t97 - t60 * t85, t24 * t60 + t35 * t50 - t38 * t66 + t54 * t58 - t55 * t59 + t196, t70 * t188 + (-t128 * t47 + t169 * t70) * t104 (-t128 * t68 + t130 * t70) * t94 + (-t128 * t48 - t187 + (-t128 * t70 - t130 * t68) * qJD(5)) * t104, t208 * t188 - t47 * t105 + t70 * t97 + (t169 * t208 + t189) * t104, t208 * t186 - t48 * t105 - t68 * t97 + (-t170 * t208 + t79) * t104, t105 * t85 + t208 * t97 (-t128 * t19 + t30) * t208 + (-t128 * t42 + t52) * t85 + t161 * t105 + t11 * t97 + t34 * t68 + t57 * t48 + (-t20 * t104 - t36 * t94) * t130 + (t104 * t185 - t12 * t105 - t191 * t208) * qJD(5), -t141 * t208 - t191 * t85 - t142 * t105 - t12 * t97 + t34 * t70 - t57 * t47 + t94 * t185 + (t20 * t128 + t169 * t36) * t104, t10 * t47 - t14 * t48 - t3 * t70 - t4 * t68 + (-t128 * t5 + t130 * t7) * t94 + (-t1 * t128 + t2 * t130 + (-t128 * t7 - t130 * t5) * qJD(5)) * t104, t2 * t14 + t7 * t4 + t1 * t10 + t5 * t3 + t13 * ((-pkin(5) * t130 - pkin(4)) * t104 + t66) + t17 * (pkin(5) * t144 + t34); 0, 0, 0, -t129 * t181, t173 * t133, 0, 0, 0, t133 * pkin(1) * t129, pkin(1) * t181 (t62 - t63) * t95 + (t64 - t61) * t92 + (-t126 * t139 - t127 * t85) * pkin(2), t61 * t63 - t62 * t64 + (-t109 * t167 + t126 * t44 - t127 * t43) * pkin(2), t119 * t85 + (-t59 - t63) * t95 - t117 * t138 + (t58 + t176) * t92, -t63 * qJD(2) + t53 * t92 + t146, -qJD(2) * t64 - t50 * t92 + t53 * t95 + 0.2e1 * t122 + t44, -t38 * t117 + t43 * t119 + t176 * t59 - t50 * t53 - t58 * t63, -t155 * t70 - t187, -t41 - t214 + (t47 + t159) * t128, t145 + t195, t140 - t194, t208 * t92, t11 * t92 + t117 * t48 - t40 * t208 + t175 * t68 + (t20 + (t32 - t171) * t208) * t128 + t205 * t130, -t117 * t47 + t192 * t208 - t12 * t92 + t175 * t70 + (-t171 * t208 + t20) * t130 - t205 * t128, t101 * t48 - t102 * t47 - t197 * t68 - t198 * t70 - t215, -t2 * t101 - t1 * t102 + t13 * (t128 * pkin(5) + t117) + t197 * t7 + t198 * t5 + (pkin(5) * t209 + t175) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t61 * t95 + t62 * t92 + t116, t210, -0.2e1 * t139, -t85 + t87, -t59 * t92 + (-qJD(4) - t58) * t95 + t134, 0, 0, 0, 0, 0, t140 + t194, t195 - t145, -t213 * t128 + t214 - t41, -t128 * t207 - t130 * t206 + t17 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 + t87, -t95 * t92, -t91 - t132, t59 * qJD(2) + t146, 0, 0, 0, 0, 0, -qJD(2) * t68 + t145, -qJD(2) * t70 + t140, t213 * t130 + (t208 * t70 - t48) * t128, -t17 * qJD(2) + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t68, -t67 + t204, -t213, t78 + (-qJD(5) + t208) * t70, t85, t12 * t208 - t36 * t70 + t135, t11 * t208 + t36 * t68 - t142, pkin(5) * t47 - t201 * t68, t201 * t7 + (-t17 * t70 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 - t204, t5 * t70 + t7 * t68 + t13;];
tauc_reg  = t8;
