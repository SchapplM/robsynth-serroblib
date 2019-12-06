% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:36
% EndTime: 2019-12-05 17:00:45
% DurationCPUTime: 2.62s
% Computational Cost: add. (2294->322), mult. (6032->437), div. (0->0), fcn. (4224->8), ass. (0->169)
t97 = sin(qJ(4));
t228 = pkin(7) * t97;
t100 = cos(qJ(4));
t179 = qJD(3) * t97;
t98 = sin(qJ(3));
t180 = qJD(2) * t98;
t74 = t100 * t180 + t179;
t101 = cos(qJ(3));
t172 = qJD(2) * t101;
t86 = -qJD(4) + t172;
t199 = t74 * t86;
t168 = t100 * qJD(3);
t72 = t97 * t180 - t168;
t202 = t72 * t86;
t167 = qJD(2) * qJD(3);
t146 = t101 * t167;
t177 = qJD(4) * t97;
t157 = t98 * t177;
t166 = qJD(3) * qJD(4);
t52 = qJD(2) * t157 + (-t146 - t166) * t100;
t170 = qJD(4) * t100;
t149 = t98 * t170;
t171 = qJD(3) * t101;
t112 = t97 * t171 + t149;
t53 = qJD(2) * t112 + t97 * t166;
t227 = (t53 - t199) * t97 + (t52 - t202) * t100;
t102 = cos(qJ(2));
t95 = sin(pkin(5));
t181 = qJD(2) * t95;
t150 = t102 * t181;
t96 = cos(pkin(5));
t111 = qJD(1) * (qJD(3) * t96 + t150);
t178 = qJD(3) * t98;
t183 = qJD(1) * t95;
t99 = sin(qJ(2));
t155 = t99 * t183;
t78 = qJD(2) * pkin(7) + t155;
t34 = t101 * t111 - t78 * t178;
t182 = qJD(1) * t96;
t153 = t98 * t182;
t58 = t101 * t78 + t153;
t49 = qJD(3) * pkin(8) + t58;
t135 = pkin(3) * t98 - pkin(8) * t101;
t77 = t135 * qJD(3);
t56 = (t77 + t155) * qJD(2);
t148 = t102 * t183;
t84 = -pkin(3) * t101 - pkin(8) * t98 - pkin(2);
t59 = t84 * qJD(2) - t148;
t145 = -t100 * t56 + t49 * t170 + t59 * t177 + t97 * t34;
t16 = t100 * t49 + t59 * t97;
t119 = -t16 * t86 - t145;
t184 = t53 * t100;
t203 = t52 * t97;
t224 = ((t100 * t74 - t72 * t97) * qJD(4) + t184 - t203) * t98 + (t100 * t72 + t74 * t97) * t171;
t187 = t101 * t97;
t223 = qJD(2) * ((t72 + t168) * t98 - t86 * t187) + t86 * t177;
t93 = t98 ^ 2;
t125 = qJD(2) * t93 - t101 * t86;
t222 = qJD(3) * (t125 * t97 + t72 * t98) - t101 * t53 - t86 * t149;
t13 = -qJ(5) * t86 + t16;
t90 = t98 * t167;
t143 = pkin(4) * t90;
t2 = -t143 + t145;
t221 = t13 * t86 + t2;
t169 = qJD(4) * t101;
t175 = t100 * t102;
t195 = t97 * t77 + t84 * t170 + (-t98 * t168 - t97 * t169) * pkin(7) - (t101 * t175 + t97 * t99) * t183;
t220 = t53 + t199;
t57 = t101 * t182 - t98 * t78;
t194 = t178 * t228 + t148 * t187 - t84 * t177 + (-t169 * pkin(7) - t155 + t77) * t100;
t215 = t74 ^ 2;
t214 = pkin(8) * t74;
t213 = pkin(8) * t86;
t35 = t111 * t98 + t78 * t171;
t5 = pkin(4) * t53 + qJ(5) * t52 - qJD(5) * t74 + t35;
t212 = t5 * t97;
t211 = t100 * t5;
t48 = -qJD(3) * pkin(3) - t57;
t17 = pkin(4) * t72 - qJ(5) * t74 + t48;
t208 = t17 * t74;
t198 = t95 * t99;
t65 = -t96 * t101 + t98 * t198;
t207 = t35 * t65;
t206 = t35 * t97;
t205 = t35 * t98;
t204 = t48 * t97;
t200 = t74 * t72;
t197 = qJ(5) * t178 - qJD(5) * t101 + t195;
t196 = -pkin(4) * t178 - t194;
t132 = pkin(4) * t97 - qJ(5) * t100;
t193 = t153 + (t132 * qJD(2) + t78) * t101 - t132 * qJD(4) + t97 * qJD(5);
t76 = t135 * qJD(2);
t29 = t100 * t57 + t97 * t76;
t176 = t100 * t101;
t62 = pkin(7) * t176 + t97 * t84;
t94 = t101 ^ 2;
t192 = t93 - t94;
t191 = qJD(2) * pkin(2);
t190 = t100 * t35;
t189 = t100 * t48;
t188 = t100 * t84;
t103 = qJD(3) ^ 2;
t186 = t103 * t98;
t104 = qJD(2) ^ 2;
t185 = t104 * t95;
t174 = t103 * t101;
t15 = t100 * t59 - t49 * t97;
t173 = qJD(5) - t15;
t165 = t97 * t213;
t164 = t100 * t213;
t163 = pkin(8) * t178;
t161 = t102 * t95 * t97;
t160 = t99 * t185;
t159 = t99 * t181;
t156 = t86 * t180;
t154 = t98 * t104 * t101;
t152 = t72 ^ 2 - t215;
t142 = t72 * t148;
t141 = t74 * t148;
t140 = t98 * t148;
t139 = t98 * t150;
t138 = t101 * t150;
t137 = qJ(5) * t90;
t136 = t98 * t146;
t79 = -t148 - t191;
t134 = -t79 - t148;
t133 = (qJD(4) * t72 - t52) * pkin(8);
t131 = pkin(4) * t100 + qJ(5) * t97;
t12 = pkin(4) * t86 + t173;
t130 = t100 * t12 - t13 * t97;
t129 = -t100 * t15 - t16 * t97;
t28 = t100 * t76 - t57 * t97;
t122 = pkin(7) + t132;
t66 = t101 * t198 + t96 * t98;
t121 = t101 * t17 + t163;
t120 = -t101 * t48 - t163;
t40 = t95 * t175 + t66 * t97;
t117 = -t100 * t34 - t59 * t170 + t49 * t177 - t97 * t56;
t41 = t100 * t66 - t161;
t38 = -qJD(3) * t65 + t138;
t8 = -qJD(4) * t161 - t100 * t159 + t66 * t170 + t38 * t97;
t9 = -qJD(4) * t40 + t38 * t100 + t97 * t159;
t116 = -t40 * t52 - t41 * t53 - t9 * t72 + t74 * t8;
t115 = qJD(2) * t134;
t39 = qJD(3) * t66 + t139;
t110 = t39 * t72 - t40 * t90 + t65 * t53 + t8 * t86;
t109 = qJD(3) * (-t134 - t191);
t108 = -t202 * t97 - t184;
t107 = -t39 * t74 + t41 * t90 + t52 * t65 - t86 * t9;
t106 = t101 * t34 + t205 + (-t101 * t57 - t58 * t98) * qJD(3);
t105 = t53 * t97 * t98 + t112 * t72;
t80 = -pkin(3) - t131;
t63 = t122 * t98;
t61 = -pkin(7) * t187 + t188;
t60 = (-t86 - t172) * t178;
t51 = -t188 + (pkin(4) + t228) * t101;
t50 = -qJ(5) * t101 + t62;
t42 = pkin(8) * t184;
t37 = pkin(4) * t74 + qJ(5) * t72;
t26 = -t52 - t202;
t23 = (t131 * qJD(4) - qJD(5) * t100) * t98 + t122 * t171;
t22 = -pkin(4) * t180 - t28;
t21 = qJ(5) * t180 + t29;
t20 = -t86 * t170 + (t86 * t176 + (-t74 + t179) * t98) * qJD(2);
t11 = -t100 * t199 - t203;
t10 = -t74 * t157 + (t74 * t171 - t52 * t98) * t100;
t7 = t86 * t157 + t52 * t101 + (t125 * t100 + t74 * t98) * qJD(3);
t1 = -qJD(5) * t86 - t117 + t137;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, -t102 * t185, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t160 + (-t39 - t139) * qJD(3), t98 * t160 + (-t38 - t138) * qJD(3), (t101 * t38 + t39 * t98 + (t101 * t65 - t66 * t98) * qJD(3)) * qJD(2), t34 * t66 + t207 + t38 * t58 - t39 * t57 + (t79 - t148) * t159, 0, 0, 0, 0, 0, 0, t110, -t107, t116, -t117 * t41 + t145 * t40 - t15 * t8 + t16 * t9 + t39 * t48 + t207, 0, 0, 0, 0, 0, 0, t110, t116, t107, t1 * t41 + t12 * t8 + t13 * t9 + t17 * t39 + t2 * t40 + t5 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, -0.2e1 * t192 * t167, t174, -0.2e1 * t136, -t186, 0, -pkin(7) * t174 + t109 * t98, pkin(7) * t186 + t101 * t109, (-t93 - t94) * qJD(2) * t148 + t106, ((-t79 - t191) * t99 + (-t101 * t58 + t57 * t98) * t102) * t183 + t106 * pkin(7), t10, -t224, t7, t105, -t222, t60, -t194 * t86 + (t145 + (pkin(7) * t72 + t204) * qJD(3)) * t101 + (-t142 + t48 * t170 + pkin(7) * t53 + t206 + (qJD(2) * t61 + t15) * qJD(3)) * t98, t195 * t86 + (-t117 + (pkin(7) * t74 + t189) * qJD(3)) * t101 + (-t141 - t48 * t177 - pkin(7) * t52 + t190 + (-qJD(2) * t62 - t16) * qJD(3)) * t98, t52 * t61 - t53 * t62 - t194 * t74 - t195 * t72 + t129 * t171 + (t100 * t145 + t117 * t97 + (-t100 * t16 + t15 * t97) * qJD(4)) * t98, -t48 * t140 - t117 * t62 - t145 * t61 + t195 * t16 + t194 * t15 + (t48 * t171 + t205) * pkin(7), t10, t7, t224, t60, t222, t105, t23 * t72 + t53 * t63 + t196 * t86 + (t17 * t179 + t2) * t101 + (-t142 + t17 * t170 + t212 + (-qJD(2) * t51 - t12) * qJD(3)) * t98, -t50 * t53 - t51 * t52 + t196 * t74 - t197 * t72 + t130 * t171 + (-t1 * t97 + t100 * t2 + (-t100 * t13 - t12 * t97) * qJD(4)) * t98, -t23 * t74 + t52 * t63 - t197 * t86 + (-t168 * t17 - t1) * t101 + (t141 + t17 * t177 - t211 + (qJD(2) * t50 + t13) * qJD(3)) * t98, t1 * t50 + t2 * t51 + t5 * t63 + (t23 - t140) * t17 + t197 * t13 + t196 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t192 * t104, 0, t154, 0, 0, t98 * t115, t101 * t115, 0, 0, t11, -t227, t20, t108, t223, t156, -pkin(3) * t53 - t190 + t28 * t86 - t58 * t72 + (t164 + t204) * qJD(4) + (t120 * t97 - t15 * t98) * qJD(2), pkin(3) * t52 - t29 * t86 + t206 - t58 * t74 + (-t165 + t189) * qJD(4) + (t100 * t120 + t16 * t98) * qJD(2), t28 * t74 + t29 * t72 - t42 + (t15 * t172 - t117 + (-t15 + t214) * qJD(4)) * t100 + (t133 - t119) * t97, -pkin(3) * t35 - t15 * t28 - t16 * t29 - t48 * t58 + (t129 * qJD(4) - t100 * t117 + t145 * t97) * pkin(8), t11, t20, t227, t156, -t223, t108, -t211 - t22 * t86 + t53 * t80 - t193 * t72 + (t17 * t97 + t164) * qJD(4) + (t12 * t98 - t121 * t97) * qJD(2), t21 * t72 - t22 * t74 - t42 + (-t12 * t172 + t1 + (t12 + t214) * qJD(4)) * t100 + (t133 + t221) * t97, t21 * t86 - t212 + t52 * t80 + t193 * t74 + (-t100 * t17 + t165) * qJD(4) + (t100 * t121 - t13 * t98) * qJD(2), -t12 * t22 - t13 * t21 + t5 * t80 - t193 * t17 + (qJD(4) * t130 + t1 * t100 + t2 * t97) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t152, t26, -t200, -t220, t90, -t48 * t74 + t119, -t15 * t86 + t48 * t72 + t117, 0, 0, t200, t26, t152, t90, t220, -t200, -t37 * t72 + t119 + 0.2e1 * t143 - t208, pkin(4) * t52 - qJ(5) * t53 + (t13 - t16) * t74 + (t12 - t173) * t72, 0.2e1 * t137 - t17 * t72 + t37 * t74 + (-0.2e1 * qJD(5) + t15) * t86 - t117, -pkin(4) * t2 + qJ(5) * t1 - t12 * t16 + t13 * t173 - t17 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 + t200, t26, -t86 ^ 2 - t215, t208 + t221;];
tauc_reg = t3;
