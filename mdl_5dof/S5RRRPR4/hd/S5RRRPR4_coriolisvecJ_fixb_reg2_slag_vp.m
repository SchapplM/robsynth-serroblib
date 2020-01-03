% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:45
% DurationCPUTime: 1.62s
% Computational Cost: add. (2284->244), mult. (3882->317), div. (0->0), fcn. (2155->6), ass. (0->164)
t126 = qJD(1) + qJD(2);
t132 = sin(qJ(3));
t129 = t132 ^ 2;
t135 = cos(qJ(3));
t130 = t135 ^ 2;
t182 = t129 + t130;
t225 = t126 * t182;
t133 = sin(qJ(2));
t198 = pkin(1) * qJD(1);
t172 = t133 * t198;
t92 = t126 * pkin(7) + t172;
t76 = t132 * t92;
t224 = qJD(4) + t76;
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t177 = qJD(5) * t134;
t178 = qJD(5) * t131;
t180 = qJD(3) * t135;
t181 = qJD(3) * t132;
t40 = t131 * t180 + t132 * t177 - t134 * t181 - t135 * t178;
t125 = qJD(3) - qJD(5);
t83 = t132 * t131 + t135 * t134;
t143 = t83 * qJD(5);
t166 = t126 * t180;
t167 = t126 * t181;
t22 = t126 * t143 - t131 * t167 - t134 * t166;
t63 = t83 * t126;
t219 = t63 * t125 + t22;
t189 = t126 * t132;
t223 = -pkin(8) * t189 + t224;
t136 = cos(qJ(2));
t197 = pkin(1) * qJD(2);
t168 = qJD(1) * t197;
t156 = t136 * t168;
t105 = t135 * t156;
t127 = qJD(3) * qJD(4);
t184 = t105 + t127;
t42 = -t92 * t181 + t184;
t103 = t132 * t156;
t47 = t92 * t180 + t103;
t221 = t47 * t132 + t42 * t135;
t213 = pkin(3) + pkin(4);
t39 = -t213 * qJD(3) + t223;
t128 = qJD(3) * qJ(4);
t188 = t126 * t135;
t77 = t135 * t92;
t55 = -pkin(8) * t188 + t77;
t46 = t128 + t55;
t15 = t131 * t39 + t134 * t46;
t26 = (pkin(8) * t126 - t92) * t181 + t184;
t29 = -pkin(8) * t166 + t47;
t5 = -t15 * qJD(5) - t131 * t26 + t134 * t29;
t220 = -t15 * t125 + t5;
t23 = t40 * t126;
t84 = -t135 * t131 + t132 * t134;
t65 = t84 * t126;
t218 = t65 * t125 + t23;
t216 = t136 * t182;
t121 = t132 * qJ(4);
t123 = t135 * pkin(3);
t215 = t121 + t123;
t183 = -t129 + t130;
t214 = 0.2e1 * t183 * t126 * qJD(3);
t212 = pkin(7) - pkin(8);
t138 = qJD(3) ^ 2;
t211 = pkin(7) * t138;
t161 = pkin(2) + t121;
t171 = t136 * t198;
t33 = t171 + (t213 * t135 + t161) * t126;
t210 = t33 * t65;
t209 = t65 * t63;
t114 = t133 * pkin(1) + pkin(7);
t208 = -pkin(8) + t114;
t118 = t132 * qJD(4);
t190 = qJ(4) * t135;
t146 = -t213 * t132 + t190;
t155 = t133 * t168;
t24 = -t155 + (qJD(3) * t146 + t118) * t126;
t207 = t24 * t83 + t33 * t40;
t41 = qJD(3) * t83 - t143;
t206 = t24 * t84 + t33 * t41;
t109 = t212 * t132;
t110 = t212 * t135;
t48 = t134 * t109 - t131 * t110;
t116 = pkin(8) * t181;
t85 = -pkin(7) * t181 + t116;
t86 = qJD(3) * t110;
t205 = qJD(5) * t48 + t131 * t86 + t134 * t85 - t83 * t171;
t49 = t131 * t109 + t134 * t110;
t204 = -qJD(5) * t49 - t131 * t85 + t134 * t86 - t84 * t171;
t95 = t134 * qJ(4) - t131 * t213;
t203 = qJD(5) * t95 + t223 * t131 + t134 * t55;
t94 = -t131 * qJ(4) - t134 * t213;
t202 = qJD(5) * t94 - t131 * t55 + t223 * t134;
t174 = t136 * t197;
t201 = t174 * t225;
t14 = -t131 * t46 + t134 * t39;
t196 = t14 * t125;
t179 = qJD(3) * t136;
t164 = t132 * t179;
t192 = t164 * t198 + t172 * t188;
t93 = -t126 * pkin(2) - t171;
t191 = t132 * t155 + t93 * t180;
t119 = t138 * t132;
t120 = t138 * t135;
t185 = t182 * t156;
t176 = -qJD(1) - t126;
t175 = pkin(2) + t215;
t173 = t133 * t197;
t4 = t131 * t29 + t134 * t26 + t39 * t177 - t46 * t178;
t170 = -t14 * t41 - t15 * t40 - t4 * t83 - t5 * t84;
t169 = t63 ^ 2 - t65 ^ 2;
t71 = pkin(3) * t181 - qJ(4) * t180 - t118;
t115 = -t136 * pkin(1) - pkin(2);
t150 = pkin(3) * t132 - t190;
t34 = t155 + (t150 * qJD(3) - t118) * t126;
t162 = -t34 - t211;
t82 = t208 * t135;
t159 = -qJD(3) * pkin(3) + qJD(4);
t158 = t125 ^ 2;
t154 = t132 * t166;
t56 = -pkin(4) * t181 - t71;
t153 = t56 + t172;
t152 = (-qJD(2) + t126) * t198;
t151 = t176 * t197;
t81 = t208 * t132;
t30 = -t131 * t82 + t134 * t81;
t31 = t131 * t81 + t134 * t82;
t60 = t159 + t76;
t67 = t77 + t128;
t149 = t132 * t60 + t135 * t67;
t78 = t115 - t215;
t148 = t126 * t78 - t174;
t57 = t71 + t173;
t147 = -t114 * t138 - t126 * t57 - t34;
t145 = t60 * t180 - t67 * t181 + t221;
t144 = t171 * t225;
t142 = -t133 * t189 + t135 * t179;
t140 = -t33 * t63 + t4;
t139 = (-t132 * t67 + t135 * t60) * qJD(3) + t221;
t124 = t126 ^ 2;
t122 = t135 * pkin(4);
t108 = t132 * t124 * t135;
t91 = -0.2e1 * t154;
t90 = 0.2e1 * t154;
t80 = t183 * t124;
t79 = t122 + t175;
t74 = t93 * t181;
t72 = t150 * t126;
t66 = t122 - t78;
t53 = t146 * t126;
t51 = qJD(3) * t82 + t132 * t174;
t50 = -t114 * t181 + t135 * t174 + t116;
t45 = -t171 + (-t161 - t123) * t126;
t44 = t56 - t173;
t38 = t45 * t181;
t36 = t41 * t125;
t35 = t40 * t125;
t11 = -t31 * qJD(5) - t131 * t50 + t134 * t51;
t10 = t30 * qJD(5) + t131 * t51 + t134 * t50;
t7 = t23 * t83 + t40 * t63;
t6 = -t22 * t84 + t41 * t65;
t1 = t22 * t83 - t23 * t84 - t40 * t65 - t41 * t63;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t151, t136 * t151, 0, 0, t90, t214, t120, t91, -t119, 0, t115 * t167 - t114 * t120 + t74 + (t176 * t135 * t133 - t164) * t197, t114 * t119 + t115 * t166 - t142 * t197 + t191, t185 + t201, ((qJD(1) * t115 + t93) * t133 + (qJD(1) * t114 + t92) * t216) * t197, t90, t120, -t214, 0, t119, t91, t135 * t147 + t148 * t181 + t38, t145 + t201, t147 * t132 + (-t148 - t45) * t180, t114 * t139 + t149 * t174 + t34 * t78 + t45 * t57, t6, t1, -t36, t7, t35, 0, -t11 * t125 + t23 * t66 + t44 * t63 + t207, t10 * t125 - t22 * t66 + t44 * t65 + t206, -t10 * t63 - t11 * t65 + t22 * t30 - t23 * t31 + t170, t10 * t15 + t11 * t14 + t24 * t66 + t30 * t5 + t31 * t4 + t33 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t152, t136 * t152, 0, 0, t90, t214, t120, t91, -t119, 0, -pkin(2) * t167 + t74 + (-t155 - t211) * t135 + t192, -pkin(2) * t166 + pkin(7) * t119 + t142 * t198 + t191, -t144 + t185, ((-pkin(2) * qJD(2) - t93) * t133 + (pkin(7) * qJD(2) - t92) * t216) * t198, t90, t120, -t214, 0, t119, t91, -t175 * t167 + t38 + (-t126 * t71 + t162) * t135 + t192, -t144 + t145, (t126 * t175 - t171 - t45) * t180 + ((-t71 + t172) * t126 + t162) * t132, -t34 * t175 + t45 * t71 + (-t133 * t45 - t136 * t149) * t198 + t139 * pkin(7), t6, t1, -t36, t7, t35, 0, -t204 * t125 + t153 * t63 + t79 * t23 + t207, t205 * t125 + t153 * t65 - t79 * t22 + t206, -t204 * t65 - t205 * t63 + t48 * t22 - t49 * t23 + t170, t204 * t14 + t205 * t15 + t153 * t33 + t24 * t79 + t4 * t49 + t5 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t80, 0, t108, 0, 0, -t93 * t189 - t103, -t93 * t188 - t105, 0, 0, -t108, 0, t80, 0, 0, t108, -t103 + (-t132 * t45 + t135 * t72) * t126, ((t67 - t128) * t132 + (t159 - t60) * t135) * t126, t105 + 0.2e1 * t127 + (t132 * t72 + t135 * t45) * t126, -t47 * pkin(3) + t42 * qJ(4) + t224 * t67 - t45 * t72 - t60 * t77, -t209, t169, t219, t209, t218, 0, t203 * t125 - t53 * t63 + t210 - t5, t202 * t125 - t53 * t65 + t140, t94 * t22 - t95 * t23 + (-t15 + t203) * t65 + (t14 - t202) * t63, -t203 * t14 + t202 * t15 - t33 * t53 + t4 * t95 + t5 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, -t124 * t129 - t138, -qJD(3) * t67 + t45 * t189 + t47, 0, 0, 0, 0, 0, 0, -t131 * t158 - t63 * t189, -t134 * t158 - t65 * t189, -t218 * t131 + t219 * t134, -t33 * t189 + t220 * t134 + (t4 + t196) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -t169, -t219, -t209, -t218, 0, -t210 + t220, -t140 - t196, 0, 0;];
tauc_reg = t2;
