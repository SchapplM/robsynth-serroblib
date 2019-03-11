% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:21
% EndTime: 2019-03-09 01:49:29
% DurationCPUTime: 2.45s
% Computational Cost: add. (3279->323), mult. (6707->464), div. (0->0), fcn. (4165->6), ass. (0->177)
t128 = sin(qJ(6));
t130 = cos(qJ(6));
t124 = sin(pkin(9));
t131 = cos(qJ(4));
t188 = qJD(1) * t131;
t163 = t124 * t188;
t125 = cos(pkin(9));
t181 = t125 * qJD(4);
t84 = -t163 + t181;
t162 = t125 * t188;
t182 = t124 * qJD(4);
t85 = t162 + t182;
t36 = t128 * t85 - t130 * t84;
t228 = t36 ^ 2;
t39 = t128 * t84 + t130 * t85;
t227 = t39 ^ 2;
t129 = sin(qJ(4));
t180 = t129 * qJD(1);
t113 = qJD(6) + t180;
t226 = t36 * t113;
t126 = -pkin(7) + qJ(2);
t185 = qJD(4) * t131;
t187 = qJD(2) * t129;
t225 = t126 * t185 + t187;
t194 = t130 * t125;
t87 = t128 * t124 - t194;
t67 = t87 * t129;
t88 = t130 * t124 + t128 * t125;
t74 = t88 * qJD(1);
t76 = t88 * qJD(6);
t210 = t129 * t74 + t76;
t224 = t210 * t113;
t183 = qJD(6) * t130;
t184 = qJD(6) * t128;
t223 = -t124 * t184 + t125 * t183;
t127 = pkin(1) + qJ(3);
t222 = qJD(1) * t127;
t122 = t129 ^ 2;
t123 = t131 ^ 2;
t221 = qJD(1) * (0.2e1 * t122 - t123);
t191 = t122 + t123;
t220 = t191 * qJD(2);
t176 = qJD(1) * qJD(4);
t160 = t129 * t176;
t151 = t124 * t160;
t152 = t125 * t160;
t208 = -t128 * t152 - t130 * t151;
t16 = qJD(6) * t39 + t208;
t170 = t124 * t180;
t209 = -t128 * t170 + t180 * t194 + t223;
t219 = -t16 * t88 - t209 * t36;
t216 = t125 * pkin(8);
t175 = t129 * t216;
t135 = (pkin(5) * t131 + t175) * qJD(1);
t148 = pkin(4) * t131 + qJ(5) * t129;
t69 = qJD(4) * t148 - t131 * qJD(5) + qJD(3);
t52 = t69 * qJD(1);
t177 = qJD(1) * qJD(2);
t116 = qJD(1) * qJ(2) + qJD(3);
t107 = -pkin(7) * qJD(1) + t116;
t193 = t131 * t107;
t55 = t129 * t177 + (qJD(5) + t193) * qJD(4);
t18 = -t124 * t55 + t125 * t52;
t12 = qJD(4) * t135 + t18;
t19 = t124 * t52 + t125 * t55;
t13 = pkin(8) * t151 + t19;
t95 = pkin(4) * t129 - qJ(5) * t131 + t127;
t62 = qJD(1) * t95 - qJD(2);
t94 = t129 * t107;
t79 = qJD(4) * qJ(5) + t94;
t25 = -t124 * t79 + t125 * t62;
t14 = pkin(5) * t180 - pkin(8) * t85 + t25;
t26 = t124 * t62 + t125 * t79;
t17 = pkin(8) * t84 + t26;
t141 = t128 * t17 - t130 * t14;
t1 = -qJD(6) * t141 + t128 * t12 + t130 * t13;
t90 = t148 * qJD(1);
t42 = -t124 * t193 + t125 * t90;
t24 = t135 + t42;
t43 = t124 * t90 + t125 * t193;
t33 = pkin(8) * t170 + t43;
t213 = pkin(8) + qJ(5);
t102 = t213 * t124;
t103 = t213 * t125;
t47 = -t102 * t130 - t103 * t128;
t218 = -qJD(5) * t87 + qJD(6) * t47 - t128 * t24 - t130 * t33;
t48 = -t102 * t128 + t103 * t130;
t217 = -qJD(5) * t88 - qJD(6) * t48 + t128 * t33 - t130 * t24;
t214 = t39 * t36;
t65 = t88 * t129;
t68 = t87 * t131;
t212 = -qJD(4) * t68 - qJD(6) * t65 - t74;
t211 = qJD(1) * t87 + qJD(6) * t67 - t185 * t88;
t207 = t124 * t84;
t206 = t124 * t85;
t205 = t125 * t84;
t204 = t125 * t85;
t186 = qJD(4) * t129;
t89 = t107 * t186;
t63 = -t131 * t177 + t89;
t203 = t63 * t124;
t202 = t63 * t125;
t201 = t63 * t131;
t195 = t126 * t129;
t51 = t124 * t95 + t125 * t195;
t200 = qJD(4) * t87;
t199 = qJD(4) * t88;
t133 = qJD(1) ^ 2;
t198 = t122 * t133;
t197 = t124 * t133;
t196 = t125 * t133;
t159 = t131 * t176;
t192 = t122 * t196 + t124 * t159;
t132 = qJD(4) ^ 2;
t189 = -t132 - t133;
t179 = t131 * qJD(2);
t108 = -qJD(2) + t222;
t178 = qJD(2) - t108;
t32 = t124 * t69 + t125 * t225;
t118 = 0.2e1 * t177;
t174 = 0.2e1 * qJD(3) * qJD(1);
t173 = t122 * t197;
t172 = t85 * t188;
t171 = t131 * t133 * t129;
t169 = t129 * t182;
t167 = t128 * t186;
t166 = t130 * t186;
t161 = t209 * t113;
t158 = -t124 * t126 + pkin(5);
t157 = pkin(5) * t124 - t126;
t156 = qJD(4) * pkin(4) - qJD(5);
t155 = -t84 + t181;
t154 = t108 + t222;
t153 = t178 * qJD(1);
t150 = t129 * t159;
t71 = -t156 - t193;
t149 = t126 * t188 - t71;
t146 = -t19 * t124 - t18 * t125;
t145 = -t18 * t124 + t19 * t125;
t144 = -t124 * t26 - t125 * t25;
t143 = t124 * t25 - t125 * t26;
t142 = t205 + t206;
t6 = t128 * t14 + t130 * t17;
t81 = t125 * t95;
t34 = t129 * t158 - t131 * t216 + t81;
t40 = -pkin(8) * t124 * t131 + t51;
t10 = -t128 * t40 + t130 * t34;
t11 = t128 * t34 + t130 * t40;
t140 = qJD(2) + t154;
t139 = qJD(1) * t155;
t138 = -t126 * t132 + t174;
t15 = -t128 * t151 + t130 * t152 - t183 * t84 + t184 * t85;
t137 = t15 * t87 - t210 * t39;
t136 = (t85 - t182) * t180;
t44 = t89 + (-pkin(5) * t169 - t179) * qJD(1);
t134 = -qJ(5) * t185 + (t156 + t71) * t129;
t2 = -qJD(6) * t6 + t12 * t130 - t128 * t13;
t120 = t125 ^ 2;
t119 = t124 ^ 2;
t117 = t123 * t133;
t115 = -pkin(5) * t125 - pkin(4);
t104 = 0.2e1 * t150;
t83 = t157 * t131;
t66 = t88 * t131;
t61 = -pkin(5) * t170 + t94;
t58 = -t157 * t186 - t179;
t57 = t125 * t69;
t50 = -t124 * t195 + t81;
t41 = -pkin(5) * t84 + t71;
t31 = -t124 * t225 + t57;
t30 = -t124 * t166 - t125 * t167 + t131 * t223;
t28 = -t124 * t167 + t125 * t166 + t131 * t76;
t21 = pkin(8) * t169 + t32;
t20 = -t124 * t187 + t57 + (t131 * t158 + t175) * qJD(4);
t4 = -qJD(6) * t11 - t128 * t21 + t130 * t20;
t3 = qJD(6) * t10 + t128 * t20 + t130 * t21;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, qJ(2) * t118, 0, 0, 0, 0, 0, 0, 0, t118, t174, t116 * qJD(2) + t108 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t127) * qJD(1), -0.2e1 * t150, 0.2e1 * (t122 - t123) * t176, -t132 * t129, t104, -t132 * t131, 0, t129 * t138 + t140 * t185, t131 * t138 - t140 * t186, -t191 * t118, t154 * qJD(3) + (qJD(1) * t126 + t107) * t220 (-t120 * t188 - t204) * t186 (-t205 + (t85 + 0.2e1 * t162) * t124) * t186 (-t125 * t221 + t131 * t85) * qJD(4) (-t119 * t188 + t207) * t186 (t124 * t221 + t131 * t84) * qJD(4), t104 (qJD(2) * t84 + t203) * t131 + (qJD(1) * t31 + t18) * t129 + ((qJD(1) * t50 + t25) * t131 + (t124 * t149 - t126 * t84) * t129) * qJD(4) (-qJD(2) * t85 + t202) * t131 + (-qJD(1) * t32 - t19) * t129 + ((-qJD(1) * t51 - t26) * t131 + (t125 * t149 + t126 * t85) * t129) * qJD(4), -t31 * t85 + t32 * t84 + t146 * t131 + ((t124 * t51 + t125 * t50) * qJD(1) - t144) * t186, -t71 * t179 + t18 * t50 + t19 * t51 + t25 * t31 + t26 * t32 + (t186 * t71 - t201) * t126, t15 * t68 - t28 * t39, t15 * t66 + t16 * t68 + t28 * t36 - t30 * t39, -t28 * t113 - t15 * t129 + (-qJD(1) * t68 + t39) * t185, t16 * t66 + t30 * t36, -t30 * t113 - t16 * t129 + (-qJD(1) * t66 - t36) * t185 (t113 + t180) * t185, t4 * t113 + t2 * t129 + t83 * t16 + t41 * t30 + t58 * t36 + t44 * t66 + (qJD(1) * t10 - t141) * t185, -t1 * t129 - t3 * t113 - t83 * t15 - t41 * t28 + t58 * t39 - t44 * t68 + (-qJD(1) * t11 - t6) * t185, -t1 * t66 + t10 * t15 - t11 * t16 - t141 * t28 + t2 * t68 - t3 * t36 - t30 * t6 - t39 * t4, t1 * t11 + t10 * t2 - t141 * t4 + t3 * t6 + t41 * t58 + t44 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t133 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t133, 0 (-qJD(3) - t116) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t159, 0.2e1 * t160, t117 + t198 (-t107 * t191 - qJD(3)) * qJD(1), 0, 0, 0, 0, 0, 0, t173 + (-t84 - t181) * t188, t172 + t192 ((-t119 - t120) * qJD(4) - t142) * t180 (t129 * t143 + t131 * t71) * qJD(1) + t146, 0, 0, 0, 0, 0, 0, t224 + (t36 + t200) * t188, t161 + (t39 + t199) * t188, t137 - t219, -t1 * t88 - t141 * t210 + t188 * t41 + t2 * t87 - t209 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t153, 0, 0, 0, 0, 0, 0, t189 * t129, t189 * t131, 0 (-t108 + t220) * qJD(1), 0, 0, 0, 0, 0, 0 (-t196 + (-t84 - t163) * qJD(4)) * t129 (t197 + (t85 - t162) * qJD(4)) * t129, t142 * t185 + (t204 - t207) * qJD(1), -t201 + t145 * t129 + t144 * qJD(1) + (t129 * t71 - t131 * t143) * qJD(4), 0, 0, 0, 0, 0, 0, -t131 * t16 + t211 * t113 + (t129 * t36 - t188 * t65) * qJD(4), t131 * t15 - t212 * t113 + (t129 * t39 + t188 * t67) * qJD(4), -t15 * t65 + t16 * t67 - t211 * t39 - t212 * t36, -t1 * t67 - t44 * t131 - t141 * t211 + t186 * t41 - t2 * t65 + t212 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t117 - t198, 0, -t171, 0, 0, t131 * t153, -t178 * t180, 0, 0, t125 * t136 (-t206 + t205 + (t119 - t120) * qJD(4)) * t180, -t172 + t192, t129 * t124 * t139, t131 * t139 - t173, -t171, t84 * t94 - t202 + (t124 * t134 - t129 * t42 - t131 * t25) * qJD(1), -t85 * t94 + t203 + (t125 * t134 + t129 * t43 + t131 * t26) * qJD(1), t42 * t85 - t43 * t84 + (qJD(5) * t84 - t180 * t25 + t19) * t125 + (qJD(5) * t85 - t180 * t26 - t18) * t124, -t63 * pkin(4) + qJ(5) * t145 - qJD(5) * t143 - t25 * t42 - t26 * t43 - t71 * t94, -t15 * t88 + t209 * t39, t137 + t219, t161 + (-t39 + t199) * t188, t16 * t87 + t210 * t36, -t224 + (t36 - t200) * t188, -t113 * t188, t115 * t16 - t61 * t36 + t44 * t87 + t210 * t41 + t217 * t113 + (qJD(4) * t47 + t141) * t188, -t115 * t15 - t61 * t39 + t44 * t88 + t209 * t41 - t218 * t113 + (-qJD(4) * t48 + t6) * t188, -t1 * t87 + t141 * t209 + t15 * t47 - t16 * t48 - t2 * t88 - t210 * t6 - t217 * t39 - t218 * t36, t1 * t48 + t44 * t115 - t141 * t217 + t2 * t47 + t218 * t6 - t41 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t155 * t180, -t84 ^ 2 - t85 ^ 2, t25 * t85 - t26 * t84 + t63, 0, 0, 0, 0, 0, 0, t39 * t113 + t16, -t15 - t226, -t227 - t228, -t141 * t39 + t36 * t6 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t227 - t228, -t15 + t226, -t214, -t208 + (-qJD(6) + t113) * t39, t159, t6 * t113 - t41 * t39 + t2, -t113 * t141 + t41 * t36 - t1, 0, 0;];
tauc_reg  = t5;
