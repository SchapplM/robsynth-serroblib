% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:58
% EndTime: 2019-03-08 21:44:04
% DurationCPUTime: 1.98s
% Computational Cost: add. (1902->288), mult. (4561->408), div. (0->0), fcn. (2994->8), ass. (0->169)
t214 = pkin(4) + pkin(8);
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t110 = cos(pkin(6));
t183 = qJD(1) * t110;
t113 = sin(qJ(2));
t109 = sin(pkin(6));
t184 = qJD(1) * t109;
t160 = t113 * t184;
t80 = qJD(2) * pkin(8) + t160;
t48 = t112 * t80 - t115 * t183;
t228 = qJD(4) + t48;
t116 = cos(qJ(2));
t182 = qJD(2) * t109;
t152 = qJD(1) * t182;
t140 = t116 * t152;
t178 = qJD(3) * t110;
t227 = qJD(1) * t178 + t140;
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t117 = -pkin(3) - pkin(9);
t181 = qJD(2) * t112;
t221 = pkin(4) * t181 + t228;
t28 = t117 * qJD(3) + t221;
t155 = t116 * t184;
t148 = -t112 * qJ(4) - pkin(2);
t66 = t117 * t115 + t148;
t39 = qJD(2) * t66 - t155;
t12 = t111 * t28 + t114 * t39;
t168 = qJD(2) * qJD(3);
t149 = t115 * t168;
t174 = qJD(3) * t115;
t26 = t112 * t227 + t80 * t174;
t19 = pkin(4) * t149 + t26;
t138 = pkin(9) * t112 - qJ(4) * t115;
t170 = t112 * qJD(4);
t122 = t138 * qJD(3) - t170;
t150 = t112 * t168;
t204 = pkin(3) * t150 + t113 * t152;
t27 = qJD(2) * t122 + t204;
t147 = -t111 * t27 + t114 * t19;
t121 = -t12 * qJD(5) + t147;
t177 = qJD(3) * t111;
t179 = qJD(2) * t115;
t69 = t114 * t179 + t177;
t44 = t69 * qJD(5) - t111 * t150;
t154 = t111 * t179;
t175 = qJD(3) * t114;
t71 = -t154 + t175;
t1 = pkin(5) * t149 + t44 * qJ(6) - t71 * qJD(6) + t121;
t172 = qJD(5) * t114;
t167 = -t111 * t19 - t114 * t27 - t28 * t172;
t173 = qJD(5) * t111;
t128 = -t39 * t173 - t167;
t45 = qJD(3) * t172 - qJD(5) * t154 - t114 * t150;
t2 = -qJ(6) * t45 - qJD(6) * t69 + t128;
t11 = -t111 * t39 + t114 * t28;
t6 = -qJ(6) * t71 + t11;
t99 = qJD(5) + t181;
t5 = pkin(5) * t99 + t6;
t7 = -qJ(6) * t69 + t12;
t226 = -(t5 * t99 - t2) * t111 + (t7 * t99 + t1) * t114;
t190 = t114 * t116;
t77 = t214 * t174;
t225 = -(-t111 * t113 + t112 * t190) * t184 + t114 * t77;
t106 = qJD(3) * qJ(4);
t49 = t112 * t183 + t115 * t80;
t42 = -t106 - t49;
t191 = t111 * t116;
t176 = qJD(3) * t112;
t101 = pkin(3) * t176;
t51 = t101 + t122;
t86 = t214 * t112;
t224 = -t111 * t77 - t114 * t51 - t86 * t172 + t66 * t173 + (t112 * t191 + t113 * t114) * t184;
t41 = -qJD(3) * pkin(3) + t228;
t210 = t69 * t99;
t223 = t44 - t210;
t209 = t71 * t99;
t222 = -t45 + t209;
t118 = qJD(3) ^ 2;
t125 = -qJ(4) * t174 - t170;
t32 = qJD(2) * t125 + t204;
t57 = t101 + t125;
t218 = qJD(2) * (-t57 + t160) - pkin(8) * t118 - t32;
t158 = t116 * t182;
t119 = qJD(2) ^ 2;
t193 = t109 * t119;
t162 = t113 * t193;
t194 = t109 * t113;
t163 = t112 * t194;
t33 = -qJD(3) * t163 + (t158 + t178) * t115;
t217 = qJD(3) * (t115 * t158 + t33) - t112 * t162;
t141 = t112 * t158;
t61 = t110 * t112 + t115 * t194;
t34 = qJD(3) * t61 + t141;
t216 = qJD(3) * (t34 + t141) + t115 * t162;
t215 = t71 ^ 2;
t213 = t5 - t6;
t196 = qJ(6) * t115;
t145 = -t66 + t196;
t212 = pkin(5) * t174 + t145 * t172 + (-qJ(6) * t176 - qJD(5) * t86 + qJD(6) * t115 - t51) * t111 + t225;
t171 = qJD(5) * t115;
t157 = t111 * t171;
t169 = t114 * qJD(6);
t211 = -t115 * t169 + (t112 * t175 + t157) * qJ(6) - t224;
t38 = pkin(4) * t179 + t49;
t102 = pkin(3) * t181;
t55 = qJD(2) * t138 + t102;
t208 = t111 * t38 + t114 * t55;
t146 = -t111 * t55 + t114 * t38;
t187 = qJ(6) - t117;
t192 = t111 * t112;
t207 = t187 * t173 - t169 - (pkin(5) * t115 - qJ(6) * t192) * qJD(2) - t146;
t83 = t187 * t114;
t206 = -qJ(6) * t114 * t181 - qJD(5) * t83 - t111 * qJD(6) - t208;
t205 = t111 * t86 + t114 * t66;
t203 = qJD(2) * pkin(2);
t202 = t114 * t44;
t201 = t114 * t99;
t200 = t115 * t71;
t199 = t117 * t99;
t105 = qJD(3) * qJD(4);
t166 = t115 * t227 - t80 * t176;
t21 = -t105 - t166;
t15 = -pkin(4) * t150 - t21;
t198 = t15 * t111;
t197 = t15 * t114;
t84 = -pkin(3) * t115 + t148;
t195 = qJD(2) * t84;
t189 = t118 * t112;
t188 = t118 * t115;
t87 = t214 * t115;
t107 = t112 ^ 2;
t108 = t115 ^ 2;
t185 = t107 - t108;
t180 = qJD(2) * t113;
t165 = t112 * t201;
t164 = t99 * t172;
t161 = t112 * t119 * t115;
t159 = t109 * t180;
t156 = t114 * t171;
t29 = t106 + t38;
t137 = -qJD(2) * t108 + t112 * t99;
t135 = t48 * qJD(3) + t166;
t134 = qJD(3) * t49 - t26;
t133 = t111 * t99;
t60 = -t110 * t115 + t163;
t130 = t109 * t190 - t111 * t60;
t35 = t109 * t191 + t114 * t60;
t129 = t112 * t29 + t117 * t174;
t81 = -t155 - t203;
t124 = qJD(3) * (t155 + t81 - t203);
t50 = -t155 + t195;
t123 = qJD(3) * (-t155 - t50 - t195);
t10 = pkin(5) * t45 + t15;
t120 = t26 * t112 - t21 * t115 + (t112 * t42 + t115 * t41) * qJD(3);
t92 = t114 * t149;
t82 = t187 * t111;
t76 = t214 * t176;
t75 = -qJ(4) * t179 + t102;
t74 = t114 * t86;
t68 = t69 ^ 2;
t40 = t50 * t181;
t25 = -t114 * t196 + t205;
t20 = t112 * pkin(5) + t111 * t145 + t74;
t18 = pkin(5) * t69 + qJD(6) + t29;
t9 = t35 * qJD(5) + t34 * t111 + t114 * t159;
t8 = t130 * qJD(5) - t111 * t159 + t34 * t114;
t3 = [0, 0, -t162, -t116 * t193, 0, 0, 0, 0, 0, -t216, -t217 (t112 * t34 + t115 * t33 + (-t112 * t61 + t115 * t60) * qJD(3)) * qJD(2), t216, t217, -t21 * t61 + t26 * t60 - t42 * t33 + t41 * t34 + (-t116 * t32 + t180 * t50) * t109, 0, 0, 0, 0, 0, t149 * t35 + t33 * t69 + t45 * t61 + t8 * t99, t130 * t149 + t33 * t71 - t44 * t61 - t9 * t99, t130 * t45 + t35 * t44 - t69 * t9 - t71 * t8, t1 * t35 + t10 * t61 - t130 * t2 + t18 * t33 + t5 * t8 + t7 * t9; 0, 0, 0, 0, 0.2e1 * t112 * t149, -0.2e1 * t185 * t168, t188, -t189, 0, -pkin(8) * t188 + t112 * t124, pkin(8) * t189 + t115 * t124 (-t107 - t108) * t140 + t120, t112 * t123 - t218 * t115, t218 * t112 + t115 * t123, t32 * t84 + t50 * t57 + (-t113 * t50 + (-t112 * t41 + t115 * t42) * t116) * t184 + t120 * pkin(8), -t71 * t156 + (t115 * t44 + t176 * t71) * t111 (-t111 * t69 + t114 * t71) * t176 + (t111 * t45 + t202 + (t111 * t71 + t114 * t69) * qJD(5)) * t115, -t99 * t156 - t44 * t112 + (t111 * t137 + t200) * qJD(3), t99 * t157 - t45 * t112 + (t114 * t137 - t115 * t69) * qJD(3) (t99 + t181) * t174, t87 * t45 - t76 * t69 + (-t111 * t51 + t225) * t99 + (-t175 * t29 + t147) * t112 + (-t112 * t12 - t205 * t99) * qJD(5) + (-t69 * t155 - t29 * t173 + t197 + ((-t111 * t66 + t74) * qJD(2) + t11) * qJD(3)) * t115, -t87 * t44 - t76 * t71 + t224 * t99 + ((qJD(3) * t29 + qJD(5) * t39) * t111 + t167) * t112 + (-t71 * t155 - t29 * t172 - t198 + (-qJD(2) * t205 - t12) * qJD(3)) * t115, t20 * t44 - t25 * t45 - t212 * t71 - t211 * t69 + (-t111 * t5 + t114 * t7) * t176 + (t1 * t111 - t114 * t2 + (t111 * t7 + t114 * t5) * qJD(5)) * t115, t2 * t25 + t1 * t20 + t10 * t87 + t211 * t7 + t212 * t5 + (-t18 * t155 + (t10 * t114 - t173 * t18) * pkin(5)) * t115 + t18 * (-pkin(5) * t114 - t214) * t176; 0, 0, 0, 0, -t161, t185 * t119, 0, 0, 0, -t181 * t81 + t134, -t179 * t81 - t135, 0, -t179 * t75 - t134 + t40, 0.2e1 * t105 + (t112 * t75 + t115 * t50) * qJD(2) + t135, -t26 * pkin(3) - t21 * qJ(4) - t228 * t42 - t41 * t49 - t50 * t75, -t133 * t71 - t202 (-t45 - t209) * t114 + (t44 + t210) * t111, -t99 * t173 + t92 + (-t192 * t99 - t200) * qJD(2), -t164 + (-t165 + (t69 - t177) * t115) * qJD(2), -t99 * t179, qJ(4) * t45 + t198 - t146 * t99 + t221 * t69 + (-t111 * t199 + t114 * t29) * qJD(5) + (-t11 * t115 + t114 * t129) * qJD(2), -qJ(4) * t44 + t197 + t208 * t99 + t221 * t71 + (-t111 * t29 - t114 * t199) * qJD(5) + (-t111 * t129 + t12 * t115) * qJD(2), -t206 * t69 - t207 * t71 - t44 * t83 + t45 * t82 - t226, -t2 * t82 - t1 * t83 + t10 * (pkin(5) * t111 + qJ(4)) + t206 * t7 + t207 * t5 + (pkin(5) * t201 + t221) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t107 * t119 - t118, qJD(3) * t42 + t26 + t40, 0, 0, 0, 0, 0, -qJD(3) * t69 - t133 * t99 + t92, -t164 - qJD(3) * t71 + (-t111 * t174 - t165) * qJD(2), t222 * t111 + t223 * t114, -t18 * qJD(3) + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t69, -t68 + t215, -t223, t222, t149, t12 * t99 - t29 * t71 + t121, t11 * t99 + t29 * t69 - t128, pkin(5) * t44 - t213 * t69, t213 * t7 + (-t18 * t71 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 - t215, t5 * t71 + t69 * t7 + t10;];
tauc_reg  = t3;
