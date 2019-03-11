% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:01
% EndTime: 2019-03-09 03:13:06
% DurationCPUTime: 1.98s
% Computational Cost: add. (2250->300), mult. (4804->387), div. (0->0), fcn. (2662->6), ass. (0->169)
t92 = sin(pkin(9)) * pkin(1) + pkin(7);
t207 = pkin(4) + t92;
t110 = sin(qJ(3));
t112 = cos(qJ(3));
t82 = t92 * qJD(1);
t56 = -t112 * qJD(2) + t110 * t82;
t217 = qJD(4) + t56;
t169 = qJD(1) * qJD(3);
t152 = t112 * t169;
t149 = pkin(5) * t152;
t109 = sin(qJ(5));
t111 = cos(qJ(5));
t171 = qJD(5) * t111;
t172 = qJD(5) * t109;
t113 = -pkin(3) - pkin(8);
t178 = qJD(1) * t110;
t181 = pkin(4) * t178 + t217;
t26 = t113 * qJD(3) + t181;
t138 = pkin(8) * t110 - qJ(4) * t112;
t170 = t110 * qJD(4);
t118 = t138 * qJD(3) - t170;
t154 = t110 * t169;
t90 = pkin(3) * t154;
t29 = t118 * qJD(1) + t90;
t168 = qJD(2) * qJD(3);
t173 = qJD(3) * t112;
t49 = t110 * t168 + t82 * t173;
t30 = pkin(4) * t152 + t49;
t93 = -cos(pkin(9)) * pkin(1) - pkin(2);
t127 = -t110 * qJ(4) + t93;
t53 = t113 * t112 + t127;
t35 = t53 * qJD(1);
t151 = -t109 * t29 + t111 * t30 - t35 * t171 - t26 * t172;
t2 = -t149 - t151;
t9 = t109 * t26 + t111 * t35;
t91 = qJD(5) + t178;
t7 = qJ(6) * t91 + t9;
t216 = t7 * t91 - t2;
t57 = t110 * qJD(2) + t112 * t82;
t46 = -qJD(3) * qJ(4) - t57;
t43 = -qJD(3) * pkin(3) + t217;
t176 = qJD(3) * t109;
t177 = qJD(1) * t112;
t75 = t111 * t177 + t176;
t203 = t75 * t91;
t39 = t75 * qJD(5) - t109 * t154;
t215 = t39 - t203;
t155 = t109 * t177;
t174 = qJD(3) * t111;
t77 = -t155 + t174;
t201 = t77 * t91;
t153 = t111 * t169;
t40 = qJD(3) * t171 - qJD(5) * t155 - t110 * t153;
t214 = -t40 + t201;
t162 = t91 * t171;
t183 = t110 * t111;
t165 = t91 * t183;
t212 = qJD(1) * (t109 * t173 + t165) + qJD(3) * t77 + t162;
t66 = t207 * t110;
t135 = t109 * t66 + t111 * t53;
t175 = qJD(3) * t110;
t97 = pkin(3) * t175;
t44 = t97 + t118;
t63 = t207 * t173;
t211 = -t135 * qJD(5) - t109 * t44 + t111 * t63;
t167 = -t109 * t30 - t111 * t29 - t26 * t171;
t122 = t35 * t172 + t167;
t146 = qJ(6) * t152;
t1 = qJD(6) * t91 - t122 + t146;
t8 = -t109 * t35 + t111 * t26;
t185 = qJD(6) - t8;
t6 = -pkin(5) * t91 + t185;
t141 = t109 * t6 + t111 * t7;
t210 = qJD(5) * t141 + t1 * t109 - t2 * t111;
t209 = t77 ^ 2;
t96 = pkin(4) * t177;
t33 = t96 - t46;
t11 = pkin(5) * t75 - qJ(6) * t77 + t33;
t206 = t11 * t77;
t103 = qJD(3) * qJD(4);
t195 = t112 * t168 - t82 * t175;
t38 = -t103 - t195;
t20 = -pkin(4) * t154 - t38;
t4 = pkin(5) * t40 + qJ(6) * t39 - qJD(6) * t77 + t20;
t205 = t4 * t109;
t204 = t4 * t111;
t202 = t77 * t75;
t140 = pkin(5) * t111 + qJ(6) * t109;
t126 = -pkin(4) - t140;
t200 = -t140 * qJD(5) + t111 * qJD(6) + t126 * t178 - t217;
t31 = t110 * t40;
t199 = t75 * t173 + t31;
t42 = t96 + t57;
t98 = pkin(3) * t178;
t58 = t138 * qJD(1) + t98;
t198 = t109 * t42 + t111 * t58;
t32 = t110 * t39;
t197 = t77 * t173 - t32;
t158 = t110 * t174;
t163 = t91 * t172;
t72 = t112 * t163;
t196 = t91 * t158 + t72;
t194 = t109 * t40;
t193 = t109 * t77;
t192 = t110 * t91;
t191 = t111 * t75;
t190 = t112 * t77;
t189 = t113 * t91;
t188 = t20 * t109;
t187 = t20 * t111;
t186 = t39 * t111;
t67 = t207 * t112;
t65 = -pkin(3) * t112 + t127;
t47 = qJD(1) * t65;
t83 = qJD(1) * t93;
t184 = t11 * qJD(3);
t114 = qJD(3) ^ 2;
t182 = t114 * t110;
t101 = t114 * t112;
t105 = t110 ^ 2;
t106 = t112 ^ 2;
t180 = t105 - t106;
t179 = qJD(1) * t106;
t166 = t109 * t189;
t164 = t111 * t189;
t115 = qJD(1) ^ 2;
t161 = t110 * t115 * t112;
t160 = t109 * t179;
t159 = t109 * t175;
t157 = t113 * t173;
t156 = t112 * t171;
t148 = t91 * t156;
t147 = t77 * t158;
t86 = t111 * t152;
t145 = t56 * qJD(3) + t195;
t144 = qJD(3) * t57 - t49;
t142 = t109 * t7 - t111 * t6;
t139 = -pkin(5) * t109 + qJ(6) * t111;
t137 = -t109 * t58 + t111 * t42;
t134 = -t109 * t53 + t111 * t66;
t133 = -0.2e1 * qJD(3) * t47;
t132 = 0.2e1 * qJD(3) * t83;
t131 = t179 - t192;
t129 = t109 * t91;
t128 = t91 * t159 - t148;
t125 = t9 * t91 + t151;
t120 = -qJ(4) * t173 - t170;
t48 = t120 * qJD(1) + t90;
t64 = t120 + t97;
t124 = qJD(1) * t64 + t114 * t92 + t48;
t121 = t109 * t63 + t111 * t44 + t66 * t171 - t53 * t172;
t117 = -qJD(3) * t75 - t91 * t129 + t86;
t116 = t49 * t110 - t38 * t112 + (t110 * t46 + t112 * t43) * qJD(3);
t81 = qJ(4) - t139;
t80 = t113 * t86;
t79 = -qJ(4) * t177 + t98;
t62 = t207 * t175;
t52 = t75 * t159;
t36 = t47 * t178;
t28 = t140 * t112 + t67;
t27 = pkin(5) * t77 + qJ(6) * t75;
t15 = -pkin(5) * t110 - t134;
t14 = qJ(6) * t110 + t135;
t13 = -pkin(5) * t177 - t137;
t12 = qJ(6) * t177 + t198;
t10 = (t139 * qJD(5) + qJD(6) * t109) * t112 + (t126 - t92) * t175;
t5 = -pkin(5) * t173 - t211;
t3 = qJ(6) * t173 + qJD(6) * t110 + t121;
t16 = [0, 0, 0, 0, 0.2e1 * t110 * t152, -0.2e1 * t180 * t169, t101, -t182, 0, -t101 * t92 + t110 * t132, t112 * t132 + t182 * t92, t116, t110 * t133 + t112 * t124, -t110 * t124 + t112 * t133, t116 * t92 + t47 * t64 + t48 * t65, -t77 * t156 + (t112 * t39 + t175 * t77) * t109, t147 - t52 + (t194 + t186 + (t191 + t193) * qJD(5)) * t112, -qJD(3) * t160 + t128 + t197, -t31 + (-t111 * t179 - t112 * t75) * qJD(3) + t196 (t91 + t178) * t173, t211 * t91 - t62 * t75 + t67 * t40 + (-t174 * t33 + t151) * t110 + (-t33 * t172 + t187 + (qJD(1) * t134 + t8) * qJD(3)) * t112, -t121 * t91 - t62 * t77 - t67 * t39 + ((qJD(3) * t33 + qJD(5) * t35) * t109 + t167) * t110 + (-t33 * t171 - t188 + (-t135 * qJD(1) - t9) * qJD(3)) * t112, t10 * t75 + t28 * t40 - t5 * t91 + (-t11 * t174 - t2) * t110 + (-t11 * t172 + t204 + (-qJD(1) * t15 - t6) * qJD(3)) * t112, -t14 * t40 - t15 * t39 - t3 * t75 + t5 * t77 + t141 * t175 + (qJD(5) * t142 - t1 * t111 - t109 * t2) * t112, -t10 * t77 + t28 * t39 + t3 * t91 + (-t11 * t176 + t1) * t110 + (t11 * t171 + t205 + (qJD(1) * t14 + t7) * qJD(3)) * t112, t1 * t14 + t10 * t11 + t15 * t2 + t28 * t4 + t3 * t7 + t5 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, -t101, 0, t182, t101, -t38 * t110 - t49 * t112 + (t110 * t43 - t112 * t46) * qJD(3), 0, 0, 0, 0, 0, -t106 * t153 + t196 + t199, t131 * t176 + t148 + t197, -t131 * t174 + t199 + t72, -t147 - t52 + (t194 - t186 + (t191 - t193) * qJD(5)) * t112, t32 + (-t160 - t190) * qJD(3) + t128 (qJD(3) * t142 + t4) * t110 + (t184 - t210) * t112; 0, 0, 0, 0, -t161, t180 * t115, 0, 0, 0, -t178 * t83 + t144, -t177 * t83 - t145, 0, -t177 * t79 - t144 + t36, 0.2e1 * t103 + (t110 * t79 + t112 * t47) * qJD(1) + t145, -t49 * pkin(3) - t38 * qJ(4) - t217 * t46 - t43 * t57 - t47 * t79, -t129 * t77 - t186 (-t40 - t201) * t111 + (t39 + t203) * t109, -t163 + t86 + (-t109 * t192 - t190) * qJD(1), -t162 + (-t165 + (t75 - t176) * t112) * qJD(1), -t91 * t177, t80 + qJ(4) * t40 + t188 - t137 * t91 + t181 * t75 + (t111 * t33 - t166) * qJD(5) + (-t8 * t112 + t183 * t33) * qJD(1), -qJ(4) * t39 + t187 + t198 * t91 + t181 * t77 + (-t109 * t33 - t164) * qJD(5) + (t9 * t112 + (-t110 * t33 - t157) * t109) * qJD(1), t205 + t13 * t91 + t81 * t40 + t80 - t200 * t75 + (t11 * t111 - t166) * qJD(5) + (t11 * t183 + t112 * t6) * qJD(1), t12 * t75 - t13 * t77 + (-t7 * t178 + t113 * t39 + t2 + (-t113 * t75 - t7) * qJD(5)) * t111 + (-t6 * t178 - t113 * t40 - t1 + (t113 * t77 - t6) * qJD(5)) * t109, -t204 - t12 * t91 + t81 * t39 + t200 * t77 + (t109 * t11 + t164) * qJD(5) + (-t112 * t7 + (t11 * t110 + t157) * t109) * qJD(1), -t200 * t11 + t113 * t210 - t7 * t12 - t6 * t13 + t4 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t105 * t115 - t114, qJD(3) * t46 + t36 + t49, 0, 0, 0, 0, 0, t117, -t212, t117, t214 * t109 + t215 * t111, t212, -t184 + t216 * t111 + (t6 * t91 + t1) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t75 ^ 2 + t209, -t215, t214, t152, -t33 * t77 + t125, t33 * t75 + t8 * t91 + t122, -t27 * t75 + t125 + 0.2e1 * t149 - t206, pkin(5) * t39 - t40 * qJ(6) + (t7 - t9) * t77 + (t6 - t185) * t75, 0.2e1 * t146 - t11 * t75 + t27 * t77 + (0.2e1 * qJD(6) - t8) * t91 - t122, -t2 * pkin(5) + t1 * qJ(6) - t11 * t27 + t185 * t7 - t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 + t202, -t215, -t91 ^ 2 - t209, t206 - t216;];
tauc_reg  = t16;
