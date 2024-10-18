% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:15
% EndTime: 2024-09-27 17:32:17
% DurationCPUTime: 1.22s
% Computational Cost: add. (2883->222), mult. (5399->315), div. (0->0), fcn. (3466->10), ass. (0->166)
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t105 = qJD(1) + qJD(2);
t101 = qJD(3) + t105;
t108 = sin(pkin(5));
t115 = cos(qJ(4));
t194 = t108 * t115;
t169 = t101 * t194;
t111 = sin(qJ(4));
t198 = t101 * t108;
t170 = t111 * t198;
t217 = -t110 * t170 + t114 * t169;
t112 = sin(qJ(3));
t187 = qJD(3) * t112;
t177 = pkin(2) * t187;
t117 = cos(qJ(2));
t113 = sin(qJ(2));
t116 = cos(qJ(3));
t190 = t113 * t116;
t131 = -t112 * t117 - t190;
t202 = pkin(1) * qJD(1);
t77 = t131 * t202;
t144 = t77 + t177;
t212 = qJD(4) + qJD(5);
t109 = cos(pkin(5));
t184 = qJD(4) * t115;
t165 = t109 * t184;
t193 = t109 * t111;
t179 = t113 * t202;
t151 = t112 * t179;
t178 = t117 * t202;
t89 = t105 * pkin(2) + t178;
t68 = t116 * t89 - t151;
t69 = -t112 * t89 - t116 * t179;
t216 = -pkin(3) * t165 + t115 * t68 + t193 * t69;
t192 = t109 * t115;
t191 = t112 * t113;
t130 = t116 * t117 - t191;
t78 = t130 * t202;
t215 = -t111 * t78 + t192 * t77 - (-t111 * t116 - t112 * t192) * qJD(3) * pkin(2);
t186 = qJD(3) * t116;
t98 = pkin(2) * t116 + pkin(3);
t214 = -t98 * t165 + t193 * t77 + (-pkin(2) * t186 + t78) * t115;
t185 = qJD(4) * t111;
t176 = pkin(4) * t185;
t213 = t108 * (t69 + t176);
t72 = (t110 * t111 - t114 * t115) * t108;
t133 = t110 * t115 + t111 * t114;
t73 = t133 * t108;
t56 = pkin(9) * t198 - t69;
t145 = pkin(10) * t198 + t56;
t210 = pkin(3) * t101;
t58 = t68 + t210;
t175 = t58 * t193;
t22 = t115 * t145 + t175;
t149 = qJD(2) * t178;
t203 = (qJD(2) + qJD(3)) * t151;
t40 = t116 * (qJD(3) * t89 + t149) - t203;
t119 = (qJD(2) * t131 - t113 * t186) * pkin(1);
t118 = qJD(1) * t119;
t41 = -t187 * t89 + t118;
t36 = t41 * t192;
t157 = -t111 * t40 + t36;
t104 = t108 ^ 2;
t195 = t104 * t115;
t200 = t115 * t56;
t211 = t41 * t195 + ((-t175 - t200) * qJD(4) + t157) * t109;
t209 = pkin(4) * t115;
t208 = pkin(10) * t108;
t39 = (-t101 * t209 - t58) * t108;
t63 = -t110 * t169 - t114 * t170;
t207 = t39 * t63;
t206 = t63 * t217;
t182 = t115 * t40 + t58 * t165 + t41 * t193;
t126 = t185 * t56 - t182;
t205 = t126 * t109;
t201 = t114 * t22;
t199 = t101 * t104;
t197 = t101 * t109;
t196 = t104 * t111;
t90 = qJD(4) + t197;
t189 = qJD(4) - t90;
t188 = t111 ^ 2 - t115 ^ 2;
t129 = t145 * t111;
t51 = t58 * t192;
t21 = t51 - t129;
t16 = t90 * pkin(4) + t21;
t134 = -t110 * t16 - t201;
t5 = -qJD(4) * t129 + t182;
t6 = -t22 * qJD(4) + t157;
t167 = -t110 * t5 + t114 * t6;
t121 = qJD(5) * t134 + t167;
t25 = (t101 * t176 - t41) * t108;
t47 = t212 * t73;
t183 = t121 * t109 + t25 * t72 + t39 * t47;
t99 = pkin(1) * t117 + pkin(2);
t52 = t99 * t186 + (qJD(2) * t130 - t113 * t187) * pkin(1);
t53 = -t187 * t99 + t119;
t75 = -pkin(1) * t191 + t116 * t99 + pkin(3);
t181 = t115 * t52 + t75 * t165 + t53 * t193;
t180 = pkin(3) * t193;
t174 = t75 * t193;
t173 = t98 * t193;
t172 = t108 * (-pkin(9) - pkin(10));
t171 = t101 * t195;
t88 = qJD(5) + t90;
t168 = -pkin(4) * t88 - t16;
t166 = t108 * t185;
t164 = -t58 - t210;
t102 = t108 * pkin(9);
t71 = pkin(1) * t190 + t112 * t99 + t102;
t163 = -t71 - t208;
t93 = pkin(2) * t112 + t102;
t162 = -t93 - t208;
t161 = -t53 * t101 - t41;
t160 = t69 * t101 - t41;
t159 = -t101 * t75 - t58;
t94 = pkin(4) * t166;
t158 = t144 * t108 + t94;
t156 = -t111 * t52 + t53 * t192;
t155 = t90 + t197;
t154 = pkin(4) * t170;
t153 = 0.2e1 * qJD(4) * t199;
t152 = t109 * t177;
t150 = t111 * t172;
t17 = qJD(5) * t110 * t22;
t46 = t212 * t72;
t146 = -(t110 * t6 - t17 + (qJD(5) * t16 + t5) * t114) * t109 - t39 * t46 + t25 * t73;
t125 = -pkin(9) * t194 - t180;
t29 = -t111 * t68 + t192 * t69;
t97 = pkin(10) * t194;
t142 = qJD(5) * (-t125 + t97) + t29 - (t115 * t172 - t180) * qJD(4);
t103 = t109 * pkin(4);
t141 = -qJD(5) * (pkin(3) * t192 + t103 + t150) - qJD(4) * t150 + t216;
t140 = -qJD(5) * (t162 * t111 + t98 * t192 + t103) - (t162 * qJD(4) - t152) * t111 + t214;
t127 = -t115 * t93 - t173;
t139 = qJD(5) * (-t127 + t97) - (t115 * t162 - t173) * qJD(4) + t215;
t138 = (-qJD(2) + t105) * t202;
t137 = pkin(1) * qJD(2) * (-qJD(1) - t105);
t136 = t163 * t111;
t135 = (-pkin(2) * t101 - t89) * qJD(3);
t128 = -t115 * t71 - t174;
t123 = t17 - t39 * t217 + (-t22 * t88 - t6) * t110;
t26 = t217 * t212;
t100 = t101 ^ 2;
t87 = (-pkin(3) - t209) * t108;
t79 = t111 * t115 * t153;
t76 = (-t98 - t209) * t108;
t60 = t188 * t153;
t59 = (-t75 - t209) * t108;
t55 = t155 * t108 * t184;
t54 = t155 * t166;
t38 = -t108 * t53 + t94;
t33 = -t128 + t97;
t28 = t192 * t75 + t103 + t136;
t27 = t101 * t47;
t20 = -t217 ^ 2 + t63 ^ 2;
t15 = -t212 * t198 * t133 - t63 * t88;
t14 = -t217 * t88 + t26;
t13 = -t109 * t27 - t47 * t88;
t12 = t109 * t26 - t46 * t88;
t11 = (t115 * t163 - t174) * qJD(4) + t156;
t10 = qJD(4) * t136 + t181;
t9 = t26 * t73 + t46 * t63;
t3 = -t217 * t46 - t26 * t72 - t27 * t73 + t47 * t63;
t1 = [0, 0, 0, 0, t113 * t137, t117 * t137, 0, -t161, -t52 * t101 - t40, t79, -t60, t55, -t54, 0, t156 * t90 + t53 * t171 + (t128 * t90 + t159 * t196) * qJD(4) + t211, -(-t185 * t71 + t181) * t90 + t205 + (t111 * t161 + t159 * t184) * t104, t9, t3, t12, t13, 0, (-t110 * t10 + t114 * t11 + (-t110 * t28 - t114 * t33) * qJD(5)) * t88 - t38 * t217 + t59 * t27 + t183, -(t114 * t10 + t110 * t11 + (-t110 * t33 + t114 * t28) * qJD(5)) * t88 - t38 * t63 + t59 * t26 + t146; 0, 0, 0, 0, t113 * t138, t117 * t138, 0, -t77 * t101 + t112 * t135 + t118, t78 * t101 + (t135 - t149) * t116 + t203, t79, -t60, t55, -t54, 0, (qJD(4) * t127 - t215) * t90 + (-t58 * t185 + (-t115 * t144 - t185 * t98) * t101) * t104 + t211, t205 + ((qJD(4) * t93 + t152) * t111 + t214) * t90 + (-t58 * t184 - t41 * t111 + (t111 * t144 - t184 * t98) * t101) * t104, t9, t3, t12, t13, 0, t76 * t27 + (t110 * t140 - t114 * t139) * t88 - t158 * t217 + t183, t76 * t26 + (t110 * t139 + t114 * t140) * t88 - t158 * t63 + t146; 0, 0, 0, 0, 0, 0, 0, -t160, t68 * t101 - t40, t79, -t60, t55, -t54, 0, -t69 * t171 - t29 * t90 + (t125 * t90 + t164 * t196) * qJD(4) + t211, t205 + (pkin(9) * t166 + t216) * t90 + (t111 * t160 + t164 * t184) * t104, t9, t3, t12, t13, 0, t87 * t27 + (t110 * t141 - t114 * t142) * t88 - t217 * t213 + t183, t87 * t26 + (t110 * t142 + t114 * t141) * t88 - t63 * t213 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t100 * t195, t188 * t104 * t100, t189 * t169, -t189 * t170, 0, t36 - t189 * t200 + (-t40 + (-t109 * t189 + t199) * t58) * t111, (-t111 * t56 + t51) * t90 + t58 * t171 + t126, t206, t20, t14, t15, 0, -(-t110 * t21 - t201) * t88 + t217 * t154 + t207 + (t110 * t168 - t201) * qJD(5) + t167, t63 * t154 + (qJD(5) * t168 + t21 * t88 - t5) * t114 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, t20, t14, t15, 0, -t134 * t88 + t121 + t207, (-t5 + (-qJD(5) + t88) * t16) * t114 + t123;];
tauc_reg = t1;
