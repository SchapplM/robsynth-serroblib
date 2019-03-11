% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:34
% EndTime: 2019-03-08 19:49:41
% DurationCPUTime: 2.31s
% Computational Cost: add. (1827->294), mult. (4329->451), div. (0->0), fcn. (3154->10), ass. (0->164)
t109 = sin(pkin(11));
t111 = cos(pkin(11));
t113 = sin(qJ(6));
t116 = cos(qJ(6));
t214 = -t109 * t113 + t116 * t111;
t216 = t214 * qJD(6);
t114 = sin(qJ(4));
t181 = qJD(2) * t114;
t103 = qJD(6) + t181;
t219 = qJD(6) - t103;
t117 = cos(qJ(4));
t180 = qJD(2) * t117;
t166 = t109 * t180;
t173 = t111 * qJD(4);
t79 = t166 - t173;
t163 = t111 * t180;
t178 = qJD(4) * t109;
t81 = t163 + t178;
t34 = t113 * t81 + t116 * t79;
t134 = t113 * t79 - t116 * t81;
t218 = t134 * t103;
t118 = cos(qJ(2));
t110 = sin(pkin(6));
t185 = qJD(1) * t110;
t115 = sin(qJ(2));
t190 = t114 * t115;
t142 = pkin(4) * t117 + qJ(5) * t114;
t65 = t142 * qJD(4) - t117 * qJD(5) + qJD(3);
t119 = -pkin(2) - pkin(8);
t175 = qJD(4) * t119;
t161 = t117 * t175;
t98 = t111 * t161;
t208 = t109 * t65 - (t109 * t118 + t111 * t190) * t185 + t98;
t217 = t111 * t65 - (-t109 * t190 + t111 * t118) * t185;
t85 = t109 * t116 + t111 * t113;
t127 = t85 * qJD(6);
t112 = cos(pkin(6));
t184 = qJD(1) * t112;
t159 = t118 * t185;
t141 = qJD(3) - t159;
t75 = t119 * qJD(2) + t141;
t215 = -t114 * t184 + t117 * t75;
t172 = qJD(2) * qJD(4);
t158 = t114 * t172;
t147 = t116 * t158;
t148 = t113 * t158;
t14 = -t134 * qJD(6) - t109 * t147 - t111 * t148;
t213 = pkin(9) * t117;
t212 = pkin(9) + qJ(5);
t160 = t115 * t185;
t149 = qJD(2) * t160;
t22 = t114 * t149 + (qJD(5) + t215) * qJD(4);
t37 = (t65 + t159) * qJD(2);
t7 = t109 * t37 + t111 * t22;
t156 = -t109 * t119 + pkin(5);
t191 = t111 * t114;
t171 = pkin(9) * t191;
t211 = (t156 * t117 + t171) * qJD(4) + t217;
t177 = qJD(4) * t114;
t162 = t109 * t177;
t153 = pkin(9) * t162;
t210 = -t153 - t208;
t102 = t117 * t184;
t49 = t114 * t75 + t102;
t43 = qJD(4) * qJ(5) + t49;
t92 = pkin(4) * t114 - qJ(5) * t117 + qJ(3);
t56 = t92 * qJD(2) + t160;
t11 = t109 * t56 + t111 * t43;
t152 = t109 * t161;
t209 = -t152 + t217;
t87 = t142 * qJD(2);
t19 = t109 * t87 + t111 * t215;
t207 = qJD(4) * t102 + t75 * t177;
t129 = t214 * t114;
t206 = qJD(2) * t129 + t216;
t128 = qJD(2) * t85;
t205 = t114 * t128 + t127;
t203 = qJD(2) * pkin(2);
t202 = t103 * t34;
t27 = -t117 * t149 + t207;
t201 = t109 * t27;
t200 = t111 * t27;
t154 = qJD(4) * pkin(4) - qJD(5);
t40 = -t154 - t215;
t198 = t114 * t40;
t183 = qJD(2) * qJ(3);
t88 = t160 + t183;
t196 = t118 * t88;
t189 = t114 * t119;
t55 = t109 * t92 + t111 * t189;
t194 = t110 * t115;
t193 = t110 * t118;
t121 = qJD(2) ^ 2;
t192 = t110 * t121;
t187 = t114 ^ 2 - t117 ^ 2;
t120 = qJD(4) ^ 2;
t186 = -t120 - t121;
t182 = qJD(2) * t110;
t179 = qJD(4) * t103;
t176 = qJD(4) * t117;
t170 = t115 * t192;
t169 = t118 * t192;
t125 = (pkin(5) * t117 + t171) * qJD(2);
t6 = -t109 * t22 + t111 * t37;
t4 = qJD(4) * t125 + t6;
t5 = qJD(2) * t153 + t7;
t168 = -t113 * t5 + t116 * t4;
t167 = t109 * t181;
t165 = t115 * t182;
t164 = t118 * t182;
t157 = t117 * t172;
t155 = pkin(5) * t109 - t119;
t10 = -t109 * t43 + t111 * t56;
t18 = -t109 * t215 + t111 * t87;
t151 = t117 * t165;
t150 = t114 * t165;
t146 = -t88 + t160;
t145 = -t109 * t6 + t111 * t7;
t144 = t113 * t4 + t116 * t5;
t8 = pkin(5) * t181 - pkin(9) * t81 + t10;
t9 = -pkin(9) * t79 + t11;
t143 = t113 * t9 - t116 * t8;
t2 = t113 * t8 + t116 * t9;
t140 = -t10 * t111 - t109 * t11;
t139 = -t10 * t109 + t11 * t111;
t77 = t111 * t92;
t32 = -t111 * t213 + t156 * t114 + t77;
t42 = -t109 * t213 + t55;
t138 = -t113 * t42 + t116 * t32;
t137 = t113 * t32 + t116 * t42;
t72 = t112 * t117 - t114 * t193;
t44 = -t109 * t72 + t111 * t194;
t45 = t109 * t194 + t111 * t72;
t136 = -t113 * t45 + t116 * t44;
t135 = t113 * t44 + t116 * t45;
t71 = t112 * t114 + t117 * t193;
t97 = t212 * t111;
t132 = qJD(5) * t109 + qJD(6) * t97 + t125 + t18;
t96 = t212 * t109;
t131 = pkin(9) * t167 - qJD(5) * t111 + qJD(6) * t96 + t19;
t130 = t146 * qJD(2);
t126 = t214 * qJD(2);
t124 = t146 - t183;
t123 = -qJ(5) * t176 + (t154 + t40) * t114;
t83 = (qJD(3) + t159) * qJD(2);
t122 = t141 * qJD(2) - t119 * t120 + t83;
t13 = -t34 * qJD(6) + t109 * t148 - t111 * t147;
t105 = -pkin(5) * t111 - pkin(4);
t86 = t141 - t203;
t78 = t155 * t117;
t68 = t155 * t177;
t64 = t214 * t117;
t63 = t85 * t117;
t54 = -t109 * t189 + t77;
t47 = t72 * qJD(4) - t151;
t46 = -t71 * qJD(4) + t150;
t31 = -pkin(5) * t167 + t49;
t26 = -t113 * t114 * t173 - t116 * t162 + t117 * t216;
t25 = -qJD(4) * t129 - t117 * t127;
t24 = t109 * t164 + t111 * t46;
t23 = -t109 * t46 + t111 * t164;
t21 = pkin(5) * t79 + t40;
t16 = (-pkin(5) * t162 - t117 * t160) * qJD(2) + t207;
t1 = [0, 0, -t170, -t169, t170, t169 (t115 * t83 + (t196 + (t86 - t159) * t115) * qJD(2)) * t110, 0, 0, 0, 0, 0, t114 * t169 + (-t47 + t151) * qJD(4), t117 * t169 + (-t46 - t150) * qJD(4), t47 * t79 + (t114 * t23 + (-t109 * t114 * t71 + t117 * t44) * qJD(4)) * qJD(2), t47 * t81 + (-t114 * t24 + (-t117 * t45 - t71 * t191) * qJD(4)) * qJD(2), -t23 * t81 - t24 * t79 + (t109 * t45 + t111 * t44) * t158, t10 * t23 + t11 * t24 + t27 * t71 + t40 * t47 + t44 * t6 + t45 * t7, 0, 0, 0, 0, 0 (-t135 * qJD(6) - t113 * t24 + t116 * t23) * t103 + t136 * t157 + t47 * t34 + t71 * t14 -(t136 * qJD(6) + t113 * t23 + t116 * t24) * t103 - t135 * t157 - t47 * t134 + t71 * t13; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), qJ(3) * t83 + qJD(3) * t88 + (-t196 + (-t86 - t203) * t115) * t185, -0.2e1 * t114 * t157, 0.2e1 * t187 * t172, -t120 * t114, -t120 * t117, 0, t122 * t114 - t124 * t176, t122 * t117 + t124 * t177 (t79 * t160 + t201 + (qJD(2) * t54 + t10) * qJD(4)) * t117 + (t6 + (-t109 * t40 + t119 * t79) * qJD(4) + (t152 + t209) * qJD(2)) * t114 (t81 * t160 + t200 + (-qJD(2) * t55 - t11) * qJD(4)) * t117 + (-t7 + (-t111 * t40 + t119 * t81) * qJD(4) + (t98 - t208) * qJD(2)) * t114, -t209 * t81 - t208 * t79 + (-t109 * t7 - t111 * t6) * t117 + ((t109 * t55 + t111 * t54) * qJD(2) - t140) * t177, t175 * t198 + t54 * t6 + t55 * t7 + (-t119 * t27 + t40 * t160) * t117 + t208 * t11 + t209 * t10, t13 * t64 - t134 * t25, -t13 * t63 + t134 * t26 - t14 * t64 - t25 * t34, t25 * t103 + t13 * t114 + (qJD(2) * t64 - t134) * t176, -t26 * t103 - t14 * t114 + (-qJD(2) * t63 - t34) * t176 (t103 + t181) * t176, t168 * t114 - t68 * t34 + t78 * t14 + t16 * t63 + t21 * t26 + (t210 * t113 + t211 * t116) * t103 + (-t137 * t103 - t2 * t114) * qJD(6) + (t34 * t160 + (t138 * qJD(2) - t143) * qJD(4)) * t117, -t144 * t114 + t68 * t134 + t78 * t13 + t16 * t64 + t21 * t25 + (-t211 * t113 + t210 * t116) * t103 + (-t103 * t138 + t114 * t143) * qJD(6) + (-t134 * t160 + (-qJD(2) * t137 - t2) * qJD(4)) * t117; 0, 0, 0, 0, 0, -t121, t130, 0, 0, 0, 0, 0, t186 * t114, t186 * t117 (-t111 * t121 + (t79 - t166) * qJD(4)) * t114 (t109 * t121 + (t81 - t163) * qJD(4)) * t114 (t109 * t81 - t111 * t79) * t176 + (t109 * t79 + t111 * t81) * qJD(2), -t117 * t27 + t145 * t114 + t140 * qJD(2) + (t139 * t117 + t198) * qJD(4), 0, 0, 0, 0, 0, -t103 * t126 + (-t85 * t179 - t14) * t117 + (-t103 * t216 + (-t85 * t180 + t34) * qJD(4)) * t114, t103 * t128 + (-t179 * t214 - t13) * t117 + (t103 * t127 + (-t117 * t126 - t134) * qJD(4)) * t114; 0, 0, 0, 0, 0, 0, 0, t117 * t121 * t114, -t187 * t121, 0, 0, 0, qJD(4) * t49 + t117 * t130 - t207, -t146 * t181, -t200 - t49 * t79 + (-t10 * t117 + t123 * t109 - t114 * t18) * qJD(2), t201 - t49 * t81 + (t11 * t117 + t123 * t111 + t114 * t19) * qJD(2), t18 * t81 + t19 * t79 + (-qJD(5) * t79 - t10 * t181 + t7) * t111 + (qJD(5) * t81 - t11 * t181 - t6) * t109, -pkin(4) * t27 + t145 * qJ(5) + t139 * qJD(5) - t10 * t18 - t11 * t19 - t40 * t49, t13 * t85 - t134 * t206, t13 * t214 + t134 * t205 - t14 * t85 - t206 * t34, t206 * t103 + (qJD(4) * t85 + t134) * t180, -t205 * t103 + (qJD(4) * t214 + t34) * t180, -t103 * t180, t105 * t14 - t16 * t214 - t31 * t34 + t205 * t21 + (t131 * t113 - t132 * t116) * t103 + ((-t113 * t97 - t116 * t96) * qJD(4) + t143) * t180, t105 * t13 + t16 * t85 + t31 * t134 + t206 * t21 + (t113 * t132 + t116 * t131) * t103 + (-(-t113 * t96 + t116 * t97) * qJD(4) + t2) * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t81 - t178) * t181 (-t79 - t173) * t181, -t79 ^ 2 - t81 ^ 2, t10 * t81 + t11 * t79 + t27, 0, 0, 0, 0, 0, t14 - t218, t13 - t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t34, t134 ^ 2 - t34 ^ 2, t13 + t202, -t14 - t218, t157, t21 * t134 - t219 * t2 + t168, t219 * t143 + t21 * t34 - t144;];
tauc_reg  = t1;
