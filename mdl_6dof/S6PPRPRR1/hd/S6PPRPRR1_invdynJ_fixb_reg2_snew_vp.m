% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:54:04
% EndTime: 2019-05-04 19:54:11
% DurationCPUTime: 2.28s
% Computational Cost: add. (11850->260), mult. (20897->407), div. (0->0), fcn. (17866->16), ass. (0->185)
t146 = cos(pkin(7));
t140 = sin(pkin(11));
t145 = cos(pkin(11));
t116 = -g(1) * t145 - g(2) * t140;
t139 = sin(pkin(12));
t144 = cos(pkin(12));
t115 = g(1) * t140 - g(2) * t145;
t135 = -g(3) + qJDD(1);
t142 = sin(pkin(6));
t147 = cos(pkin(6));
t166 = t115 * t147 + t135 * t142;
t74 = -t116 * t139 + t144 * t166;
t203 = t146 * t74;
t141 = sin(pkin(7));
t94 = -t115 * t142 + t135 * t147 + qJDD(2);
t204 = t141 * t94;
t208 = t203 + t204;
t157 = qJD(3) ^ 2;
t150 = sin(qJ(6));
t151 = sin(qJ(5));
t188 = qJD(3) * qJD(5);
t126 = t151 * t188;
t154 = cos(qJ(5));
t129 = t154 * qJDD(3);
t109 = t129 - t126;
t102 = -qJDD(6) + t109;
t153 = cos(qJ(6));
t191 = qJD(3) * t151;
t103 = -qJD(5) * t153 + t150 * t191;
t105 = qJD(5) * t150 + t153 * t191;
t197 = t105 * t103;
t161 = -t102 - t197;
t207 = t150 * t161;
t206 = t153 * t161;
t123 = qJD(3) * t154 - qJD(6);
t185 = t154 * t188;
t187 = t151 * qJDD(3);
t108 = t185 + t187;
t182 = -qJDD(5) * t153 + t108 * t150;
t63 = (qJD(6) + t123) * t105 + t182;
t100 = t103 ^ 2;
t101 = t105 ^ 2;
t121 = t123 ^ 2;
t138 = sin(pkin(13));
t143 = cos(pkin(13));
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t75 = t116 * t144 + t139 * t166;
t52 = -t152 * t75 + t155 * t208;
t159 = qJDD(3) * pkin(3) + t52;
t53 = t152 * t208 + t155 * t75;
t51 = -pkin(3) * t157 + t53;
t35 = t138 * t159 + t143 * t51;
t205 = t141 * t74;
t202 = t146 * t94;
t156 = qJD(5) ^ 2;
t181 = -pkin(5) * t154 - pkin(10) * t151;
t33 = -pkin(4) * t157 + qJDD(3) * pkin(9) + t35;
t183 = t157 * t181 + t33;
t61 = t202 - t205;
t160 = qJDD(4) + t61;
t58 = t154 * t160;
t22 = -qJDD(5) * pkin(5) - pkin(10) * t156 + t151 * t183 - t58;
t201 = t150 * t22;
t77 = t102 - t197;
t200 = t150 * t77;
t199 = t153 * t22;
t198 = t153 * t77;
t196 = t123 * t150;
t195 = t123 * t153;
t194 = t144 * t146;
t122 = t151 * t157 * t154;
t117 = qJDD(5) + t122;
t193 = t151 * t117;
t118 = qJDD(5) - t122;
t192 = t154 * t118;
t190 = qJD(6) - t123;
t186 = t154 * t197;
t184 = t138 * t51 - t143 * t159;
t158 = t151 * t160;
t23 = -t156 * pkin(5) + qJDD(5) * pkin(10) + t154 * t183 + t158;
t178 = -t109 + t126;
t179 = t108 + t185;
t32 = -qJDD(3) * pkin(4) - pkin(9) * t157 + t184;
t29 = pkin(5) * t178 - pkin(10) * t179 + t32;
t16 = t150 * t23 - t153 * t29;
t17 = t150 * t29 + t153 * t23;
t8 = t150 * t16 + t153 * t17;
t25 = t151 * t33 - t58;
t26 = t154 * t33 + t158;
t14 = t151 * t25 + t154 * t26;
t6 = t151 * t22 + t154 * t8;
t7 = t150 * t17 - t153 * t16;
t2 = t138 * t6 - t143 * t7;
t3 = t138 * t7 + t143 * t6;
t180 = t152 * t3 + t155 * t2;
t10 = t138 * t14 - t143 * t32;
t11 = t138 * t32 + t14 * t143;
t177 = t10 * t155 + t11 * t152;
t20 = t138 * t35 - t143 * t184;
t21 = t138 * t184 + t143 * t35;
t176 = t152 * t21 + t155 * t20;
t162 = -qJDD(5) * t150 - t108 * t153;
t81 = -qJD(6) * t103 - t162;
t97 = t103 * t123;
t67 = t81 - t97;
t55 = t150 * t67 - t153 * t63;
t76 = t100 + t101;
t41 = -t151 * t76 + t154 * t55;
t54 = -t150 * t63 - t153 * t67;
t30 = t138 * t41 - t143 * t54;
t31 = t138 * t54 + t143 * t41;
t175 = t152 * t31 + t155 * t30;
t84 = -t121 - t100;
t57 = t153 * t84 - t207;
t64 = -t105 * t190 - t182;
t43 = -t151 * t64 + t154 * t57;
t56 = t150 * t84 + t206;
t36 = t138 * t43 - t143 * t56;
t37 = t138 * t56 + t143 * t43;
t174 = t152 * t37 + t155 * t36;
t88 = -t101 - t121;
t60 = -t150 * t88 + t198;
t68 = t103 * t190 + t162;
t45 = -t151 * t68 + t154 * t60;
t59 = t153 * t88 + t200;
t38 = t138 * t45 - t143 * t59;
t39 = t138 * t59 + t143 * t45;
t173 = t152 * t39 + t155 * t38;
t172 = t152 * t53 + t155 * t52;
t110 = t129 - 0.2e1 * t126;
t134 = t154 ^ 2;
t132 = t134 * t157;
t120 = -t132 - t156;
t92 = t120 * t154 - t193;
t70 = t110 * t143 + t138 * t92;
t72 = -t110 * t138 + t143 * t92;
t171 = t152 * t72 + t155 * t70;
t107 = 0.2e1 * t185 + t187;
t133 = t151 ^ 2;
t130 = t133 * t157;
t119 = -t130 - t156;
t93 = -t119 * t151 - t192;
t71 = -t107 * t143 + t138 * t93;
t73 = t107 * t138 + t143 * t93;
t170 = t152 * t73 + t155 * t71;
t113 = (t133 + t134) * qJDD(3);
t114 = t130 + t132;
t86 = t113 * t138 + t114 * t143;
t87 = t113 * t143 - t114 * t138;
t169 = t152 * t87 + t155 * t86;
t111 = qJDD(3) * t143 - t138 * t157;
t112 = -qJDD(3) * t138 - t143 * t157;
t168 = t111 * t155 + t112 * t152;
t167 = -t111 * t152 + t112 * t155;
t165 = -pkin(4) + t181;
t164 = qJDD(3) * t155 - t152 * t157;
t163 = -qJDD(3) * t152 - t155 * t157;
t99 = t164 * t141;
t98 = t163 * t141;
t96 = -t101 + t121;
t95 = t100 - t121;
t91 = -t118 * t151 + t119 * t154;
t90 = t117 * t154 + t120 * t151;
t89 = t101 - t100;
t83 = t168 * t141;
t82 = t167 * t141;
t80 = -qJD(6) * t105 - t182;
t66 = t81 + t97;
t62 = t169 * t141;
t47 = t141 * t170 + t146 * t91;
t46 = t141 * t171 + t146 * t90;
t44 = t151 * t60 + t154 * t68;
t42 = t151 * t57 + t154 * t64;
t40 = t151 * t55 + t154 * t76;
t27 = t141 * t172 + t146 * t61;
t19 = t141 * t173 + t146 * t44;
t18 = t141 * t174 + t146 * t42;
t13 = t151 * t26 - t154 * t25;
t12 = t141 * t175 + t146 * t40;
t9 = t146 * (qJDD(4) + t202) + (t176 - t203) * t141;
t5 = t151 * t8 - t154 * t22;
t4 = t13 * t146 + t141 * t177;
t1 = t141 * t180 + t146 * t5;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t94 + (t139 * t75 + t144 * t74) * t142, 0, 0, 0, 0, 0, 0, t147 * t99 + (t139 * t163 + t164 * t194) * t142, t147 * t98 + (-t139 * t164 + t163 * t194) * t142, 0, t147 * t27 + (t139 * (-t152 * t52 + t155 * t53) + t144 * (-t141 * t61 + t146 * t172)) * t142, 0, 0, 0, 0, 0, 0, t147 * t83 + (t139 * t167 + t168 * t194) * t142, t147 * t82 + (-t139 * t168 + t167 * t194) * t142, 0, t147 * t9 + (t139 * (-t152 * t20 + t155 * t21) + t144 * (-t141 * (qJDD(4) - t205) + (t176 - t204) * t146)) * t142, 0, 0, 0, 0, 0, 0, t147 * t46 + (t139 * (-t152 * t70 + t155 * t72) + t144 * (-t141 * t90 + t146 * t171)) * t142, t147 * t47 + (t139 * (-t152 * t71 + t155 * t73) + t144 * (-t141 * t91 + t146 * t170)) * t142, t147 * t62 + (t139 * (-t152 * t86 + t155 * t87) + t169 * t194) * t142, t147 * t4 + (t139 * (-t10 * t152 + t11 * t155) + t144 * (-t13 * t141 + t146 * t177)) * t142, 0, 0, 0, 0, 0, 0, t147 * t18 + (t139 * (-t152 * t36 + t155 * t37) + t144 * (-t141 * t42 + t146 * t174)) * t142, t147 * t19 + (t139 * (-t152 * t38 + t155 * t39) + t144 * (-t141 * t44 + t146 * t173)) * t142, t147 * t12 + (t139 * (-t152 * t30 + t155 * t31) + t144 * (-t141 * t40 + t146 * t175)) * t142, t147 * t1 + (t139 * (-t152 * t2 + t155 * t3) + t144 * (-t141 * t5 + t146 * t180)) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, 0, 0, 0, 0, 0, t99, t98, 0, t27, 0, 0, 0, 0, 0, 0, t83, t82, 0, t9, 0, 0, 0, 0, 0, 0, t46, t47, t62, t4, 0, 0, 0, 0, 0, 0, t18, t19, t12, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t52, -t53, 0, 0, 0, 0, 0, 0, 0, qJDD(3), pkin(3) * t111 - t184, pkin(3) * t112 - t35, 0, pkin(3) * t20, t179 * t151, t107 * t154 + t110 * t151, t193 + t154 * (-t130 + t156), -t178 * t154, t151 * (t132 - t156) + t192, 0, pkin(3) * t70 + pkin(4) * t110 + pkin(9) * t92 - t154 * t32, pkin(3) * t71 - pkin(4) * t107 + pkin(9) * t93 + t151 * t32, pkin(3) * t86 + pkin(4) * t114 + pkin(9) * t113 + t14, pkin(3) * t10 - pkin(4) * t32 + pkin(9) * t14, t151 * (t105 * t196 + t153 * t81) - t186, t151 * (-t150 * t66 + t153 * t64) - t154 * t89, t151 * (-t150 * t96 + t206) - t154 * t67, t151 * (-t103 * t195 - t150 * t80) + t186, t151 * (t153 * t95 + t200) + t154 * t63, t154 * t102 + t151 * (t103 * t153 - t105 * t150) * t123, t151 * (-pkin(10) * t56 + t201) + t154 * (-pkin(5) * t56 + t16) - pkin(4) * t56 + pkin(9) * t43 + pkin(3) * t36, t151 * (-pkin(10) * t59 + t199) + t154 * (-pkin(5) * t59 + t17) - pkin(4) * t59 + pkin(9) * t45 + pkin(3) * t38, pkin(3) * t30 + pkin(9) * t41 - t151 * t7 + t165 * t54, pkin(3) * t2 + pkin(9) * t6 + t165 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, 0, 0, 0, 0, 0, 0, t90, t91, 0, t13, 0, 0, 0, 0, 0, 0, t42, t44, t40, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t130 - t132, t187, t122, t129, qJDD(5), -t25, -t26, 0, 0, -t105 * t195 + t150 * t81, t150 * t64 + t153 * t66, t153 * t96 + t207, -t103 * t196 + t153 * t80, t150 * t95 - t198, (t103 * t150 + t105 * t153) * t123, pkin(5) * t64 + pkin(10) * t57 - t199, pkin(5) * t68 + pkin(10) * t60 + t201, pkin(5) * t76 + pkin(10) * t55 + t8, -pkin(5) * t22 + pkin(10) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t89, t67, -t197, -t63, -t102, -t16, -t17, 0, 0;];
tauJ_reg  = t15;
