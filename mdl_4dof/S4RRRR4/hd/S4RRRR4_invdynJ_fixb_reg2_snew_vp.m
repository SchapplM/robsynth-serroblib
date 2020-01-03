% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:20
% DurationCPUTime: 1.80s
% Computational Cost: add. (5618->274), mult. (11698->373), div. (0->0), fcn. (7975->8), ass. (0->173)
t145 = sin(qJ(4));
t146 = sin(qJ(3));
t147 = sin(qJ(2));
t150 = cos(qJ(3));
t151 = cos(qJ(2));
t120 = (t146 * t151 + t147 * t150) * qJD(1);
t142 = qJD(2) + qJD(3);
t149 = cos(qJ(4));
t105 = t120 * t145 - t142 * t149;
t107 = t120 * t149 + t142 * t145;
t86 = t107 * t105;
t136 = t147 * qJDD(1);
t181 = qJD(1) * qJD(2);
t173 = t151 * t181;
t125 = t136 + t173;
t137 = t151 * qJDD(1);
t174 = t147 * t181;
t126 = t137 - t174;
t168 = t146 * t125 - t126 * t150;
t88 = -qJD(3) * t120 - t168;
t87 = qJDD(4) - t88;
t209 = -t86 + t87;
t213 = t145 * t209;
t185 = qJD(1) * t147;
t118 = -qJD(1) * t150 * t151 + t146 * t185;
t101 = t120 * t118;
t141 = qJDD(2) + qJDD(3);
t207 = -t101 + t141;
t212 = t146 * t207;
t211 = t149 * t209;
t210 = t150 * t207;
t113 = t142 * t118;
t163 = t150 * t125 + t146 * t126;
t89 = -qJD(3) * t118 + t163;
t208 = -t89 + t113;
t144 = t151 ^ 2;
t154 = qJD(1) ^ 2;
t148 = sin(qJ(1));
t152 = cos(qJ(1));
t171 = t148 * g(1) - g(2) * t152;
t161 = qJDD(1) * pkin(1) + t171;
t162 = qJD(2) * pkin(2) - pkin(6) * t185;
t91 = t126 * pkin(2) + (pkin(6) * t144 + pkin(5)) * t154 - t162 * t185 + t161;
t29 = t208 * pkin(7) + (t120 * t142 - t88) * pkin(3) - t91;
t205 = t142 ^ 2;
t167 = g(1) * t152 + g(2) * t148;
t193 = qJDD(1) * pkin(5);
t157 = -pkin(1) * t154 - t167 + t193;
t110 = -t147 * g(3) + t151 * t157;
t139 = t144 * t154;
t85 = -pkin(2) * t139 + t126 * pkin(6) - qJD(2) * t162 + t110;
t196 = t150 * t85;
t156 = t147 * t157;
t187 = t147 * t154;
t202 = t125 * pkin(6);
t206 = qJDD(2) * pkin(2) - t156 + (pkin(2) * t187 + pkin(6) * t181 - g(3)) * t151 - t202;
t61 = t146 * t206 + t196;
t99 = pkin(3) * t118 - pkin(7) * t120;
t37 = -pkin(3) * t205 + t141 * pkin(7) - t118 * t99 + t61;
t11 = t145 * t37 - t149 * t29;
t12 = t145 * t29 + t149 * t37;
t6 = t145 * t11 + t149 * t12;
t115 = qJD(4) + t118;
t170 = -t141 * t149 + t145 * t89;
t53 = (qJD(4) - t115) * t107 + t170;
t102 = t105 ^ 2;
t103 = t107 ^ 2;
t114 = t115 ^ 2;
t116 = t118 ^ 2;
t117 = t120 ^ 2;
t60 = t146 * t85 - t150 * t206;
t36 = -pkin(3) * t141 - pkin(7) * t205 + t120 * t99 + t60;
t204 = -pkin(3) * t36 + pkin(7) * t6;
t203 = pkin(3) * t146;
t33 = t145 * t36;
t63 = t86 + t87;
t201 = t145 * t63;
t200 = t146 * t91;
t97 = t101 + t141;
t199 = t146 * t97;
t26 = t146 * t61 - t150 * t60;
t198 = t147 * t26;
t34 = t149 * t36;
t197 = t149 * t63;
t195 = t150 * t91;
t194 = t150 * t97;
t192 = t115 * t145;
t191 = t115 * t149;
t190 = t142 * t146;
t189 = t142 * t150;
t131 = t151 * t187;
t128 = qJDD(2) + t131;
t188 = t147 * t128;
t186 = t151 * (qJDD(2) - t131);
t184 = qJD(3) + t142;
t182 = qJD(4) + t115;
t82 = -t103 - t114;
t40 = -t145 * t82 - t197;
t164 = -t145 * t141 - t149 * t89;
t58 = t105 * t182 + t164;
t180 = pkin(3) * t58 + pkin(7) * t40 + t33;
t73 = -t114 - t102;
t32 = t149 * t73 - t213;
t54 = -t107 * t182 - t170;
t179 = pkin(3) * t54 + pkin(7) * t32 - t34;
t178 = t146 * t86;
t177 = t150 * t86;
t176 = -pkin(3) * t150 - pkin(2);
t68 = -qJD(4) * t105 - t164;
t94 = t115 * t105;
t57 = t68 + t94;
t25 = t145 * t57 - t149 * t53;
t70 = t102 + t103;
t172 = pkin(3) * t70 + pkin(7) * t25 + t6;
t27 = t146 * t60 + t150 * t61;
t109 = t151 * g(3) + t156;
t169 = t147 * t109 + t110 * t151;
t165 = -t11 * t149 + t12 * t145;
t160 = (-qJD(3) + t142) * t120 - t168;
t153 = qJD(2) ^ 2;
t143 = t147 ^ 2;
t138 = t143 * t154;
t127 = t137 - 0.2e1 * t174;
t124 = t136 + 0.2e1 * t173;
t121 = t154 * pkin(5) + t161;
t112 = -t117 + t205;
t111 = t116 - t205;
t108 = -t117 - t205;
t100 = t117 - t116;
t95 = -t205 - t116;
t93 = -t103 + t114;
t92 = t102 - t114;
t90 = -t116 - t117;
t84 = t103 - t102;
t81 = -t108 * t146 - t194;
t80 = t108 * t150 - t199;
t79 = t113 + t89;
t77 = -t118 * t184 + t163;
t74 = t120 * t184 + t168;
t72 = t150 * t95 - t212;
t71 = t146 * t95 + t210;
t67 = -qJD(4) * t107 - t170;
t66 = (-t105 * t149 + t107 * t145) * t115;
t65 = (-t105 * t145 - t107 * t149) * t115;
t56 = t68 - t94;
t50 = -t107 * t192 + t149 * t68;
t49 = t107 * t191 + t145 * t68;
t48 = t105 * t191 - t145 * t67;
t47 = t105 * t192 + t149 * t67;
t46 = t146 * t79 + t150 * t160;
t45 = t146 * t160 - t150 * t79;
t44 = t149 * t92 - t201;
t43 = -t145 * t93 + t211;
t42 = t145 * t92 + t197;
t41 = t149 * t93 + t213;
t39 = t149 * t82 - t201;
t31 = t145 * t73 + t211;
t24 = -t145 * t56 + t149 * t54;
t23 = -t145 * t53 - t149 * t57;
t22 = t145 * t54 + t149 * t56;
t20 = -t146 * t58 + t150 * t40;
t19 = t146 * t40 + t150 * t58;
t18 = -t146 * t54 + t150 * t32;
t17 = t146 * t32 + t150 * t54;
t16 = -pkin(7) * t39 + t34;
t15 = -t146 * t70 + t150 * t25;
t14 = t146 * t25 + t150 * t70;
t13 = -pkin(7) * t31 + t33;
t8 = -pkin(3) * t39 + t12;
t7 = -pkin(3) * t31 + t11;
t2 = t146 * t6 - t150 * t36;
t1 = -pkin(7) * t23 - t165;
t3 = [0, 0, 0, 0, 0, qJDD(1), t171, t167, 0, 0, (t125 + t173) * t147, t124 * t151 + t127 * t147, t188 + t151 * (-t138 + t153), (t126 - t174) * t151, t147 * (t139 - t153) + t186, 0, t151 * t121 + pkin(1) * t127 + pkin(5) * (t151 * (-t139 - t153) - t188), -t147 * t121 - pkin(1) * t124 + pkin(5) * (-t186 - t147 * (-t138 - t153)), pkin(1) * (t138 + t139) + (t143 + t144) * t193 + t169, pkin(1) * t121 + pkin(5) * t169, t147 * (-t120 * t190 + t150 * t89) + t151 * (t120 * t189 + t146 * t89), t147 * (t146 * t208 - t150 * t74) + t151 * (-t146 * t74 - t150 * t208), t147 * (-t112 * t146 + t210) + t151 * (t112 * t150 + t212), t147 * (t118 * t189 - t146 * t88) + t151 * (t118 * t190 + t150 * t88), t147 * (t111 * t150 - t199) + t151 * (t111 * t146 + t194), (t147 * (-t118 * t150 + t120 * t146) + t151 * (-t118 * t146 - t120 * t150)) * t142, t147 * (-pkin(6) * t71 - t200) + t151 * (-pkin(2) * t74 + pkin(6) * t72 + t195) - pkin(1) * t74 + pkin(5) * (-t147 * t71 + t151 * t72), t147 * (-pkin(6) * t80 - t195) + t151 * (-pkin(2) * t77 + pkin(6) * t81 - t200) - pkin(1) * t77 + pkin(5) * (-t147 * t80 + t151 * t81), t147 * (-pkin(6) * t45 - t26) + t151 * (-pkin(2) * t90 + pkin(6) * t46 + t27) - pkin(1) * t90 + pkin(5) * (-t147 * t45 + t151 * t46), -pkin(6) * t198 + t151 * (pkin(2) * t91 + pkin(6) * t27) + pkin(1) * t91 + pkin(5) * (t151 * t27 - t198), t147 * (t150 * t50 + t178) + t151 * (t146 * t50 - t177), t147 * (t146 * t84 + t150 * t24) + t151 * (t146 * t24 - t150 * t84), t147 * (t146 * t57 + t150 * t43) + t151 * (t146 * t43 - t150 * t57), t147 * (t150 * t48 - t178) + t151 * (t146 * t48 + t177), t147 * (-t146 * t53 + t150 * t44) + t151 * (t146 * t44 + t150 * t53), t147 * (t146 * t87 + t150 * t66) + t151 * (t146 * t66 - t150 * t87), t147 * (-pkin(6) * t17 + t13 * t150 - t146 * t7) + t151 * (-pkin(2) * t31 + pkin(6) * t18 + t13 * t146 + t150 * t7) - pkin(1) * t31 + pkin(5) * (-t147 * t17 + t151 * t18), t147 * (-pkin(6) * t19 - t146 * t8 + t150 * t16) + t151 * (-pkin(2) * t39 + pkin(6) * t20 + t146 * t16 + t150 * t8) - pkin(1) * t39 + pkin(5) * (-t147 * t19 + t151 * t20), t147 * (-pkin(6) * t14 + t1 * t150) + t151 * (pkin(6) * t15 + t146 * t1) + pkin(5) * (-t14 * t147 + t15 * t151) + (t147 * t203 + t151 * t176 - pkin(1)) * t23, (t147 * (-pkin(7) * t150 + t203) + t151 * (-pkin(7) * t146 + t176) - pkin(1)) * t165 + (pkin(5) + pkin(6)) * (-t147 * t2 + t151 * (t146 * t36 + t150 * t6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t138 - t139, t136, t131, t137, qJDD(2), -t109, -t110, 0, 0, t101, t100, t79, -t101, t160, t141, pkin(2) * t71 - t60, -t196 - t146 * (pkin(6) * t173 - t109 - t202) + (-t128 * t146 + t80) * pkin(2), pkin(2) * t45, pkin(2) * t26, t49, t22, t41, t47, t42, t65, pkin(2) * t17 + t179, pkin(2) * t19 + t180, pkin(2) * t14 + t172, pkin(2) * t2 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t100, t79, -t101, t160, t141, -t60, -t61, 0, 0, t49, t22, t41, t47, t42, t65, t179, t180, t172, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t84, t57, -t86, -t53, t87, -t11, -t12, 0, 0;];
tauJ_reg = t3;
