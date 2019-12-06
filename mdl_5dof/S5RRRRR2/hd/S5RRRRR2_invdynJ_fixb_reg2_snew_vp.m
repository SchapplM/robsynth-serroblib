% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:43
% EndTime: 2019-12-05 18:53:50
% DurationCPUTime: 1.88s
% Computational Cost: add. (5757->233), mult. (8036->334), div. (0->0), fcn. (5980->10), ass. (0->185)
t160 = sin(qJ(3));
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t210 = sin(qJ(1));
t211 = cos(qJ(1));
t171 = t211 * g(1) + t210 * g(2);
t139 = -qJD(1) ^ 2 * pkin(1) - t171;
t161 = sin(qJ(2));
t165 = cos(qJ(2));
t170 = t210 * g(1) - t211 * g(2);
t168 = qJDD(1) * pkin(1) + t170;
t111 = t165 * t139 + t161 * t168;
t164 = cos(qJ(3));
t103 = t164 * g(3) + t160 * t111;
t155 = qJD(1) + qJD(2);
t151 = t155 ^ 2;
t143 = t164 * t151 * t160;
t185 = qJDD(3) + t143;
t167 = pkin(2) * t185 - t103;
t104 = -t160 * g(3) + t164 * t111;
t157 = t164 ^ 2;
t149 = t157 * t151;
t166 = qJD(3) ^ 2;
t215 = -t149 - t166;
t95 = t215 * pkin(2) + t104;
t65 = t159 * t95 - t163 * t167;
t62 = t163 * t65;
t201 = t163 * t95;
t66 = t159 * t167 + t201;
t218 = -t159 * t66 + t62;
t224 = t160 * t218;
t158 = sin(qJ(5));
t193 = t160 * t163;
t131 = (t164 * t159 + t193) * t155;
t153 = qJDD(1) + qJDD(2);
t144 = t160 * t153;
t189 = qJD(3) * t155;
t179 = t164 * t189;
t135 = t144 + t179;
t146 = t164 * t153;
t180 = t160 * t189;
t136 = t146 - t180;
t176 = t159 * t135 - t163 * t136;
t87 = -t131 * qJD(4) - t176;
t84 = qJDD(5) - t87;
t154 = qJD(3) + qJD(4);
t162 = cos(qJ(5));
t116 = t158 * t131 - t162 * t154;
t118 = t162 * t131 + t158 * t154;
t86 = t118 * t116;
t219 = t84 - t86;
t223 = t158 * t219;
t194 = t160 * t159;
t129 = (-t164 * t163 + t194) * t155;
t105 = t131 * t129;
t152 = qJDD(3) + qJDD(4);
t217 = -t105 + t152;
t222 = t159 * t217;
t221 = t162 * t219;
t220 = t163 * t217;
t190 = t160 * t103 + t164 * t104;
t216 = t136 - t180;
t126 = qJD(5) + t129;
t172 = t163 * t135 + t159 * t136;
t88 = -t129 * qJD(4) + t172;
t177 = -t162 * t152 + t158 * t88;
t49 = (qJD(5) - t126) * t118 + t177;
t112 = t116 ^ 2;
t113 = t118 ^ 2;
t125 = t126 ^ 2;
t127 = t129 ^ 2;
t128 = t131 ^ 2;
t150 = t154 ^ 2;
t110 = t161 * t139 - t165 * t168;
t90 = -t216 * pkin(2) + t110;
t41 = t158 * t66 - t162 * t90;
t42 = t158 * t90 + t162 * t66;
t16 = t158 * t42 - t162 * t41;
t173 = -t158 * t152 - t162 * t88;
t70 = -t116 * qJD(5) - t173;
t94 = t126 * t116;
t53 = t70 + t94;
t20 = -t158 * t49 - t162 * t53;
t214 = -t16 * t193 + t164 * (-pkin(2) * t20 - t159 * t16);
t205 = t159 * t65;
t79 = -t125 - t112;
t28 = t158 * t79 + t221;
t213 = t160 * (t158 * t62 - t159 * t41) + t164 * (-pkin(2) * t28 + t158 * t205 + t163 * t41);
t58 = t84 + t86;
t206 = t158 * t58;
t81 = -t113 - t125;
t30 = t162 * t81 - t206;
t61 = t162 * t65;
t212 = t160 * (-t159 * t42 + t162 * t62) + t164 * (-pkin(2) * t30 + t159 * t61 + t163 * t42);
t209 = t164 * pkin(2);
t178 = t163 * t66 + t205;
t89 = -t127 - t128;
t208 = t224 + t164 * (-pkin(2) * t89 + t178);
t18 = t158 * t41 + t162 * t42;
t188 = qJD(4) + t154;
t76 = -t188 * t129 + t172;
t207 = t90 * t193 + t164 * (-pkin(2) * t76 + t159 * t90);
t98 = t105 + t152;
t203 = t159 * t98;
t202 = t162 * t58;
t200 = t163 * t98;
t199 = t126 * t158;
t198 = t126 * t162;
t197 = t154 * t159;
t196 = t154 * t163;
t195 = t160 * t185;
t192 = t161 * t153;
t191 = t164 * (qJDD(3) - t143);
t186 = qJD(5) + t126;
t184 = t16 * t209;
t183 = t90 * t209;
t182 = t159 * t86;
t181 = t163 * t86;
t73 = t188 * t131 + t176;
t175 = t90 * t194 + t164 * (-pkin(2) * t73 - t163 * t90);
t169 = (-qJD(4) + t154) * t131 - t176;
t156 = t160 ^ 2;
t147 = t156 * t151;
t137 = t146 - 0.2e1 * t180;
t134 = t144 + 0.2e1 * t179;
t124 = t154 * t129;
t123 = -t128 + t150;
t122 = t127 - t150;
t121 = -t128 - t150;
t120 = t195 + t164 * (-t147 + t166);
t119 = t160 * (t149 - t166) + t191;
t115 = (t135 + t179) * t160;
t114 = t216 * t164;
t109 = t165 * t110;
t108 = t164 * t110;
t107 = t160 * t110;
t106 = t164 * t134 + t160 * t137;
t102 = t128 - t127;
t96 = -t150 - t127;
t93 = -t113 + t125;
t92 = t112 - t125;
t83 = t113 - t112;
t80 = t163 * t121 - t203;
t78 = t124 + t88;
t77 = -t124 + t88;
t72 = t112 + t113;
t71 = t159 * t96 + t220;
t69 = -t118 * qJD(5) - t177;
t68 = (-t116 * t162 + t118 * t158) * t126;
t67 = (-t116 * t158 - t118 * t162) * t126;
t64 = (t160 * (-t129 * t163 + t131 * t159) + t164 * (-t129 * t159 - t131 * t163)) * t154;
t60 = t158 * t65;
t54 = t186 * t116 + t173;
t52 = t70 - t94;
t50 = -t186 * t118 - t177;
t48 = t160 * (t163 * t122 - t203) + t164 * (t159 * t122 + t200);
t47 = t160 * (-t159 * t123 + t220) + t164 * (t163 * t123 + t222);
t46 = -t118 * t199 + t162 * t70;
t45 = t118 * t198 + t158 * t70;
t44 = t116 * t198 - t158 * t69;
t43 = t116 * t199 + t162 * t69;
t40 = t159 * t169 - t163 * t78;
t37 = t162 * t92 - t206;
t36 = -t158 * t93 + t221;
t35 = t158 * t92 + t202;
t34 = t162 * t93 + t223;
t33 = t160 * (-t131 * t197 + t163 * t88) + t164 * (t131 * t196 + t159 * t88);
t32 = t160 * (t129 * t196 - t159 * t87) + t164 * (t129 * t197 + t163 * t87);
t31 = -t158 * t81 - t202;
t29 = t162 * t79 - t223;
t22 = t158 * t53 - t162 * t49;
t21 = -t158 * t52 + t162 * t50;
t19 = t158 * t50 + t162 * t52;
t15 = t160 * (-t159 * t77 - t163 * t73) + t164 * (-t159 * t73 + t163 * t77);
t13 = t160 * (t159 * t84 + t163 * t68) + t164 * (t159 * t68 - t163 * t84);
t12 = t159 * t31 + t163 * t54;
t11 = t159 * t29 + t163 * t50;
t10 = t159 * t22 + t163 * t72;
t9 = t160 * (t163 * t46 + t182) + t164 * (t159 * t46 - t181);
t8 = t160 * (t163 * t44 - t182) + t164 * (t159 * t44 + t181);
t7 = t159 * t18 - t62;
t3 = t160 * (-t159 * t49 + t163 * t37) + t164 * (t159 * t37 + t163 * t49);
t2 = t160 * (t159 * t53 + t163 * t36) + t164 * (t159 * t36 - t163 * t53);
t1 = t160 * (t159 * t83 + t163 * t21) + t164 * (t159 * t21 - t163 * t83);
t4 = [0, 0, 0, 0, 0, qJDD(1), t170, t171, 0, 0, 0, 0, 0, 0, 0, t153, pkin(1) * (-t161 * t151 + t165 * t153) - t110, pkin(1) * (-t165 * t151 - t192) - t111, 0, pkin(1) * (t161 * t111 - t109), t115, t106, t120, t114, t119, 0, -t108 + pkin(1) * (t161 * (t164 * t215 - t195) + t165 * t137), t107 + pkin(1) * (t161 * (-t191 - t160 * (-t147 - t166)) - t165 * t134), pkin(1) * (t165 * (t147 + t149) + (t156 + t157) * t192) + t190, pkin(1) * (t161 * t190 - t109), t33, t15, t47, t32, t48, t64, pkin(1) * (t161 * (t164 * (t163 * t96 - t222) - t160 * t71) - t165 * t73) + t175, pkin(1) * (t161 * (t164 * (-t159 * t121 - t200) - t160 * t80) - t165 * t76) + t207, pkin(1) * (t161 * (t164 * (t159 * t78 + t163 * t169) - t160 * t40) - t165 * t89) + t208, -t183 + pkin(1) * (t161 * (t164 * t178 + t224) - t165 * t90), t9, t1, t2, t8, t3, t13, pkin(1) * (t161 * (t164 * (-t159 * t50 + t163 * t29) - t160 * t11) - t165 * t28) + t213, pkin(1) * (t161 * (t164 * (-t159 * t54 + t163 * t31) - t160 * t12) - t165 * t30) + t212, pkin(1) * (t161 * (t164 * (-t159 * t72 + t163 * t22) - t160 * t10) - t165 * t20) + t214, -t184 + pkin(1) * (t161 * (t164 * (t163 * t18 + t205) - t160 * t7) - t165 * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t110, -t111, 0, 0, t115, t106, t120, t114, t119, 0, -t108, t107, t190, 0, t33, t15, t47, t32, t48, t64, t175, t207, t208, -t183, t9, t1, t2, t8, t3, t13, t213, t212, t214, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t147 - t149, t144, t143, t146, qJDD(3), -t103, -t104, 0, 0, t105, t102, t78, -t105, t169, t152, pkin(2) * t71 - t65, -t201 + t159 * t103 + (-t159 * t185 + t80) * pkin(2), pkin(2) * t40, -pkin(2) * t218, t45, t19, t34, t43, t35, t67, pkin(2) * t11 - t61, pkin(2) * t12 + t60, pkin(2) * t10 + t18, pkin(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t102, t78, -t105, t169, t152, -t65, -t66, 0, 0, t45, t19, t34, t43, t35, t67, -t61, t60, t18, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t83, t53, -t86, -t49, t84, -t41, -t42, 0, 0;];
tauJ_reg = t4;
