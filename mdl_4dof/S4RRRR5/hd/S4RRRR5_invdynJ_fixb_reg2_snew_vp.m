% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:14
% EndTime: 2019-12-31 17:28:19
% DurationCPUTime: 1.80s
% Computational Cost: add. (6214->271), mult. (12594->364), div. (0->0), fcn. (8484->8), ass. (0->171)
t145 = sin(qJ(4));
t147 = sin(qJ(2));
t173 = qJD(1) * qJD(2);
t134 = t147 * t173;
t151 = cos(qJ(2));
t172 = t151 * qJDD(1);
t123 = -t134 + t172;
t116 = -qJDD(3) + t123;
t112 = -qJDD(4) + t116;
t146 = sin(qJ(3));
t150 = cos(qJ(3));
t177 = qJD(1) * t147;
t117 = -t150 * qJD(2) + t146 * t177;
t119 = t146 * qJD(2) + t150 * t177;
t149 = cos(qJ(4));
t97 = t149 * t117 + t145 * t119;
t99 = -t145 * t117 + t149 * t119;
t68 = t99 * t97;
t204 = -t112 - t68;
t207 = t145 * t204;
t206 = t149 * t204;
t153 = qJD(1) ^ 2;
t135 = t147 * qJDD(1);
t169 = t151 * t173;
t122 = t135 + t169;
t165 = -t150 * qJDD(2) + t146 * t122;
t92 = -t119 * qJD(3) - t165;
t157 = -t146 * qJDD(2) - t150 * t122;
t93 = -t117 * qJD(3) - t157;
t52 = -t97 * qJD(4) + t145 * t92 + t149 * t93;
t132 = t151 * qJD(1) - qJD(3);
t128 = -qJD(4) + t132;
t84 = t97 * t128;
t205 = t52 + t84;
t108 = t117 * t132;
t203 = -t93 + t108;
t184 = t119 * t117;
t154 = -t116 - t184;
t202 = t146 * t154;
t201 = t150 * t154;
t167 = t145 * t93 - t149 * t92;
t37 = (qJD(4) + t128) * t99 + t167;
t72 = (qJD(3) + t132) * t119 + t165;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t114 = t117 ^ 2;
t115 = t119 ^ 2;
t127 = t128 ^ 2;
t130 = t132 ^ 2;
t200 = qJD(2) ^ 2;
t148 = sin(qJ(1));
t152 = cos(qJ(1));
t168 = t148 * g(1) - t152 * g(2);
t110 = qJDD(1) * pkin(1) + t153 * pkin(5) + t168;
t160 = -t123 + t134;
t161 = t122 + t169;
t71 = pkin(2) * t160 - pkin(6) * t161 - t110;
t162 = t152 * g(1) + t148 * g(2);
t185 = qJDD(1) * pkin(5);
t111 = -t153 * pkin(1) - t162 + t185;
t163 = -t151 * pkin(2) - t147 * pkin(6);
t164 = t153 * t163 + t111;
t196 = t147 * g(3);
t81 = -t200 * pkin(2) + qJDD(2) * pkin(6) + t151 * t164 - t196;
t54 = t146 * t81 - t150 * t71;
t28 = t154 * pkin(3) + t203 * pkin(7) - t54;
t105 = -t132 * pkin(3) - t119 * pkin(7);
t55 = t146 * t71 + t150 * t81;
t29 = -t114 * pkin(3) + t92 * pkin(7) + t132 * t105 + t55;
t13 = t145 * t29 - t149 * t28;
t14 = t145 * t28 + t149 * t29;
t6 = -t149 * t13 + t145 * t14;
t199 = pkin(3) * t6;
t40 = t52 - t84;
t18 = -t145 * t37 - t149 * t40;
t198 = pkin(3) * t18;
t197 = t146 * t6;
t195 = t150 * t6;
t194 = t151 * g(3);
t80 = -qJDD(2) * pkin(2) - t200 * pkin(6) + t147 * t164 + t194;
t44 = -t92 * pkin(3) - t114 * pkin(7) + t119 * t105 + t80;
t193 = t145 * t44;
t59 = t112 - t68;
t192 = t145 * t59;
t191 = t146 * t80;
t87 = t116 - t184;
t190 = t146 * t87;
t189 = t149 * t44;
t188 = t149 * t59;
t187 = t150 * t80;
t186 = t150 * t87;
t183 = t128 * t145;
t182 = t128 * t149;
t181 = t132 * t146;
t180 = t132 * t150;
t131 = t151 * t153 * t147;
t179 = t147 * (qJDD(2) + t131);
t178 = t151 * (qJDD(2) - t131);
t176 = qJD(3) - t132;
t171 = t151 * t68;
t170 = t151 * t184;
t7 = t145 * t13 + t149 * t14;
t26 = t146 * t54 + t150 * t55;
t103 = t147 * t111 + t194;
t104 = t151 * t111 - t196;
t166 = t147 * t103 + t151 * t104;
t159 = t146 * t55 - t150 * t54;
t158 = -pkin(1) + t163;
t64 = -t127 - t95;
t30 = t145 * t64 + t206;
t156 = pkin(3) * t30 - t13;
t78 = -t96 - t127;
t42 = t149 * t78 + t192;
t155 = pkin(3) * t42 - t14;
t142 = t151 ^ 2;
t141 = t147 ^ 2;
t139 = t142 * t153;
t137 = t141 * t153;
t124 = -0.2e1 * t134 + t172;
t121 = t135 + 0.2e1 * t169;
t107 = -t115 + t130;
t106 = t114 - t130;
t101 = t115 - t114;
t100 = -t115 - t130;
t94 = -t130 - t114;
t86 = t114 + t115;
t83 = -t96 + t127;
t82 = t95 - t127;
t77 = t117 * t176 + t157;
t75 = t108 + t93;
t73 = -t119 * t176 - t165;
t67 = t96 - t95;
t66 = -t146 * t100 + t186;
t65 = t150 * t100 + t190;
t63 = t150 * t94 - t202;
t62 = t146 * t94 + t201;
t58 = (-t145 * t99 + t149 * t97) * t128;
t57 = (t145 * t97 + t149 * t99) * t128;
t56 = -t95 - t96;
t51 = -t99 * qJD(4) - t167;
t50 = -t146 * t203 - t150 * t72;
t48 = t149 * t82 + t192;
t47 = -t145 * t83 + t206;
t46 = t145 * t82 - t188;
t45 = t149 * t83 + t207;
t43 = -t145 * t78 + t188;
t36 = (qJD(4) - t128) * t99 + t167;
t35 = t149 * t52 + t99 * t183;
t34 = t145 * t52 - t99 * t182;
t33 = -t145 * t51 - t97 * t182;
t32 = t149 * t51 - t97 * t183;
t31 = t149 * t64 - t207;
t24 = -pkin(7) * t42 + t189;
t23 = -t146 * t42 + t150 * t43;
t22 = t146 * t43 + t150 * t42;
t21 = -pkin(7) * t30 + t193;
t20 = t145 * t40 - t149 * t37;
t19 = -t145 * t205 - t149 * t36;
t17 = -t145 * t36 + t149 * t205;
t16 = -t146 * t30 + t150 * t31;
t15 = t146 * t31 + t150 * t30;
t11 = -pkin(3) * t205 + pkin(7) * t43 + t193;
t10 = -pkin(3) * t36 + pkin(7) * t31 - t189;
t9 = -t146 * t18 + t150 * t20;
t8 = t146 * t20 + t150 * t18;
t5 = -pkin(3) * t44 + pkin(7) * t7;
t4 = -pkin(7) * t18 - t6;
t3 = -pkin(3) * t56 + pkin(7) * t20 + t7;
t2 = t150 * t7 - t197;
t1 = t146 * t7 + t195;
t12 = [0, 0, 0, 0, 0, qJDD(1), t168, t162, 0, 0, t161 * t147, t151 * t121 + t147 * t124, t179 + t151 * (-t137 + t200), -t160 * t151, t147 * (t139 - t200) + t178, 0, t151 * t110 + pkin(1) * t124 + pkin(5) * (t151 * (-t139 - t200) - t179), -t147 * t110 - pkin(1) * t121 + pkin(5) * (-t178 - t147 * (-t137 - t200)), pkin(1) * (t137 + t139) + (t141 + t142) * t185 + t166, pkin(1) * t110 + pkin(5) * t166, t147 * (t119 * t181 + t150 * t93) - t170, t147 * (-t146 * t75 + t150 * t73) - t151 * t101, t147 * (-t146 * t107 + t201) + t151 * t203, t147 * (-t117 * t180 - t146 * t92) + t170, t147 * (t150 * t106 + t190) + t151 * t72, t151 * t116 + t147 * (t117 * t150 - t119 * t146) * t132, t147 * (-pkin(6) * t62 + t191) + t151 * (-pkin(2) * t62 + t54) - pkin(1) * t62 + pkin(5) * (-t147 * t73 + t151 * t63), t147 * (-pkin(6) * t65 + t187) + t151 * (-pkin(2) * t65 + t55) - pkin(1) * t65 + pkin(5) * (-t147 * t77 + t151 * t66), -t147 * t159 + pkin(5) * (-t147 * t86 + t151 * t50) + t158 * (-t146 * t72 + t150 * t203), pkin(5) * (t147 * t80 + t151 * t26) + t158 * t159, t147 * (-t146 * t34 + t150 * t35) - t171, t147 * (-t146 * t17 + t150 * t19) - t151 * t67, t147 * (-t146 * t45 + t150 * t47) - t151 * t40, t147 * (-t146 * t32 + t150 * t33) + t171, t147 * (-t146 * t46 + t150 * t48) + t151 * t37, t147 * (-t146 * t57 + t150 * t58) + t151 * t112, t147 * (-pkin(6) * t15 - t146 * t10 + t150 * t21) + t151 * (-pkin(2) * t15 - t156) - pkin(1) * t15 + pkin(5) * (t147 * t36 + t151 * t16), t147 * (-pkin(6) * t22 - t146 * t11 + t150 * t24) + t151 * (-pkin(2) * t22 - t155) - pkin(1) * t22 + pkin(5) * (t147 * t205 + t151 * t23), t147 * (-pkin(6) * t8 - t146 * t3 + t150 * t4) + t151 * (-pkin(2) * t8 - t198) - pkin(1) * t8 + pkin(5) * (t147 * t56 + t151 * t9), t147 * (-pkin(6) * t1 - pkin(7) * t195 - t146 * t5) + t151 * (-pkin(2) * t1 - t199) - pkin(1) * t1 + pkin(5) * (t147 * t44 + t151 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t137 - t139, t135, t131, t172, qJDD(2), -t103, -t104, 0, 0, -t119 * t180 + t146 * t93, t146 * t73 + t150 * t75, t150 * t107 + t202, -t117 * t181 + t150 * t92, t146 * t106 - t186, (t117 * t146 + t119 * t150) * t132, pkin(2) * t73 + pkin(6) * t63 - t187, pkin(2) * t77 + pkin(6) * t66 + t191, pkin(2) * t86 + pkin(6) * t50 + t26, -pkin(2) * t80 + pkin(6) * t26, t146 * t35 + t150 * t34, t146 * t19 + t150 * t17, t146 * t47 + t150 * t45, t146 * t33 + t150 * t32, t146 * t48 + t150 * t46, t146 * t58 + t150 * t57, -pkin(2) * t36 + pkin(6) * t16 + t150 * t10 + t146 * t21, -pkin(2) * t205 + pkin(6) * t23 + t150 * t11 + t146 * t24, -pkin(2) * t56 + pkin(6) * t9 + t146 * t4 + t150 * t3, -pkin(2) * t44 + pkin(6) * t2 - pkin(7) * t197 + t150 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t101, -t203, -t184, -t72, -t116, -t54, -t55, 0, 0, t68, t67, t40, -t68, -t37, -t112, t156, t155, t198, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, t40, -t68, -t37, -t112, -t13, -t14, 0, 0;];
tauJ_reg = t12;
