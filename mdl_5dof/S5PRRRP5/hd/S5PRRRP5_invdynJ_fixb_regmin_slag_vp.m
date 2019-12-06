% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:21
% EndTime: 2019-12-05 16:49:26
% DurationCPUTime: 1.39s
% Computational Cost: add. (1315->220), mult. (2952->299), div. (0->0), fcn. (2092->10), ass. (0->138)
t100 = cos(qJ(4));
t101 = cos(qJ(3));
t99 = sin(qJ(2));
t154 = t99 * qJD(1);
t72 = qJD(2) * pkin(6) + t154;
t132 = pkin(7) * qJD(2) + t72;
t40 = t132 * t101;
t97 = sin(qJ(4));
t29 = t97 * t40;
t98 = sin(qJ(3));
t39 = t132 * t98;
t32 = qJD(3) * pkin(3) - t39;
t136 = t100 * t32 - t29;
t57 = t100 * t98 + t97 * t101;
t47 = t57 * qJD(2);
t161 = t47 * qJ(5);
t197 = t161 - t136;
t145 = t101 * qJDD(2);
t148 = t98 * qJDD(2);
t91 = qJD(3) + qJD(4);
t196 = t91 * t57;
t17 = qJD(2) * t196 - t100 * t145 + t148 * t97;
t102 = cos(qJ(2));
t95 = sin(pkin(8));
t96 = cos(pkin(8));
t129 = g(1) * t96 + g(2) * t95;
t120 = t129 * t102;
t176 = g(3) * t99;
t195 = t120 + t176;
t146 = qJD(2) * qJD(3);
t133 = t101 * t146;
t194 = -t133 - t148;
t172 = g(3) * t102;
t191 = t129 * t99;
t193 = t191 - t172;
t153 = qJDD(1) - g(3);
t192 = t102 * t153 + t191;
t182 = pkin(6) + pkin(7);
t66 = t182 * t98;
t67 = t182 * t101;
t167 = t100 * t67 - t97 * t66;
t155 = t100 * t101;
t169 = t97 * t98;
t123 = t155 - t169;
t42 = t123 * t99;
t117 = t123 * t102;
t150 = qJD(4) * t100;
t157 = qJD(4) * t97;
t140 = qJD(3) * t182;
t62 = t98 * t140;
t63 = t101 * t140;
t190 = qJD(1) * t117 + t100 * t62 + t66 * t150 + t157 * t67 + t97 * t63;
t149 = t102 * qJD(1);
t189 = -t167 * qJD(4) - t100 * t63 + t57 * t149 + t97 * t62;
t163 = t102 * t96;
t164 = t102 * t95;
t94 = qJ(3) + qJ(4);
t86 = sin(t94);
t87 = cos(t94);
t187 = t86 * t176 - g(2) * (-t164 * t86 - t96 * t87) - g(1) * (-t163 * t86 + t95 * t87);
t103 = qJD(3) ^ 2;
t144 = t102 * qJDD(1);
t147 = qJD(1) * qJD(2);
t81 = t99 * t147;
t186 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t103 + (t129 + t147) * t99 + t144 - t172 - t81;
t104 = qJD(2) ^ 2;
t143 = t102 * qJDD(2);
t185 = (t103 + t104) * t99 - t143;
t151 = qJD(3) * t101;
t53 = qJDD(2) * pkin(6) + t99 * qJDD(1) + t102 * t147;
t19 = qJDD(3) * pkin(3) + t194 * pkin(7) - t72 * t151 - t98 * t53;
t139 = t98 * t146;
t158 = qJD(3) * t98;
t20 = -t72 * t158 + t101 * t53 + (-t139 + t145) * pkin(7);
t184 = -(qJD(4) * t32 + t20) * t100 + t40 * t157 - t97 * t19;
t183 = t47 ^ 2;
t8 = t91 * pkin(4) - t197;
t181 = t8 + t197;
t88 = t101 * pkin(3);
t175 = pkin(2) + t88;
t174 = -qJ(5) * t196 + qJD(5) * t123 - t190;
t126 = t91 * t169;
t24 = -t100 * t151 - t101 * t150 + t126;
t173 = t24 * qJ(5) - t57 * qJD(5) + t189;
t138 = qJD(2) * t155;
t160 = qJD(2) * t98;
t142 = t97 * t160;
t45 = -t138 + t142;
t54 = -qJD(2) * t175 - t149;
t23 = t45 * pkin(4) + qJD(5) + t54;
t171 = t23 * t47;
t170 = t47 * t45;
t168 = -t100 * t39 - t29;
t65 = pkin(4) * t87 + t88;
t92 = t98 ^ 2;
t166 = -t101 ^ 2 + t92;
t165 = qJD(2) * pkin(2);
t31 = t100 * t40;
t162 = t45 * qJ(5);
t159 = qJD(2) * t99;
t135 = t97 * t39 - t31;
t134 = -t100 * t66 - t97 * t67;
t131 = -qJD(4) * t138 + t194 * t100 - t97 * t145;
t128 = g(1) * t95 - g(2) * t96;
t125 = -t97 * t32 - t31;
t119 = pkin(3) * t158 - t154;
t116 = -qJDD(3) * t99 - 0.2e1 * t102 * t146;
t115 = pkin(3) * t139 - qJDD(2) * t175 + t81;
t112 = qJD(4) * t125 + t100 * t19 - t97 * t20;
t110 = t17 * pkin(4) + qJDD(5) + t115;
t73 = -t149 - t165;
t109 = -pkin(6) * qJDD(3) + (t149 + t73 - t165) * qJD(3);
t108 = -t73 * qJD(2) + t195 - t53;
t106 = -g(1) * (-t163 * t87 - t95 * t86) - g(2) * (-t164 * t87 + t96 * t86) + t54 * t45 + t87 * t176 + t184;
t105 = -t54 * t47 + t112 + t187;
t90 = -qJ(5) - t182;
t89 = qJDD(3) + qJDD(4);
t84 = t100 * pkin(3) + pkin(4);
t64 = -t98 * pkin(3) - pkin(4) * t86;
t61 = pkin(2) + t65;
t44 = t45 ^ 2;
t41 = t57 * t99;
t26 = t115 - t144;
t22 = qJ(5) * t123 + t167;
t21 = -t57 * qJ(5) + t134;
t18 = -t44 + t183;
t16 = qJD(2) * t126 + t131;
t14 = -t102 * t47 - t91 * t42;
t13 = qJD(2) * t117 - t196 * t99;
t12 = -t161 + t168;
t11 = t135 + t162;
t10 = -t125 - t162;
t7 = t47 * t91 - t17;
t6 = -t131 + (-t142 + t45) * t91;
t5 = t110 - t144;
t2 = -t17 * qJ(5) - t45 * qJD(5) - t184;
t1 = t89 * pkin(4) + t16 * qJ(5) - t47 * qJD(5) + t112;
t3 = [t153, 0, -t104 * t99 + t143, -qJDD(2) * t99 - t104 * t102, 0, 0, 0, 0, 0, -t185 * t101 + t116 * t98, t116 * t101 + t185 * t98, 0, 0, 0, 0, 0, -t102 * t17 + t14 * t91 + t159 * t45 - t41 * t89, t102 * t16 - t13 * t91 + t159 * t47 - t42 * t89, -t13 * t45 - t14 * t47 - t41 * t16 - t42 * t17, -t1 * t41 + t10 * t13 - t5 * t102 + t8 * t14 + t159 * t23 + t2 * t42 - g(3); 0, qJDD(2), t192, -t153 * t99 + t120, t92 * qJDD(2) + 0.2e1 * t98 * t133, 0.2e1 * t98 * t145 - 0.2e1 * t146 * t166, qJDD(3) * t98 + t103 * t101, qJDD(3) * t101 - t103 * t98, 0, t186 * t101 + t109 * t98, t109 * t101 - t186 * t98, -t16 * t57 - t47 * t24, -t123 * t16 - t57 * t17 - t196 * t47 + t24 * t45, -t24 * t91 + t57 * t89, t123 * t89 - t196 * t91, 0, t119 * t45 - t26 * t123 + t134 * t89 - t175 * t17 + t189 * t91 + t193 * t87 + t196 * t54, t119 * t47 + t175 * t16 - t167 * t89 + t190 * t91 - t193 * t86 - t54 * t24 + t26 * t57, -t1 * t57 - t10 * t196 + t123 * t2 + t21 * t16 - t22 * t17 - t173 * t47 - t174 * t45 + t8 * t24 - t195, t2 * t22 + t1 * t21 + t5 * (-pkin(4) * t123 - t175) - g(3) * (t102 * t61 - t99 * t90) + t173 * t8 + (pkin(4) * t196 + t119) * t23 + t174 * t10 + t129 * (t102 * t90 + t61 * t99); 0, 0, 0, 0, -t98 * t104 * t101, t166 * t104, t148, t145, qJDD(3), -t101 * t128 + t108 * t98, t101 * t108 + t128 * t98, t170, t18, t6, t7, t89, -t135 * t91 + (t100 * t89 - t157 * t91 - t160 * t45) * pkin(3) + t105, t168 * t91 + (-t150 * t91 - t160 * t47 - t97 * t89) * pkin(3) + t106, t84 * t16 + (t10 + t11) * t47 + (t12 - t8) * t45 + (-t17 * t97 + (-t100 * t45 + t47 * t97) * qJD(4)) * pkin(3), t1 * t84 - t10 * t12 - t8 * t11 - pkin(4) * t171 - g(1) * (t163 * t64 + t95 * t65) - g(2) * (t164 * t64 - t96 * t65) - t64 * t176 + (-t23 * t160 + t2 * t97 + (t10 * t100 - t8 * t97) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t18, t6, t7, t89, -t125 * t91 + t105, t136 * t91 + t106, pkin(4) * t16 - t181 * t45, t181 * t10 + (t1 - t171 + t187) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 - t183, t10 * t45 + t8 * t47 + t110 - t192;];
tau_reg = t3;
