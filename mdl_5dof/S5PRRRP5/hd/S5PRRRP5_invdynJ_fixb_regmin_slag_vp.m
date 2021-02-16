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
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 16:33:42
% EndTime: 2021-01-15 16:33:50
% DurationCPUTime: 1.65s
% Computational Cost: add. (1752->261), mult. (3879->336), div. (0->0), fcn. (2776->10), ass. (0->159)
t200 = cos(qJ(4));
t108 = sin(qJ(4));
t111 = cos(qJ(3));
t110 = sin(qJ(2));
t169 = t110 * qJD(1);
t81 = qJD(2) * pkin(6) + t169;
t150 = pkin(7) * qJD(2) + t81;
t45 = t150 * t111;
t33 = t108 * t45;
t109 = sin(qJ(3));
t44 = t150 * t109;
t36 = qJD(3) * pkin(3) - t44;
t149 = t200 * t36 - t33;
t62 = t108 * t111 + t200 * t109;
t52 = t62 * qJD(2);
t183 = t52 * qJ(5);
t11 = t149 - t183;
t146 = qJDD(2) * t200;
t164 = t109 * qJDD(2);
t102 = qJD(3) + qJD(4);
t217 = t102 * t62;
t19 = qJD(2) * t217 + t108 * t164 - t111 * t146;
t112 = cos(qJ(2));
t106 = sin(pkin(8));
t107 = cos(pkin(8));
t142 = g(1) * t107 + g(2) * t106;
t135 = t142 * t112;
t196 = g(3) * t110;
t216 = t135 + t196;
t100 = qJDD(3) + qJDD(4);
t154 = qJD(4) * t200;
t199 = pkin(3) * t102;
t215 = -t108 * pkin(3) * t100 - t154 * t199;
t168 = t112 * qJD(1);
t206 = pkin(6) + pkin(7);
t73 = t206 * t109;
t74 = t206 * t111;
t191 = -t108 * t73 + t200 * t74;
t157 = qJD(3) * t206;
t67 = t109 * t157;
t68 = t111 * t157;
t214 = -t191 * qJD(4) + t108 * t67 + t62 * t168 - t200 * t68;
t156 = t200 * t111;
t178 = t108 * t109;
t133 = t156 - t178;
t130 = t112 * t133;
t170 = qJD(4) * t108;
t213 = qJD(1) * t130 + t108 * t68 + t73 * t154 + t74 * t170 + t200 * t67;
t166 = qJD(2) * qJD(3);
t152 = t111 * t166;
t167 = qJD(1) * qJD(2);
t58 = qJDD(2) * pkin(6) + t110 * qJDD(1) + t112 * t167;
t21 = -t111 * t81 * qJD(3) + qJDD(3) * pkin(3) - t109 * t58 + (-t152 - t164) * pkin(7);
t153 = t109 * t166;
t163 = t111 * qJDD(2);
t171 = qJD(3) * t109;
t22 = -t81 * t171 + t111 * t58 + (-t153 + t163) * pkin(7);
t212 = -t108 * t22 + t200 * t21;
t139 = t102 * t178;
t143 = qJD(2) * t156;
t145 = -t102 * t143 - t108 * t163 - t109 * t146;
t18 = qJD(2) * t139 + t145;
t187 = t18 * qJ(5);
t96 = t100 * pkin(4);
t211 = t187 + t96;
t47 = t133 * t110;
t179 = t107 * t112;
t180 = t106 * t112;
t105 = qJ(3) + qJ(4);
t97 = sin(t105);
t98 = cos(t105);
t210 = t97 * t196 - g(2) * (-t107 * t98 - t97 * t180) - g(1) * (t106 * t98 - t97 * t179);
t113 = qJD(3) ^ 2;
t162 = t112 * qJDD(1);
t195 = g(3) * t112;
t91 = t110 * t167;
t209 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t113 + (t142 + t167) * t110 + t162 - t195 - t91;
t176 = qJDD(1) - g(3);
t208 = t142 * t110 + t176 * t112;
t207 = t52 ^ 2;
t99 = t111 * pkin(3);
t203 = pkin(2) + t99;
t202 = -qJ(5) * t217 + qJD(5) * t133 - t213;
t27 = -qJD(3) * t156 - t111 * t154 + t139;
t201 = t27 * qJ(5) - t62 * qJD(5) + t214;
t173 = qJD(2) * t109;
t155 = t108 * t173;
t50 = -t143 + t155;
t194 = t52 * t50;
t10 = t102 * pkin(4) + t11;
t193 = -t11 + t10;
t192 = -t200 * t44 - t33;
t72 = pkin(4) * t98 + t99;
t190 = qJD(2) * pkin(2);
t188 = t110 * t98;
t186 = t19 * qJ(5);
t185 = t50 * qJ(5);
t184 = t50 * t102;
t182 = t52 * t102;
t151 = t50 * pkin(4) + qJD(5);
t59 = -qJD(2) * t203 - t168;
t26 = t151 + t59;
t177 = qJD(5) + t26;
t103 = t109 ^ 2;
t175 = -t111 ^ 2 + t103;
t114 = qJD(2) ^ 2;
t174 = t113 + t114;
t172 = qJD(2) * t110;
t165 = qJDD(3) * t109;
t161 = t112 * qJDD(2);
t160 = pkin(3) * t171;
t159 = pkin(3) * t173;
t35 = t200 * t45;
t148 = t108 * t44 - t35;
t147 = -t108 * t74 - t200 * t73;
t25 = pkin(4) * t217 + t160;
t144 = t25 - t169;
t141 = g(1) * t106 - g(2) * t107;
t138 = t142 * t188 - t98 * t195;
t136 = -t108 * t36 - t35;
t132 = pkin(3) * t153 - qJDD(2) * t203 + t91;
t128 = -t102 * t155 - t145;
t127 = t19 * pkin(4) + qJDD(5) + t132;
t82 = -t168 - t190;
t126 = -pkin(6) * qJDD(3) + (t168 + t82 - t190) * qJD(3);
t125 = t136 * qJD(4) + t212;
t123 = t108 * t21 + t36 * t154 - t45 * t170 + t200 * t22;
t122 = t97 * t195 + (-qJD(1) * t52 - t142 * t97) * t110;
t121 = -t82 * qJD(2) + t216 - t58;
t120 = t125 + t210;
t119 = -g(1) * (-t106 * t97 - t98 * t179) - g(2) * (t107 * t97 - t98 * t180) - t123 + g(3) * t188;
t117 = -t59 * t52 + t120;
t116 = t59 * t50 + t119;
t115 = t177 * t50 + t119 + t186;
t101 = -qJ(5) - t206;
t94 = t200 * pkin(3) + pkin(4);
t71 = -t109 * pkin(3) - pkin(4) * t97;
t66 = pkin(2) + t72;
t49 = t50 ^ 2;
t46 = t62 * t110;
t39 = -pkin(4) * t133 - t203;
t31 = t52 * pkin(4) + t159;
t29 = t132 - t162;
t24 = qJ(5) * t133 + t191;
t23 = -t62 * qJ(5) + t147;
t20 = -t49 + t207;
t16 = -t102 * t47 - t112 * t52;
t15 = qJD(2) * t130 - t110 * t217;
t14 = -t183 + t192;
t13 = t148 + t185;
t12 = -t136 - t185;
t9 = t182 - t19;
t8 = t128 + t184;
t7 = t127 - t162;
t4 = -t46 * t100 + t16 * t102 - t112 * t19 + t50 * t172;
t3 = -t47 * t100 - t15 * t102 + t112 * t18 + t52 * t172;
t2 = -t50 * qJD(5) + t123 - t186;
t1 = -t52 * qJD(5) + t125 + t211;
t5 = [t176, 0, -t114 * t110 + t161, -qJDD(2) * t110 - t114 * t112, 0, 0, 0, 0, 0, (-0.2e1 * t153 + t163) * t112 + (-t174 * t111 - t165) * t110, (-qJDD(3) * t110 - 0.2e1 * t112 * t166) * t111 + (t174 * t110 - t161) * t109, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t15 * t50 - t16 * t52 - t46 * t18 - t47 * t19, -t1 * t46 + t10 * t16 - t7 * t112 + t12 * t15 + t26 * t172 + t2 * t47 - g(3); 0, qJDD(2), t208, -t176 * t110 + t135, t103 * qJDD(2) + 0.2e1 * t109 * t152, 0.2e1 * t109 * t163 - 0.2e1 * t175 * t166, t113 * t111 + t165, qJDD(3) * t111 - t113 * t109, 0, t126 * t109 + t209 * t111, -t209 * t109 + t126 * t111, -t18 * t62 - t52 * t27, -t133 * t18 - t62 * t19 - t217 * t52 + t27 * t50, t62 * t100 - t27 * t102, t100 * t133 - t102 * t217, 0, t147 * t100 - t203 * t19 + t59 * t217 - t29 * t133 + t138 + (t160 - t169) * t50 + t214 * t102, -t191 * t100 + t213 * t102 + t52 * t160 + t18 * t203 - t59 * t27 + t29 * t62 + t122, t23 * t100 + t201 * t102 - t133 * t7 + t144 * t50 + t39 * t19 + t217 * t26 + t138, -t24 * t100 - t202 * t102 - t39 * t18 + t25 * t52 - t26 * t27 + t7 * t62 + t122, -t1 * t62 + t10 * t27 - t12 * t217 + t133 * t2 + t23 * t18 - t24 * t19 - t201 * t52 - t202 * t50 - t216, t2 * t24 + t1 * t23 + t7 * t39 - g(3) * (-t110 * t101 + t112 * t66) + t144 * t26 + t202 * t12 + t201 * t10 + t142 * (t101 * t112 + t110 * t66); 0, 0, 0, 0, -t109 * t114 * t111, t175 * t114, t164, t163, qJDD(3), t121 * t109 - t141 * t111, t141 * t109 + t121 * t111, t194, t20, t8, t9, t100, -t148 * t102 + (t200 * t100 - t102 * t170 - t50 * t173) * pkin(3) + t117, t192 * t102 - t52 * t159 + t116 + t215, t94 * t100 - t13 * t102 - t31 * t50 - t177 * t52 + (-t35 + (-t36 - t199) * t108) * qJD(4) + t210 + t211 + t212, t14 * t102 - t31 * t52 + t115 + t215, t94 * t18 + (t12 + t13) * t52 + (-t10 + t14) * t50 + (-t108 * t19 + (t108 * t52 - t200 * t50) * qJD(4)) * pkin(3), t1 * t94 - t12 * t14 - t10 * t13 - t26 * t31 - g(1) * (t106 * t72 + t71 * t179) - g(2) * (-t107 * t72 + t71 * t180) - t71 * t196 + (t2 * t108 + (-t10 * t108 + t200 * t12) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t20, t8, t9, t100, -t136 * t102 + t117, t149 * t102 + t116, t187 + t12 * t102 + 0.2e1 * t96 + (-t151 - t26) * t52 + t120, -t207 * pkin(4) + t11 * t102 + t115, t18 * pkin(4) - t193 * t50, t193 * t12 + (-t26 * t52 + t1 + t210) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 + t182, t128 - t184, -t49 - t207, t10 * t52 + t12 * t50 + t127 - t208;];
tau_reg = t5;
