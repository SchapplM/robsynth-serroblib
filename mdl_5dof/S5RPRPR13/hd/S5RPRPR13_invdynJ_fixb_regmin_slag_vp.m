% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:55
% EndTime: 2019-12-31 18:33:00
% DurationCPUTime: 1.79s
% Computational Cost: add. (1720->271), mult. (4062->334), div. (0->0), fcn. (3037->10), ass. (0->148)
t105 = cos(qJ(5));
t100 = cos(pkin(8));
t103 = sin(qJ(3));
t186 = cos(qJ(3));
t99 = sin(pkin(8));
t74 = t103 * t100 + t186 * t99;
t201 = t74 * qJD(1);
t204 = qJD(5) + t201;
t102 = sin(qJ(5));
t206 = t102 * t204;
t172 = t103 * t99;
t150 = qJD(1) * t172;
t149 = t186 * t100;
t138 = qJD(1) * t149;
t143 = qJDD(1) * t186;
t153 = t100 * qJDD(1);
t152 = qJD(3) * t138 + t103 * t153 + t99 * t143;
t36 = qJD(3) * t150 - t152;
t34 = -qJDD(5) + t36;
t123 = -t105 * t34 - t204 * t206;
t157 = t102 * qJD(3);
t66 = -t138 + t150;
t42 = -t105 * t66 + t157;
t207 = t204 * t42;
t169 = qJDD(1) * pkin(1);
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t202 = g(1) * t104 - g(2) * t106;
t124 = t202 - qJDD(2) + t169;
t135 = g(1) * t106 + g(2) * t104;
t98 = pkin(8) + qJ(3);
t92 = sin(t98);
t93 = cos(t98);
t116 = -g(3) * t92 - t135 * t93;
t147 = qJD(3) * t186;
t159 = qJD(3) * t103;
t155 = qJD(1) * qJD(2);
t178 = pkin(6) + qJ(2);
t193 = t178 * qJDD(1) + t155;
t49 = t193 * t99;
t50 = t193 * t100;
t79 = t178 * t99;
t75 = qJD(1) * t79;
t80 = t178 * t100;
t76 = qJD(1) * t80;
t142 = t103 * t49 + t75 * t147 + t76 * t159 - t186 * t50;
t205 = t116 - t142;
t177 = -t103 * t76 - t186 * t75;
t199 = qJD(4) - t177;
t156 = t99 * qJDD(1);
t133 = -t100 * t143 + t103 * t156;
t203 = 0.2e1 * t201 * qJD(3) + t133;
t200 = -qJD(5) + t204;
t198 = qJ(2) * qJDD(1);
t88 = g(3) * t93;
t197 = -t135 * t92 + t88;
t41 = -t103 * t79 + t186 * t80;
t25 = t74 * qJD(2) + t41 * qJD(3);
t40 = t103 * t80 + t186 * t79;
t196 = -t25 * qJD(3) - t40 * qJDD(3) + t202 * t93;
t24 = (qJD(2) * t99 + qJD(3) * t80) * t103 - qJD(2) * t149 + t79 * t147;
t195 = t24 * qJD(3) - t41 * qJDD(3) - t202 * t92;
t187 = t66 * pkin(4);
t39 = -t103 * t75 + t186 * t76;
t32 = -qJD(3) * qJ(4) - t39;
t17 = -t32 - t187;
t190 = pkin(3) + pkin(7);
t194 = t190 * t34 + (t17 - t39 + t187) * t204;
t192 = t66 ^ 2;
t191 = t201 ^ 2;
t71 = t74 * qJD(3);
t37 = qJD(1) * t71 + t133;
t188 = t37 * pkin(3);
t90 = t100 * pkin(2) + pkin(1);
t126 = -t74 * qJ(4) - t90;
t73 = -t149 + t172;
t20 = t190 * t73 + t126;
t182 = t20 * t34;
t44 = t105 * qJD(3) + t102 * t66;
t181 = t44 * t66;
t180 = t66 * t42;
t179 = t201 * t66;
t175 = t100 ^ 2 + t99 ^ 2;
t158 = qJD(5) * t105;
t10 = -qJD(5) * t157 + t105 * qJDD(3) + t102 * t37 + t66 * t158;
t174 = t10 * t105;
t171 = t66 * qJ(4);
t170 = qJD(5) * t73;
t168 = qJDD(3) * pkin(3);
t167 = t104 * t102;
t166 = t104 * t105;
t165 = t106 * t102;
t164 = t106 * t105;
t163 = t39 * qJD(3);
t161 = pkin(4) * t201 + t199;
t154 = qJDD(3) * qJ(4);
t151 = t73 * t158;
t146 = t175 * qJD(1) ^ 2;
t77 = -t90 * qJDD(1) + qJDD(2);
t114 = t36 * qJ(4) + t77;
t112 = -qJD(4) * t201 + t114;
t1 = t190 * t37 + t112;
t16 = -t190 * qJD(3) + t161;
t145 = qJD(5) * t16 + t1;
t78 = -t90 * qJD(1) + qJD(2);
t115 = -qJ(4) * t201 + t78;
t13 = t190 * t66 + t115;
t122 = t103 * t50 + t76 * t147 - t75 * t159 + t186 * t49;
t119 = qJDD(4) + t122;
t5 = -t36 * pkin(4) - t190 * qJDD(3) + t119;
t144 = -qJD(5) * t13 + t5;
t140 = t102 * qJDD(3) - t105 * t37;
t139 = 0.2e1 * t175;
t137 = t93 * pkin(3) + t92 * qJ(4);
t132 = t204 * t71 - t34 * t73;
t3 = t102 * t16 + t105 * t13;
t70 = -t100 * t147 + t99 * t159;
t130 = t70 * qJ(4) - t74 * qJD(4);
t125 = t137 + t90;
t26 = t74 * pkin(4) + t40;
t8 = -qJD(3) * qJD(4) + t142 - t154;
t6 = -t37 * pkin(4) - t8;
t121 = t17 * t71 + t26 * t34 + t6 * t73;
t120 = -t105 * t204 ^ 2 + t102 * t34;
t117 = t124 + t169;
t113 = t139 * t155 - t135;
t111 = t6 + (t204 * t190 + t171) * t204 + t116;
t23 = t66 * pkin(3) + t115;
t110 = t201 * t23 + t119 + t197;
t61 = -t92 * t167 + t164;
t60 = t92 * t166 + t165;
t59 = t92 * t165 + t166;
t58 = t92 * t164 - t167;
t52 = qJD(3) * t66;
t35 = t73 * pkin(3) + t126;
t33 = pkin(3) * t201 + t171;
t31 = -qJD(3) * pkin(3) + t199;
t27 = -t73 * pkin(4) + t41;
t19 = t71 * pkin(3) + t130;
t15 = -t70 * pkin(4) + t25;
t14 = -t71 * pkin(4) - t24;
t12 = t190 * t71 + t130;
t11 = t44 * qJD(5) + t140;
t9 = t119 - t168;
t7 = t112 + t188;
t4 = t105 * t5;
t2 = -t102 * t13 + t105 * t16;
t18 = [qJDD(1), t202, t135, t117 * t100, -t117 * t99, t139 * t198 + t113, t124 * pkin(1) + (t175 * t198 + t113) * qJ(2), -t201 * t70 - t36 * t74, -t201 * t71 + t36 * t73 - t74 * t37 + t70 * t66, -t70 * qJD(3) + t74 * qJDD(3), -t71 * qJD(3) - t73 * qJDD(3), 0, -t90 * t37 + t78 * t71 + t77 * t73 + t196, t90 * t36 - t78 * t70 + t77 * t74 + t195, t201 * t25 + t24 * t66 - t31 * t70 + t32 * t71 - t36 * t40 - t37 * t41 + t73 * t8 + t74 * t9 - t135, -t19 * t66 - t23 * t71 - t35 * t37 - t7 * t73 - t196, -t19 * t201 + t23 * t70 + t35 * t36 - t7 * t74 - t195, t23 * t19 + t32 * t24 + t31 * t25 + t7 * t35 + t9 * t40 - t8 * t41 + (-g(1) * t178 - g(2) * t125) * t106 + (g(1) * t125 - g(2) * t178) * t104, t44 * t151 + (t10 * t73 + t44 * t71) * t102, (-t102 * t42 + t105 * t44) * t71 + (t174 - t102 * t11 + (-t102 * t44 - t105 * t42) * qJD(5)) * t73, t10 * t74 + t102 * t132 + t151 * t204 - t44 * t70, t105 * t132 - t11 * t74 - t170 * t206 + t42 * t70, -t204 * t70 - t34 * t74, -g(1) * t61 - g(2) * t59 + t27 * t11 + t14 * t42 - t2 * t70 + t4 * t74 + (-t1 * t74 - t12 * t204 + t182) * t102 + (t15 * t204 - t121) * t105 + ((-t102 * t26 - t105 * t20) * t204 - t3 * t74 + t17 * t102 * t73) * qJD(5), g(1) * t60 - g(2) * t58 + t27 * t10 + t14 * t44 + t3 * t70 + (-(qJD(5) * t26 + t12) * t204 + t182 - t145 * t74 + t17 * t170) * t105 + (-(-qJD(5) * t20 + t15) * t204 - t144 * t74 + t121) * t102; 0, 0, 0, -t153, t156, -t146, -qJ(2) * t146 - t124, 0, 0, 0, 0, 0, t203, (-t66 - t150) * qJD(3) + t152, -t191 - t192, -t203, t36 + t52, t188 - t32 * t66 + (-qJD(4) - t31) * t201 + t114 - t202, 0, 0, 0, 0, 0, t120 + t180, -t123 + t181; 0, 0, 0, 0, 0, 0, 0, t179, t191 - t192, (t66 - t150) * qJD(3) + t152, -t133, qJDD(3), -t201 * t78 - t122 + t163 - t197, qJD(3) * t177 + t78 * t66 - t205, pkin(3) * t36 - qJ(4) * t37 + (-t32 - t39) * t201 + (t31 - t199) * t66, t33 * t66 + t110 - t163 - 0.2e1 * t168, 0.2e1 * t154 - t23 * t66 + t33 * t201 + (0.2e1 * qJD(4) - t177) * qJD(3) + t205, -t9 * pkin(3) - g(3) * t137 - t8 * qJ(4) - t199 * t32 - t23 * t33 - t31 * t39 + t135 * (pkin(3) * t92 - qJ(4) * t93), -t206 * t44 + t174, (-t204 * t44 - t11) * t105 + (-t10 + t207) * t102, t123 + t181, t120 - t180, t204 * t66, qJ(4) * t11 + t111 * t102 + t105 * t194 + t161 * t42 + t2 * t66, qJ(4) * t10 - t194 * t102 + t105 * t111 + t161 * t44 - t3 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 + t52, qJDD(3) - t179, -qJD(3) ^ 2 - t191, t32 * qJD(3) + t110 - t168, 0, 0, 0, 0, 0, -qJD(3) * t42 + t123, -qJD(3) * t44 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t42 ^ 2 + t44 ^ 2, t10 + t207, t200 * t44 - t140, -t34, -g(1) * t58 - g(2) * t60 - t102 * t1 + t105 * t88 - t17 * t44 + t200 * t3 + t4, g(1) * t59 - g(2) * t61 + t17 * t42 + t2 * t204 - t145 * t105 + (-t144 - t88) * t102;];
tau_reg = t18;
