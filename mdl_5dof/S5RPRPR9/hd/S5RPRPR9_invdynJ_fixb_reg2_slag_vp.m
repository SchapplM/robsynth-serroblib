% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR9
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:42
% EndTime: 2019-12-31 18:24:46
% DurationCPUTime: 2.21s
% Computational Cost: add. (2174->331), mult. (4280->412), div. (0->0), fcn. (2441->10), ass. (0->185)
t114 = sin(qJ(3));
t178 = qJD(1) * qJD(3);
t163 = t114 * t178;
t117 = cos(qJ(3));
t172 = t117 * qJDD(1);
t233 = -t163 + t172;
t111 = sin(pkin(8));
t90 = pkin(1) * t111 + pkin(6);
t235 = qJD(2) * qJD(3) + qJDD(1) * t90;
t102 = t117 * qJD(2);
t71 = t90 * qJD(1);
t40 = t114 * t71 - t102;
t234 = -qJD(4) - t40;
t199 = qJ(4) * t117;
t144 = pkin(7) * t114 - t199;
t179 = t114 * qJD(4);
t127 = qJD(3) * t144 - t179;
t104 = t114 * qJ(4);
t161 = -pkin(2) - t104;
t112 = cos(pkin(8));
t217 = pkin(1) * t112;
t225 = pkin(3) + pkin(7);
t129 = -t117 * t225 + t161 - t217;
t85 = pkin(3) * t163;
t10 = qJD(1) * t127 + qJDD(1) * t129 + t85;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t193 = (pkin(4) * qJD(1) + t71) * t114 + qJD(4) - t102;
t23 = -qJD(3) * t225 + t193;
t25 = t129 * qJD(1);
t6 = t113 * t23 + t116 * t25;
t162 = t117 * t178;
t173 = t114 * qJDD(1);
t135 = t162 + t173;
t184 = qJD(3) * t117;
t158 = -t117 * qJDD(2) + t114 * t235 + t71 * t184;
t145 = qJDD(4) + t158;
t8 = pkin(4) * t135 - qJDD(3) * t225 + t145;
t2 = -qJD(5) * t6 - t113 * t10 + t116 * t8;
t190 = qJD(1) * t114;
t86 = qJD(5) + t190;
t232 = t6 * t86 + t2;
t30 = -qJD(3) * pkin(3) - t234;
t187 = qJD(3) * t113;
t189 = qJD(1) * t117;
t60 = t116 * t189 + t187;
t197 = qJD(3) * t60;
t57 = qJDD(5) + t135;
t47 = t116 * t57;
t231 = t197 - t47;
t143 = t60 * t86;
t21 = t60 * qJD(5) - t116 * qJDD(3) + t113 * t233;
t230 = t21 - t143;
t165 = t113 * t189;
t185 = qJD(3) * t116;
t62 = -t165 + t185;
t213 = t62 * t86;
t22 = -qJD(5) * t165 + qJDD(3) * t113 + (qJD(3) * qJD(5) + t233) * t116;
t229 = -t22 + t213;
t41 = qJD(2) * t114 + t117 * t71;
t32 = -qJD(3) * qJ(4) - t41;
t200 = pkin(1) * qJDD(1);
t108 = qJ(1) + pkin(8);
t97 = sin(t108);
t98 = cos(t108);
t153 = g(1) * t98 + g(2) * t97;
t95 = pkin(4) * t189;
t24 = -t32 + t95;
t228 = -t225 * t57 + t24 * t86;
t181 = qJD(5) * t117;
t166 = t116 * t181;
t186 = qJD(3) * t114;
t227 = -t113 * (-t117 * t57 + t186 * t86) + t86 * t166;
t5 = -t113 * t25 + t116 * t23;
t1 = qJD(5) * t5 + t116 * t10 + t113 * t8;
t147 = t113 * t5 - t116 * t6;
t226 = -qJD(5) * t147 + t1 * t113 + t2 * t116;
t224 = g(1) * t97;
t221 = g(2) * t98;
t220 = t5 * t86;
t218 = pkin(4) + t90;
t216 = pkin(7) * t117;
t215 = g(3) * t114;
t107 = g(3) * t117;
t105 = t117 * pkin(3);
t214 = t62 * t60;
t212 = -t114 * t21 + t184 * t62;
t167 = t113 * t181;
t168 = t114 * t185;
t211 = (t167 + t168) * t86;
t210 = t113 * t57;
t209 = t113 * t62;
t208 = t114 * t22;
t207 = t114 * t97;
t206 = t114 * t98;
t205 = t116 * t21;
t204 = t116 * t60;
t203 = t117 * t98;
t201 = t22 * t113;
t91 = -pkin(2) - t217;
t72 = qJD(1) * t91;
t198 = qJD(3) * t24;
t196 = qJDD(3) * pkin(3);
t195 = t113 * t114;
t194 = t114 * t116;
t192 = t105 + t104;
t109 = t114 ^ 2;
t110 = t117 ^ 2;
t191 = t109 - t110;
t183 = qJD(5) * t113;
t182 = qJD(5) * t116;
t70 = qJDD(1) * t91;
t180 = qJDD(3) * t90;
t176 = qJDD(1) * t109;
t175 = qJDD(1) * t110;
t174 = qJDD(3) * qJ(4);
t171 = g(1) * t206 + g(2) * t207 - t107;
t118 = cos(qJ(1));
t170 = pkin(1) * t118 + pkin(2) * t98 + pkin(6) * t97;
t169 = t60 * t186;
t52 = t218 * t117;
t115 = sin(qJ(1));
t164 = -pkin(1) * t115 + pkin(6) * t98;
t159 = -t114 * qJDD(2) - t117 * t235 + t71 * t186;
t156 = t62 * t168;
t154 = t114 * t162;
t152 = pkin(3) * t203 + t104 * t98 + t170;
t151 = g(1) * t115 - g(2) * t118;
t120 = qJD(3) ^ 2;
t150 = t120 * t90 + t221;
t148 = t113 * t6 + t116 * t5;
t48 = t91 - t192;
t39 = t48 - t216;
t51 = t218 * t114;
t17 = t113 * t51 + t116 * t39;
t16 = -t113 * t39 + t116 * t51;
t142 = t113 * t86;
t141 = t161 - t105;
t137 = -t153 + (t175 + t176) * t90;
t136 = -qJ(4) * t184 - t179;
t14 = t145 - t196;
t134 = t141 - t217;
t133 = -t150 - 0.2e1 * t70;
t132 = 0.2e1 * qJD(3) * t72 - t180;
t33 = t134 * qJD(1);
t131 = t180 + (-qJD(1) * t48 - t33) * qJD(3);
t130 = qJD(3) * t41 - t158 + t171;
t128 = qJD(3) * qJD(4) - t159 + t174;
t15 = qJD(1) * t136 + qJDD(1) * t134 + t85;
t94 = pkin(3) * t186;
t46 = t136 + t94;
t126 = qJD(1) * t46 + qJDD(1) * t48 + t15 + t150;
t125 = t14 * t114 + t128 * t117 + (t114 * t32 + t117 * t30) * qJD(3);
t124 = t158 * t114 - t159 * t117 + (-t114 * t41 + t117 * t40) * qJD(3);
t9 = pkin(4) * t233 + t128;
t123 = qJD(5) * t225 * t86 - t153 * t117 - t215 + t9;
t121 = qJD(1) ^ 2;
t96 = pkin(3) * t190;
t84 = t114 * t121 * t117;
t78 = t117 * t224;
t75 = t98 * t199;
t73 = t97 * t199;
t68 = t191 * t121;
t67 = qJDD(3) * t114 + t117 * t120;
t66 = qJDD(3) * t117 - t114 * t120;
t63 = -qJ(4) * t189 + t96;
t50 = -0.2e1 * t154 + t175;
t49 = 0.2e1 * t154 + t176;
t45 = qJD(3) * t52;
t44 = t218 * t186;
t42 = qJD(1) * t144 + t96;
t38 = t116 * t98 - t195 * t97;
t37 = t113 * t98 + t194 * t97;
t36 = t116 * t97 + t195 * t98;
t35 = -t113 * t97 + t194 * t98;
t34 = t113 * t169;
t31 = t94 + t127;
t29 = 0.2e1 * t114 * t172 - 0.2e1 * t178 * t191;
t28 = t41 + t95;
t26 = t33 * t190;
t12 = t113 * t28 + t116 * t42;
t11 = -t113 * t42 + t116 * t28;
t4 = -qJD(5) * t17 - t113 * t31 + t116 * t45;
t3 = qJD(5) * t16 + t113 * t45 + t116 * t31;
t7 = [0, 0, 0, 0, 0, qJDD(1), t151, g(1) * t118 + g(2) * t115, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t112 * t200 - t221 + t224, -0.2e1 * t111 * t200 + t153, 0, (t151 + (t111 ^ 2 + t112 ^ 2) * t200) * pkin(1), t49, t29, t67, t50, t66, 0, t114 * t132 + t117 * t133 + t78, t132 * t117 + (-t133 - t224) * t114, t124 + t137, t70 * t91 - g(1) * (-pkin(2) * t97 + t164) - g(2) * t170 + t124 * t90, 0, -t67, -t66, t49, t29, t50, t125 + t137, t114 * t131 + t117 * t126 - t78, t131 * t117 + (-t126 + t224) * t114, -g(1) * t164 - g(2) * t152 + t125 * t90 - t141 * t224 + t15 * t48 + t33 * t46, -t62 * t166 + (t117 * t21 + t186 * t62) * t113, t156 - t34 + (t201 + t205 + (t204 + t209) * qJD(5)) * t117, t212 - t227, -t60 * t167 + (t117 * t22 - t169) * t116, -t208 + (-t197 - t47) * t117 + t211, t114 * t57 + t184 * t86, -g(1) * t38 - g(2) * t36 + t16 * t57 + t22 * t52 + t4 * t86 - t44 * t60 + (-t185 * t24 + t2) * t114 + (qJD(3) * t5 + t116 * t9 - t183 * t24) * t117, g(1) * t37 - g(2) * t35 - t17 * t57 - t21 * t52 - t3 * t86 - t44 * t62 + (t187 * t24 - t1) * t114 + (-qJD(3) * t6 - t113 * t9 - t182 * t24) * t117, t16 * t21 - t17 * t22 - t3 * t60 - t4 * t62 + t78 - t147 * t186 + (qJD(5) * t148 - t1 * t116 + t113 * t2 - t221) * t117, t1 * t17 + t6 * t3 + t2 * t16 + t5 * t4 + t9 * t52 - t24 * t44 - g(1) * (pkin(4) * t98 + t164) - g(2) * (pkin(7) * t203 + t152) + (-g(1) * (t141 - t216) - g(2) * pkin(4)) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t66, -t67, 0, -t159 * t114 - t158 * t117 - g(3) + (t114 * t40 + t117 * t41) * qJD(3), 0, 0, 0, 0, 0, 0, 0, -t66, t67, t128 * t114 - t14 * t117 - g(3) + (t114 * t30 - t117 * t32) * qJD(3), 0, 0, 0, 0, 0, 0, t117 * t231 + t208 + t211, t212 + t227, -t156 - t34 + (t201 - t205 + (t204 - t209) * qJD(5)) * t117, -g(3) + (qJD(3) * t148 + t9) * t114 + (t198 - t226) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t68, t173, t84, t172, qJDD(3), -t190 * t72 + t130, t215 - qJD(3) * t40 + (-qJD(1) * t72 + t153) * t117 + t159, 0, 0, qJDD(3), -t173, -t172, -t84, t68, t84, (-pkin(3) * t114 + t199) * qJDD(1), -t189 * t63 + qJDD(4) - t130 - 0.2e1 * t196 + t26, 0.2e1 * t174 + (qJD(1) * t63 - g(3)) * t114 + (0.2e1 * qJD(4) + t40) * qJD(3) + (qJD(1) * t33 - t153) * t117 - t159, t128 * qJ(4) - t14 * pkin(3) - t33 * t63 - t30 * t41 - g(1) * (-pkin(3) * t206 + t75) - g(2) * (-pkin(3) * t207 + t73) - g(3) * t192 + t234 * t32, -t142 * t62 - t205, (-t22 - t213) * t116 + (t21 + t143) * t113, -t86 * t183 + t47 + (-t117 * t62 - t195 * t86) * qJD(1), t116 * t143 + t201, -t86 * t182 - t210 + (t117 * t60 - t194 * t86) * qJD(1), -t86 * t189, qJ(4) * t22 - t11 * t86 + t123 * t113 + t116 * t228 - t5 * t189 + t193 * t60, -qJ(4) * t21 - t113 * t228 + t123 * t116 + t12 * t86 + t6 * t189 + t193 * t62, t11 * t62 + t12 * t60 + (-t6 * t190 - t225 * t21 - t2 + (t225 * t60 - t6) * qJD(5)) * t116 + (t5 * t190 + t225 * t22 - t1 + (-t225 * t62 + t5) * qJD(5)) * t113 + t171, t9 * qJ(4) - t6 * t12 - t5 * t11 - g(1) * t75 - g(2) * t73 - g(3) * (t192 + t216) + t193 * t24 + (t114 * t153 - t226) * t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, qJDD(3) + t84, -t109 * t121 - t120, qJD(3) * t32 + t14 - t171 + t26, 0, 0, 0, 0, 0, 0, -t142 * t86 - t231, -t116 * t86 ^ 2 - qJD(3) * t62 - t210, t113 * t229 + t116 * t230, -t198 + t232 * t116 + (t1 - t220) * t113 - t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, -t60 ^ 2 + t62 ^ 2, -t230, -t214, t229, t57, -g(1) * t35 - g(2) * t37 + t107 * t116 - t24 * t62 + t232, g(1) * t36 - g(2) * t38 + t24 * t60 + t220 + (-qJD(5) * t23 - t10) * t116 + (qJD(5) * t25 - t107 - t8) * t113, 0, 0;];
tau_reg = t7;
