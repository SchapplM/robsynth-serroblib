% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:18
% EndTime: 2019-12-31 18:30:24
% DurationCPUTime: 2.48s
% Computational Cost: add. (2966->313), mult. (7159->417), div. (0->0), fcn. (5586->14), ass. (0->164)
t140 = pkin(8) + qJ(3);
t133 = sin(t140);
t135 = cos(t140);
t148 = sin(qJ(1));
t150 = cos(qJ(1));
t179 = g(1) * t150 + g(2) * t148;
t159 = -g(3) * t135 + t179 * t133;
t147 = sin(qJ(3));
t221 = cos(qJ(3));
t142 = sin(pkin(8));
t216 = pkin(6) + qJ(2);
t117 = t216 * t142;
t112 = qJD(1) * t117;
t144 = cos(pkin(8));
t119 = t216 * t144;
t113 = qJD(1) * t119;
t76 = -t147 * t112 + t221 * t113;
t230 = t76 * qJD(3);
t193 = qJD(1) * qJD(2);
t222 = t216 * qJDD(1) + t193;
t92 = t222 * t142;
t93 = t222 * t144;
t166 = t147 * t93 + t221 * t92 + t230;
t30 = -qJDD(3) * pkin(3) + qJDD(4) + t166;
t154 = t30 - t159;
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t111 = t221 * t142 + t147 * t144;
t103 = t111 * qJD(1);
t141 = sin(pkin(9));
t143 = cos(pkin(9));
t82 = -t143 * qJD(3) + t141 * t103;
t84 = t141 * qJD(3) + t143 * t103;
t37 = t146 * t84 + t149 * t82;
t187 = t221 * t144;
t124 = qJD(1) * t187;
t202 = t147 * t142;
t186 = qJD(1) * t202;
t101 = -t124 + t186;
t96 = qJD(5) + t101;
t229 = t37 * t96;
t171 = t146 * t82 - t149 * t84;
t228 = t171 * t96;
t110 = t149 * t141 + t146 * t143;
t105 = t110 * qJD(5);
t208 = t110 * t101 + t105;
t207 = qJDD(1) * pkin(1);
t225 = g(1) * t148 - g(2) * t150;
t170 = -qJDD(2) + t207 + t225;
t75 = -t221 * t112 - t147 * t113;
t227 = t225 * t133;
t167 = t225 * t135;
t226 = -qJD(5) + t96;
t224 = qJ(2) * qJDD(1);
t108 = t146 * t141 - t149 * t143;
t209 = t96 * t108;
t107 = t111 * qJD(3);
t183 = qJDD(1) * t221;
t192 = t142 * qJDD(1);
t175 = -t144 * t183 + t147 * t192;
t73 = qJD(1) * t107 + t175;
t70 = qJDD(5) + t73;
t223 = -t110 * t70 + t209 * t96;
t218 = g(3) * t133;
t157 = -t179 * t135 - t218;
t191 = t144 * qJDD(1);
t190 = qJD(3) * t124 + t142 * t183 + t147 * t191;
t72 = -qJD(3) * t186 + t190;
t57 = -t143 * qJDD(3) + t141 * t72;
t58 = t141 * qJDD(3) + t143 * t72;
t9 = -t171 * qJD(5) + t146 * t58 + t149 * t57;
t97 = t101 ^ 2;
t220 = pkin(7) * t143;
t215 = pkin(7) + qJ(4);
t127 = t144 * pkin(2) + pkin(1);
t114 = -t127 * qJDD(1) + qJDD(2);
t21 = t73 * pkin(3) - t72 * qJ(4) - t103 * qJD(4) + t114;
t168 = -t147 * t92 + t221 * t93;
t27 = qJDD(3) * qJ(4) + (qJD(4) + t75) * qJD(3) + t168;
t7 = t141 * t21 + t143 * t27;
t163 = t187 - t202;
t106 = t163 * qJD(3);
t43 = t107 * pkin(3) - t106 * qJ(4) - t111 * qJD(4);
t164 = -t221 * t117 - t147 * t119;
t52 = qJD(2) * t163 + qJD(3) * t164;
t19 = t141 * t43 + t143 * t52;
t115 = -t127 * qJD(1) + qJD(2);
t51 = t101 * pkin(3) - t103 * qJ(4) + t115;
t68 = qJD(3) * qJ(4) + t76;
t25 = t141 * t51 + t143 * t68;
t69 = t103 * pkin(3) + t101 * qJ(4);
t32 = t141 * t69 + t143 * t75;
t71 = -pkin(3) * t163 - t111 * qJ(4) - t127;
t81 = -t147 * t117 + t221 * t119;
t34 = t141 * t71 + t143 * t81;
t214 = t103 * t37;
t212 = t141 * t73;
t211 = t143 * t73;
t210 = t171 * t103;
t206 = t101 * t141;
t205 = t106 * t141;
t204 = t111 * t141;
t203 = t111 * t143;
t139 = pkin(9) + qJ(5);
t132 = sin(t139);
t201 = t148 * t132;
t134 = cos(t139);
t200 = t148 * t134;
t199 = t150 * t132;
t198 = t150 * t134;
t66 = -qJD(3) * pkin(3) + qJD(4) - t75;
t197 = -qJD(4) + t66;
t196 = t142 ^ 2 + t144 ^ 2;
t195 = qJD(5) * t146;
t194 = qJD(5) * t149;
t6 = -t141 * t27 + t143 * t21;
t2 = t73 * pkin(4) - t58 * pkin(7) + t6;
t5 = -t57 * pkin(7) + t7;
t189 = -t146 * t5 + t149 * t2;
t18 = -t141 * t52 + t143 * t43;
t24 = -t141 * t68 + t143 * t51;
t31 = -t141 * t75 + t143 * t69;
t33 = -t141 * t81 + t143 * t71;
t182 = t196 * qJD(1) ^ 2;
t181 = -t108 * t70 - t208 * t96;
t180 = 0.2e1 * t196;
t177 = -t7 * t141 - t6 * t143;
t176 = t146 * t2 + t149 * t5;
t11 = t101 * pkin(4) - t84 * pkin(7) + t24;
t14 = -t82 * pkin(7) + t25;
t3 = t149 * t11 - t146 * t14;
t4 = t146 * t11 + t149 * t14;
t174 = -t24 * t141 + t25 * t143;
t20 = -pkin(4) * t163 - pkin(7) * t203 + t33;
t26 = -pkin(7) * t204 + t34;
t173 = -t146 * t26 + t149 * t20;
t172 = t146 * t20 + t149 * t26;
t169 = t135 * pkin(3) + t133 * qJ(4) + t127;
t8 = -t146 * t57 + t149 * t58 - t82 * t194 - t84 * t195;
t116 = t215 * t141;
t162 = pkin(7) * t206 - t143 * qJD(4) + qJD(5) * t116 + t32;
t118 = t215 * t143;
t161 = t103 * pkin(4) + t141 * qJD(4) + qJD(5) * t118 + t101 * t220 + t31;
t160 = t170 + t207;
t156 = t66 * t106 + t30 * t111 - t179;
t153 = t180 * t193 - t179;
t53 = qJD(2) * t111 + qJD(3) * t81;
t128 = -t143 * pkin(4) - pkin(3);
t90 = t135 * t198 + t201;
t89 = -t135 * t199 + t200;
t88 = -t135 * t200 + t199;
t87 = t135 * t201 + t198;
t63 = t108 * t111;
t62 = t110 * t111;
t54 = pkin(4) * t204 - t164;
t42 = -pkin(4) * t206 + t76;
t36 = t82 * pkin(4) + t66;
t35 = pkin(4) * t205 + t53;
t29 = t110 * t106 + t194 * t203 - t195 * t204;
t28 = -t111 * t105 - t108 * t106;
t13 = -pkin(7) * t205 + t19;
t12 = t57 * pkin(4) + t30;
t10 = t107 * pkin(4) - t106 * t220 + t18;
t1 = [qJDD(1), t225, t179, t160 * t144, -t160 * t142, t180 * t224 + t153, t170 * pkin(1) + (t196 * t224 + t153) * qJ(2), t103 * t106 + t72 * t111, -t106 * t101 - t103 * t107 - t111 * t73 + t163 * t72, t106 * qJD(3) + t111 * qJDD(3), -t107 * qJD(3) + qJDD(3) * t163, 0, -t53 * qJD(3) + qJDD(3) * t164 + t115 * t107 - t114 * t163 - t127 * t73 + t167, -t52 * qJD(3) - t81 * qJDD(3) + t115 * t106 + t114 * t111 - t127 * t72 - t227, t18 * t101 + t24 * t107 + t141 * t156 + t143 * t167 - t163 * t6 - t164 * t57 + t33 * t73 + t53 * t82, -t19 * t101 - t25 * t107 - t141 * t167 + t156 * t143 + t163 * t7 - t164 * t58 - t34 * t73 + t53 * t84, -t18 * t84 - t19 * t82 - t33 * t58 - t34 * t57 + t227 + t177 * t111 + (-t141 * t25 - t143 * t24) * t106, t24 * t18 + t25 * t19 - t30 * t164 + t6 * t33 + t7 * t34 + t66 * t53 + (-g(1) * t216 - g(2) * t169) * t150 + (g(1) * t169 - g(2) * t216) * t148, -t171 * t28 - t8 * t63, t171 * t29 - t28 * t37 - t8 * t62 + t63 * t9, -t107 * t171 - t163 * t8 + t28 * t96 - t63 * t70, -t37 * t107 + t163 * t9 - t29 * t96 - t62 * t70, t96 * t107 - t163 * t70, (t149 * t10 - t146 * t13) * t96 + t173 * t70 - t189 * t163 + t3 * t107 + t35 * t37 + t54 * t9 + t12 * t62 + t36 * t29 - g(1) * t88 - g(2) * t90 + (t163 * t4 - t172 * t96) * qJD(5), -(t146 * t10 + t149 * t13) * t96 - t172 * t70 + t176 * t163 - t4 * t107 - t35 * t171 + t54 * t8 - t12 * t63 + t36 * t28 - g(1) * t87 - g(2) * t89 + (t163 * t3 - t173 * t96) * qJD(5); 0, 0, 0, -t191, t192, -t182, -qJ(2) * t182 - t170, 0, 0, 0, 0, 0, 0.2e1 * t103 * qJD(3) + t175, (-t101 - t186) * qJD(3) + t190, -t103 * t82 - t141 * t97 + t211, -t103 * t84 - t143 * t97 - t212, -t141 * t57 - t143 * t58 + (t141 * t84 - t143 * t82) * t101, t174 * t101 - t66 * t103 - t177 - t225, 0, 0, 0, 0, 0, t181 - t214, t210 + t223; 0, 0, 0, 0, 0, 0, 0, t103 * t101, t103 ^ 2 - t97, (t101 - t186) * qJD(3) + t190, -t175, qJDD(3), -t115 * t103 + t159 - t166 + t230, t115 * t101 - t157 - t168, -qJ(4) * t212 - pkin(3) * t57 - t24 * t103 - t76 * t82 + (t197 * t141 - t31) * t101 - t154 * t143, -qJ(4) * t211 - pkin(3) * t58 + t25 * t103 - t76 * t84 + (t197 * t143 + t32) * t101 + t154 * t141, t31 * t84 + t32 * t82 + (-qJ(4) * t57 - qJD(4) * t82 - t101 * t24 + t7) * t143 + (qJ(4) * t58 + qJD(4) * t84 - t101 * t25 - t6) * t141 + t157, -t24 * t31 - t25 * t32 - t66 * t76 + t174 * qJD(4) - t154 * pkin(3) + (-t6 * t141 + t7 * t143 + t157) * qJ(4), t8 * t110 + t171 * t209, -t8 * t108 - t110 * t9 + t171 * t208 + t209 * t37, t210 - t223, t181 + t214, -t96 * t103, (-t149 * t116 - t146 * t118) * t70 + t128 * t9 + t12 * t108 - t3 * t103 - t42 * t37 + (t146 * t162 - t149 * t161) * t96 + t208 * t36 + t159 * t134, -(-t146 * t116 + t149 * t118) * t70 + t128 * t8 + t12 * t110 + t4 * t103 + t42 * t171 + (t146 * t161 + t149 * t162) * t96 - t209 * t36 - t159 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t84 + t57, -t101 * t82 + t58, -t82 ^ 2 - t84 ^ 2, t24 * t84 + t25 * t82 + t154, 0, 0, 0, 0, 0, t9 - t228, t8 - t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171 * t37, t171 ^ 2 - t37 ^ 2, t8 + t229, -t9 - t228, t70, -g(1) * t89 + g(2) * t87 + t132 * t218 + t171 * t36 + t226 * t4 + t189, g(1) * t90 - g(2) * t88 + t134 * t218 + t226 * t3 + t36 * t37 - t176;];
tau_reg = t1;
