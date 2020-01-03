% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:59
% EndTime: 2019-12-31 19:55:04
% DurationCPUTime: 1.82s
% Computational Cost: add. (3240->285), mult. (7596->355), div. (0->0), fcn. (5543->12), ass. (0->166)
t149 = sin(qJ(4));
t146 = sin(pkin(8));
t147 = cos(pkin(8));
t150 = sin(qJ(2));
t152 = cos(qJ(2));
t101 = t146 * t152 + t147 * t150;
t174 = t101 * qJD(2);
t166 = qJD(1) * t174;
t182 = t146 * t150 - t147 * t152;
t172 = t182 * qJDD(1);
t156 = -t172 - t166;
t226 = cos(qJ(4));
t93 = t182 * qJD(1);
t94 = t101 * qJD(1);
t178 = t149 * t93 - t226 * t94;
t200 = qJD(1) * qJD(2);
t192 = t152 * t200;
t193 = t150 * t200;
t62 = qJDD(1) * t101 - t146 * t193 + t147 * t192;
t160 = qJD(4) * t178 - t149 * t62 + t226 * t156;
t142 = qJD(2) + qJD(4);
t213 = t178 * t142;
t242 = t160 - t213;
t194 = qJD(4) * t226;
t201 = qJD(4) * t149;
t176 = t149 * t156 - t93 * t194 - t94 * t201 + t226 * t62;
t52 = -t149 * t94 - t226 * t93;
t212 = t52 * t142;
t7 = t176 - t212;
t220 = t52 ^ 2;
t240 = t178 ^ 2;
t241 = -t220 + t240;
t138 = t152 * pkin(2);
t131 = t138 + pkin(1);
t110 = -qJD(1) * t131 + qJD(3);
t68 = t93 * pkin(3) + t110;
t20 = -pkin(4) * t52 + qJ(5) * t178 + t68;
t239 = t20 * t52;
t238 = t68 * t52;
t219 = t178 * t52;
t237 = t178 * t68;
t143 = qJ(2) + pkin(8);
t137 = qJ(4) + t143;
t128 = sin(t137);
t129 = cos(t137);
t205 = t129 * pkin(4) + t128 * qJ(5);
t236 = pkin(3) * cos(t143) + t138 + t205;
t140 = qJDD(2) + qJDD(4);
t136 = t140 * pkin(4);
t231 = qJDD(5) - t136;
t235 = -t178 * t20 + t231;
t27 = -pkin(4) * t178 - t52 * qJ(5);
t130 = t147 * pkin(2) + pkin(3);
t223 = t146 * pkin(2);
t196 = t149 * t223;
t228 = t93 * pkin(7);
t218 = qJ(3) + pkin(6);
t118 = t218 * t150;
t106 = qJD(1) * t118;
t119 = t218 * t152;
t107 = qJD(1) * t119;
t207 = t147 * t107;
t66 = t146 * t106 - t207;
t41 = t66 + t228;
t227 = t94 * pkin(7);
t96 = t146 * t107;
t67 = -t147 * t106 - t96;
t42 = t67 - t227;
t234 = qJD(4) * t196 - t130 * t194 + t149 * t41 + t226 * t42;
t204 = t149 * t130 + t226 * t223;
t132 = t140 * qJ(5);
t135 = t142 * qJD(5);
t233 = t132 + t135;
t75 = pkin(3) * t182 - t131;
t151 = sin(qJ(1));
t153 = cos(qJ(1));
t232 = g(1) * t151 - g(2) * t153;
t69 = -t147 * t118 - t146 * t119;
t45 = -t101 * pkin(7) + t69;
t70 = -t146 * t118 + t147 * t119;
t46 = -pkin(7) * t182 + t70;
t26 = t149 * t45 + t226 * t46;
t179 = -t149 * t46 + t226 * t45;
t173 = t182 * qJD(2);
t188 = qJD(2) * t218;
t90 = t152 * qJD(3) - t150 * t188;
t91 = -t150 * qJD(3) - t152 * t188;
t43 = -t146 * t90 + t147 * t91;
t34 = pkin(7) * t173 + t43;
t44 = t146 * t91 + t147 * t90;
t35 = -pkin(7) * t174 + t44;
t4 = t179 * qJD(4) + t149 * t34 + t226 * t35;
t230 = t128 * t232 + t26 * t140 + t4 * t142;
t224 = g(3) * t152;
t222 = t150 * pkin(2);
t59 = qJDD(2) * pkin(2) + t91 * qJD(1) - qJDD(1) * t118;
t65 = t90 * qJD(1) + qJDD(1) * t119;
t29 = t146 * t59 + t147 * t65;
t217 = qJD(5) - t234;
t216 = t204 * qJD(4) - t149 * t42 + t226 * t41;
t215 = qJD(2) * pkin(2);
t100 = -t106 + t215;
t61 = t146 * t100 + t207;
t60 = t147 * t100 - t96;
t38 = qJD(2) * pkin(3) - t227 + t60;
t40 = t61 - t228;
t18 = t149 * t38 + t226 * t40;
t214 = t18 * t142;
t211 = t128 * t151;
t210 = t128 * t153;
t209 = t129 * t151;
t208 = t129 * t153;
t17 = -t149 * t40 + t226 * t38;
t206 = qJD(5) - t17;
t144 = t150 ^ 2;
t202 = -t152 ^ 2 + t144;
t199 = t150 * qJDD(1);
t198 = t152 * qJDD(1);
t197 = pkin(2) * t193 + qJDD(3);
t134 = t150 * t215;
t71 = t94 * pkin(3) + qJD(1) * t222;
t191 = -pkin(4) * t128 - pkin(3) * sin(t143) - t222;
t28 = -t146 * t65 + t147 * t59;
t16 = qJDD(2) * pkin(3) - t62 * pkin(7) + t28;
t19 = pkin(7) * t156 + t29;
t187 = t149 * t16 + t226 * t19 + t38 * t194 - t40 * t201;
t186 = t149 * t19 - t226 * t16 + t40 * t194 + t38 * t201;
t185 = g(1) * t153 + g(2) * t151;
t180 = pkin(1) + t236;
t177 = -0.2e1 * pkin(1) * t200 - pkin(6) * qJDD(2);
t175 = g(1) * t208 + g(2) * t209 + g(3) * t128 - t187;
t171 = t226 * t130 - t196;
t169 = -qJDD(1) * t131 + t197;
t168 = g(1) * t210 + g(2) * t211 - g(3) * t129 - t186;
t167 = t226 * t182;
t165 = t176 + t212;
t164 = t17 * t142 + t175;
t5 = t26 * qJD(4) + t149 * t35 - t226 * t34;
t163 = g(1) * t209 - g(2) * t208 + t140 * t179 - t5 * t142;
t72 = pkin(3) * t174 + t134;
t154 = qJD(2) ^ 2;
t162 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t154 + t232;
t155 = qJD(1) ^ 2;
t161 = pkin(1) * t155 - pkin(6) * qJDD(1) + t185;
t64 = t226 * t101 - t149 * t182;
t159 = -t216 * t142 + t168;
t158 = -t168 + t235;
t157 = -t160 - t213;
t39 = pkin(3) * t166 + t75 * qJDD(1) + t197;
t3 = -pkin(4) * t160 - qJ(5) * t176 + qJD(5) * t178 + t39;
t141 = -pkin(7) - t218;
t112 = qJ(5) * t208;
t111 = qJ(5) * t209;
t88 = -pkin(4) - t171;
t87 = qJ(5) + t204;
t63 = t149 * t101 + t167;
t31 = t64 * qJD(4) - t149 * t173 + t226 * t174;
t30 = t101 * t201 + t142 * t167 + t149 * t174;
t24 = t63 * pkin(4) - t64 * qJ(5) + t75;
t23 = t27 + t71;
t9 = t142 * qJ(5) + t18;
t8 = -t142 * pkin(4) + t206;
t6 = t31 * pkin(4) + t30 * qJ(5) - t64 * qJD(5) + t72;
t2 = t186 + t231;
t1 = t187 + t233;
t10 = [qJDD(1), t232, t185, t144 * qJDD(1) + 0.2e1 * t150 * t192, 0.2e1 * t150 * t198 - 0.2e1 * t202 * t200, qJDD(2) * t150 + t154 * t152, qJDD(2) * t152 - t154 * t150, 0, t150 * t177 + t152 * t162, -t150 * t162 + t152 * t177, -t44 * t93 - t29 * t182 - t43 * t94 - t69 * t62 - t28 * t101 - t70 * t172 + (-t101 * t61 + t182 * t60 - t70 * t94) * qJD(2) - t185, t29 * t70 + t61 * t44 + t28 * t69 + t60 * t43 - t169 * t131 + t110 * t134 - g(1) * (-t151 * t131 + t153 * t218) - g(2) * (t153 * t131 + t151 * t218), t176 * t64 + t178 * t30, t160 * t64 - t176 * t63 + t178 * t31 - t30 * t52, t64 * t140 - t30 * t142, -t63 * t140 - t31 * t142, 0, -t160 * t75 + t68 * t31 + t39 * t63 - t52 * t72 + t163, t176 * t75 - t178 * t72 - t68 * t30 + t39 * t64 - t230, -t160 * t24 + t20 * t31 + t3 * t63 - t52 * t6 + t163, -t1 * t63 + t160 * t26 - t176 * t179 - t178 * t5 + t2 * t64 - t8 * t30 - t9 * t31 + t4 * t52 - t185, -t176 * t24 + t178 * t6 + t20 * t30 - t3 * t64 + t230, t1 * t26 - t2 * t179 + t20 * t6 + t3 * t24 + t9 * t4 + t8 * t5 + (g(1) * t141 - g(2) * t180) * t153 + (g(1) * t180 + g(2) * t141) * t151; 0, 0, 0, -t150 * t155 * t152, t202 * t155, t199, t198, qJDD(2), t150 * t161 - t224, g(3) * t150 + t152 * t161, (t61 + t66) * t94 - (-t67 + t60) * t93 + (-t147 * t62 + ((-t192 - t199) * t146 + (-t193 + t198) * t147) * t146) * pkin(2), -t60 * t66 - t61 * t67 + (-t224 + t29 * t146 + t147 * t28 + (-qJD(1) * t110 + t185) * t150) * pkin(2), t219, t241, t7, t242, t140, t140 * t171 + t52 * t71 + t159 + t237, -t204 * t140 + t234 * t142 + t178 * t71 + t175 - t238, -t88 * t140 + t23 * t52 + t159 - t235, t87 * t160 + t176 * t88 + (t217 - t8) * t52 + (-t216 - t9) * t178, t87 * t140 + t217 * t142 - t178 * t23 - t175 + t233 + t239, t1 * t87 + t2 * t88 - t20 * t23 - g(1) * (t153 * t191 + t112) - g(2) * (t151 * t191 + t111) - g(3) * t236 + t217 * t9 + t216 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 ^ 2 - t94 ^ 2, t60 * t94 + t61 * t93 + t169 - t232, 0, 0, 0, 0, 0, t157, t165, t157, -t220 - t240, -t165, t178 * t8 - t9 * t52 - t232 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, t241, t7, t242, t140, t168 + t214 + t237, t164 - t238, t27 * t52 + t136 - t158 + t214, -pkin(4) * t176 + t160 * qJ(5) - (-t18 + t9) * t178 - (t8 - t206) * t52, -t178 * t27 + 0.2e1 * t132 + 0.2e1 * t135 - t164 + t239, t1 * qJ(5) - t2 * pkin(4) - t20 * t27 - t8 * t18 - g(1) * (-pkin(4) * t210 + t112) - g(2) * (-pkin(4) * t211 + t111) - g(3) * t205 + t206 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 + t219, t7, -t142 ^ 2 - t240, -t9 * t142 + t158;];
tau_reg = t10;
