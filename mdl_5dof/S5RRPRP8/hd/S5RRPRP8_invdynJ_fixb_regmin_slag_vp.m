% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:19
% EndTime: 2021-01-15 20:52:28
% DurationCPUTime: 1.99s
% Computational Cost: add. (1819->294), mult. (3877->335), div. (0->0), fcn. (2452->6), ass. (0->168)
t138 = sin(qJ(4));
t141 = cos(qJ(4));
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t163 = t139 * t138 + t142 * t141;
t157 = t163 * qJD(4);
t195 = qJD(1) * qJD(2);
t186 = t139 * t195;
t194 = t142 * qJDD(1);
t239 = t186 - t194;
t121 = t139 * qJDD(1);
t185 = t142 * t195;
t240 = t185 + t121;
t11 = qJD(1) * t157 - t239 * t138 - t240 * t141;
t130 = qJD(2) - qJD(4);
t55 = t163 * qJD(1);
t211 = t55 * t130;
t235 = t11 + t211;
t196 = qJD(4) * t141;
t197 = qJD(4) * t138;
t198 = qJD(2) * t142;
t242 = t138 * t198 + t139 * t196 - t142 * t197;
t140 = sin(qJ(1));
t143 = cos(qJ(1));
t232 = g(1) * t143 + g(2) * t140;
t201 = qJD(1) * t139;
t117 = pkin(6) * t201;
t175 = -qJD(2) * pkin(2) + qJD(3);
t82 = t117 + t175;
t200 = qJD(1) * t142;
t118 = pkin(6) * t200;
t133 = qJD(2) * qJ(3);
t85 = t118 + t133;
t164 = t139 * t85 - t142 * t82;
t113 = pkin(6) * t194;
t131 = qJDD(2) * qJ(3);
t132 = qJD(2) * qJD(3);
t40 = -pkin(6) * t186 + t113 + t131 + t132;
t112 = pkin(6) * t121;
t188 = pkin(6) * t185 + qJDD(3) + t112;
t207 = qJDD(2) * pkin(2);
t47 = t188 - t207;
t241 = -t164 * qJD(2) + t47 * t139 + t40 * t142 - t232;
t71 = -t142 * t138 + t139 * t141;
t238 = g(3) * t163 - t232 * t71;
t12 = t242 * qJD(1) + t163 * qJDD(1) - t141 * t186;
t237 = -t12 * qJ(5) - t55 * qJD(5);
t57 = -t138 * t200 + t141 * t201;
t209 = t57 * t130;
t234 = t12 + t209;
t228 = t57 ^ 2;
t53 = t55 ^ 2;
t233 = -t53 + t228;
t144 = -pkin(2) - pkin(3);
t174 = -t138 * qJ(3) + t141 * t144;
t126 = g(1) * t140;
t225 = g(2) * t143;
t183 = t126 - t225;
t184 = -t55 * pkin(4) - qJD(5);
t60 = -qJD(1) * pkin(1) - pkin(2) * t200 - qJ(3) * t201;
t39 = pkin(3) * t200 - t60;
t20 = -t184 + t39;
t230 = t20 * t55 - t237;
t27 = -t240 * pkin(7) + t144 * qJDD(2) + t188;
t28 = t239 * pkin(7) + t40;
t190 = t144 * qJD(2);
t77 = pkin(7) * t201 - t117;
t44 = qJD(3) + t190 - t77;
t79 = -pkin(7) * t200 + t118;
t59 = t133 + t79;
t178 = t138 * t28 - t141 * t27 + t59 * t196 + t44 * t197;
t151 = -t178 + t238;
t179 = t138 * t27 + t141 * t28 + t44 * t196 - t59 * t197;
t150 = g(3) * t71 + t163 * t232 - t179;
t208 = pkin(6) * qJDD(2);
t167 = t142 * pkin(2) + t139 * qJ(3);
t83 = pkin(1) + t167;
t229 = (qJD(1) * t83 - t60) * qJD(2) + t208;
t227 = pkin(6) - pkin(7);
t182 = -t138 * t59 + t141 * t44;
t210 = t57 * qJ(5);
t7 = t182 - t210;
t6 = -t130 * pkin(4) + t7;
t226 = -t7 + t6;
t128 = qJDD(2) - qJDD(4);
t224 = t128 * pkin(4);
t222 = t39 * t55;
t221 = t39 * t57;
t220 = t57 * t55;
t165 = -t138 * t44 - t141 * t59;
t212 = t55 * qJ(5);
t8 = -t165 - t212;
t219 = t8 * t130;
t216 = t138 * t79 + t141 * t77;
t15 = t210 + t216;
t45 = t141 * qJD(3) + t174 * qJD(4);
t218 = t45 - t15;
t181 = -t138 * t77 + t141 * t79;
t14 = t181 - t212;
t81 = t141 * qJ(3) + t138 * t144;
t46 = -t138 * qJD(3) - t81 * qJD(4);
t217 = t46 - t14;
t87 = t227 * t139;
t88 = t227 * t142;
t215 = t138 * t87 + t141 * t88;
t213 = t11 * qJ(5);
t136 = qJDD(1) * pkin(1);
t146 = qJD(1) ^ 2;
t205 = t139 * t146;
t122 = t139 * qJD(3);
t204 = qJ(3) * t198 + t122;
t134 = t139 ^ 2;
t135 = t142 ^ 2;
t202 = t134 - t135;
t199 = qJD(2) * t139;
t193 = g(3) * t142 - t232 * t139;
t192 = t142 * t205;
t168 = pkin(2) * t194 + t240 * qJ(3) + qJD(1) * t122 + t136;
t171 = t139 * t190;
t13 = pkin(3) * t194 + qJD(1) * t171 + t168;
t5 = t12 * pkin(4) + qJDD(5) + t13;
t191 = -t5 + t225;
t187 = -t13 + t225;
t180 = -t138 * t88 + t141 * t87;
t102 = t141 * pkin(4) - t144;
t103 = t138 * pkin(4) + qJ(3);
t176 = t102 * t142 + t103 * t139;
t173 = t130 ^ 2;
t172 = -t112 - t193;
t145 = qJD(2) ^ 2;
t170 = pkin(6) * t145 + t225;
t166 = -t139 * pkin(2) + t142 * qJ(3);
t107 = qJ(3) * t200;
t48 = t144 * t201 + t107;
t160 = -0.2e1 * pkin(1) * t195 - t208;
t78 = t227 * t199;
t80 = qJD(2) * t88;
t159 = t138 * t80 - t141 * t78 + t87 * t196 - t88 * t197;
t67 = t142 * pkin(3) + t83;
t36 = t171 + t204;
t156 = -t170 + 0.2e1 * t136;
t154 = -t215 * qJD(4) + t138 * t78 + t141 * t80;
t25 = pkin(2) * t186 - t168;
t50 = pkin(2) * t199 - t204;
t153 = -qJD(1) * t50 + qJDD(1) * t83 - t170 - t25;
t149 = -t46 * t130 - t151;
t148 = t81 * t128 + t45 * t130 - t150;
t129 = qJ(5) - t227;
t100 = t142 * t126;
t76 = -pkin(4) + t174;
t75 = pkin(2) * t201 - t107;
t52 = t163 * t126;
t51 = t71 * t126;
t38 = pkin(1) + t176;
t31 = pkin(4) * t163 + t67;
t30 = t163 * qJD(2) - t157;
t29 = -t141 * t199 + t242;
t26 = -t57 * pkin(4) + t48;
t19 = -qJ(5) * t163 + t215;
t18 = -t71 * qJ(5) + t180;
t17 = t138 * t128 - t141 * t173 - t57 * t201;
t16 = -t141 * t128 - t138 * t173 - t55 * t201;
t9 = t29 * pkin(4) + t36;
t4 = -t30 * qJ(5) - t71 * qJD(5) + t154;
t3 = -t29 * qJ(5) - qJD(5) * t163 + t159;
t2 = t179 + t237;
t1 = -t57 * qJD(5) - t178 + t213 - t224;
t10 = [qJDD(1), t183, t232, t134 * qJDD(1) + 0.2e1 * t139 * t185, 0.2e1 * t139 * t194 - 0.2e1 * t195 * t202, qJDD(2) * t139 + t145 * t142, qJDD(2) * t142 - t145 * t139, 0, t139 * t160 + t142 * t156 + t100, t160 * t142 + (-t156 - t126) * t139, -t229 * t139 + t153 * t142 + t100, (t134 + t135) * qJDD(1) * pkin(6) + t241, t229 * t142 + (t153 + t126) * t139, t60 * t50 + (-t25 + t183) * t83 + t241 * pkin(6), -t11 * t71 + t57 * t30, t11 * t163 - t71 * t12 - t57 * t29 - t30 * t55, -t71 * t128 - t30 * t130, t128 * t163 + t29 * t130, 0, t67 * t12 - t128 * t180 - t130 * t154 - t163 * t187 + t39 * t29 + t36 * t55 + t52, -t67 * t11 + t215 * t128 + t159 * t130 - t187 * t71 + t39 * t30 + t36 * t57 + t51, t31 * t12 - t18 * t128 - t4 * t130 - t163 * t191 + t20 * t29 + t9 * t55 + t52, -t31 * t11 + t19 * t128 + t3 * t130 - t191 * t71 + t20 * t30 + t9 * t57 + t51, -t1 * t71 + t18 * t11 - t19 * t12 - t163 * t2 - t8 * t29 - t3 * t55 - t6 * t30 - t4 * t57 + t232, t2 * t19 + t8 * t3 + t1 * t18 + t6 * t4 + t5 * t31 + t20 * t9 - g(1) * (-t129 * t143 - t38 * t140) - g(2) * (-t129 * t140 + t38 * t143); 0, 0, 0, -t192, t202 * t146, t121, t194, qJDD(2), pkin(1) * t205 + t172, g(3) * t139 - t113 + (pkin(1) * t146 + t232) * t142, 0.2e1 * t207 - qJDD(3) + (-t139 * t60 + t142 * t75) * qJD(1) + t172, t166 * qJDD(1) + ((t85 - t133) * t139 + (t175 - t82) * t142) * qJD(1), t113 + 0.2e1 * t131 + 0.2e1 * t132 + (qJD(1) * t75 - g(3)) * t139 + (qJD(1) * t60 - t232) * t142, pkin(6) * qJD(1) * t164 - t47 * pkin(2) - g(3) * t167 + t40 * qJ(3) + t85 * qJD(3) - t166 * t232 - t60 * t75, -t220, -t233, t235, t234, t128, -t128 * t174 + t130 * t181 - t48 * t55 + t149 + t221, -t216 * t130 - t48 * t57 + t148 - t222, -t213 + t14 * t130 - t26 * t55 + (qJD(5) + t20) * t57 + (pkin(4) - t76) * t128 + t149, -t15 * t130 - t26 * t57 + t148 - t230, t76 * t11 - t81 * t12 + (-t8 - t217) * t57 + (t6 - t218) * t55, t2 * t81 + t1 * t76 - t20 * t26 - g(3) * t176 + t218 * t8 + t217 * t6 - t232 * (-t102 * t139 + t103 * t142); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t192, t121, -t134 * t146 - t145, -t85 * qJD(2) + t201 * t60 + t193 + t47, 0, 0, 0, 0, 0, t16, t17, t16, t17, -t234 * t138 + t235 * t141, -t20 * t201 + (t1 - t219) * t141 + (t130 * t6 + t2) * t138 + t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t233, -t235, -t234, -t128, t130 * t165 + t151 - t221, -t130 * t182 + t150 + t222, -0.2e1 * t224 + t213 - t219 + (t184 - t20) * t57 + t151, -t228 * pkin(4) - t7 * t130 + t150 + t230, t11 * pkin(4) - t226 * t55, t226 * t8 + (-t20 * t57 + t1 + t238) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 - t209, -t11 + t211, -t53 - t228, t8 * t55 + t6 * t57 + t183 + t5;];
tau_reg = t10;
