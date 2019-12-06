% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:59
% EndTime: 2019-12-05 17:45:09
% DurationCPUTime: 1.96s
% Computational Cost: add. (2104->182), mult. (5030->271), div. (0->0), fcn. (5497->8), ass. (0->171)
t246 = qJD(4) + qJD(5);
t135 = sin(qJ(5));
t131 = sin(pkin(9));
t138 = cos(qJ(4));
t205 = t138 * t131;
t133 = cos(pkin(9));
t136 = sin(qJ(4));
t206 = t136 * t133;
t111 = t205 + t206;
t132 = sin(pkin(8));
t92 = t111 * t132;
t222 = t135 * t92;
t137 = cos(qJ(5));
t204 = t138 * t133;
t207 = t136 * t131;
t240 = t204 - t207;
t95 = t240 * t132;
t91 = t137 * t95;
t243 = t91 - t222;
t59 = t135 * t95 + t137 * t92;
t245 = -t243 ^ 2 + t59 ^ 2;
t254 = t245 * qJD(1);
t208 = t131 * t132;
t173 = pkin(3) * t208 + qJ(2) * t132;
t75 = pkin(4) * t92 + t173;
t253 = t75 * t59;
t252 = t243 * t75;
t134 = cos(pkin(8));
t172 = qJD(1) * t134;
t20 = t59 * t172;
t251 = qJD(1) * t59;
t250 = qJD(4) * t59;
t51 = t243 * t172;
t249 = t59 * qJD(5);
t50 = -t111 * t137 - t135 * t240;
t143 = t50 * t134;
t94 = t111 * t134;
t217 = t137 * t94;
t96 = t240 * t134;
t220 = t135 * t96;
t148 = -t220 / 0.2e1 - t217 / 0.2e1;
t142 = t143 / 0.2e1 + t148;
t233 = t142 * qJD(1);
t248 = t246 * t50 - t233;
t49 = t111 * t135 - t137 * t240;
t144 = t49 * t134;
t216 = t137 * t96;
t221 = t135 * t94;
t149 = t221 / 0.2e1 - t216 / 0.2e1;
t141 = t144 / 0.2e1 + t149;
t234 = t141 * qJD(1);
t247 = t246 * t49 - t234;
t182 = t243 * qJD(1);
t244 = qJD(2) * t141 + qJD(3) * t59;
t231 = t133 ^ 2;
t232 = t131 ^ 2;
t239 = t231 + t232;
t139 = -t143 / 0.2e1 + t148;
t238 = qJD(2) * t139;
t140 = -t144 / 0.2e1 + t149;
t237 = qJD(2) * t140;
t235 = qJD(2) * t142;
t161 = t91 / 0.2e1;
t230 = pkin(4) * t95;
t228 = t134 * pkin(4);
t227 = pkin(4) * qJD(4);
t226 = pkin(4) * qJD(5);
t112 = -pkin(2) * t134 - qJ(3) * t132 - pkin(1);
t106 = t133 * t112;
t72 = -t133 * t132 * pkin(6) + t106 + (-qJ(2) * t131 - pkin(3)) * t134;
t210 = qJ(2) * t134;
t88 = t112 * t131 + t133 * t210;
t77 = -pkin(6) * t208 + t88;
t42 = t136 * t77 - t138 * t72;
t30 = -pkin(7) * t95 - t42;
t29 = t30 - t228;
t152 = t228 / 0.2e1 - t29 / 0.2e1;
t150 = t30 / 0.2e1 + t152;
t1 = t150 * t135;
t225 = t1 * qJD(1);
t224 = t135 * t30;
t43 = -t136 * t72 - t138 * t77;
t31 = -pkin(7) * t92 - t43;
t223 = t135 * t31;
t219 = t137 * t30;
t218 = t137 * t31;
t3 = t150 * t137;
t215 = t3 * qJD(1);
t11 = -t218 - t224;
t5 = t11 * t134 - t230 * t59 - t252;
t214 = t5 * qJD(1);
t12 = t219 - t223;
t6 = t12 * t134 + t230 * t243 - t253;
t213 = t6 * qJD(1);
t9 = -t137 * t29 + t223;
t7 = -t134 * t9 - t253;
t212 = t7 * qJD(1);
t10 = -t135 * t29 - t218;
t8 = t10 * t134 - t252;
t211 = t8 * qJD(1);
t16 = -t134 * t42 - t173 * t92;
t202 = t16 * qJD(1);
t17 = t134 * t43 - t173 * t95;
t201 = t17 * qJD(1);
t21 = (-t217 - t220) * t134 - t132 * t59;
t200 = t21 * qJD(1);
t22 = (t216 - t221) * t134 + t132 * t243;
t199 = t22 * qJD(1);
t45 = 0.2e1 * t161 - t222;
t193 = t45 * qJD(1);
t129 = t132 ^ 2;
t87 = -t131 * t210 + t106;
t46 = t129 * qJ(2) + (-t131 * t87 + t133 * t88) * t134;
t192 = t46 * qJD(1);
t47 = (t131 * t88 + t133 * t87) * t132;
t191 = t47 * qJD(1);
t48 = t92 ^ 2 - t95 ^ 2;
t190 = t48 * qJD(1);
t54 = t161 - t91 / 0.2e1;
t185 = t54 * qJD(1);
t55 = -t132 * t92 - t134 * t94;
t184 = t55 * qJD(1);
t56 = t132 * t95 + t134 * t96;
t183 = t56 * qJD(1);
t181 = t243 * qJD(4);
t146 = -t206 / 0.2e1 - t205 / 0.2e1;
t68 = (-t111 / 0.2e1 + t146) * t134;
t179 = t68 * qJD(1);
t145 = -t204 / 0.2e1 + t207 / 0.2e1;
t69 = (-t240 / 0.2e1 + t145) * t134;
t178 = t69 * qJD(1);
t177 = t92 * qJD(1);
t176 = t95 * qJD(1);
t175 = t95 * qJD(4);
t97 = (0.1e1 / 0.2e1 + t232 / 0.2e1 + t231 / 0.2e1) * t132;
t174 = t97 * qJD(1);
t121 = t134 ^ 2 + t129;
t171 = qJD(3) * t134;
t170 = qJD(4) * t134;
t107 = t239 * t129;
t169 = t107 * qJD(1);
t108 = t121 * t131;
t168 = t108 * qJD(1);
t109 = t121 * t133;
t167 = t109 * qJD(1);
t117 = t121 * qJ(2);
t166 = t117 * qJD(1);
t165 = t121 * qJD(1);
t164 = t59 * t182;
t163 = t243 * t251;
t162 = t92 * t176;
t160 = t132 * t172;
t159 = t132 * t171;
t79 = t95 * t172;
t158 = -t79 + t175;
t157 = qJD(5) * t54 + t51;
t156 = pkin(4) * t246;
t155 = -qJD(5) + t172;
t154 = t131 * t160;
t153 = t133 * t160;
t151 = -t249 - t250;
t147 = qJD(5) * t45 + t181 - t51;
t98 = (0.1e1 / 0.2e1 - t239 / 0.2e1) * t132;
t78 = t92 * t172;
t71 = (t111 / 0.2e1 + t146) * t134;
t70 = (t240 / 0.2e1 + t145) * t134;
t57 = -qJD(4) * t92 + t78;
t14 = -t249 - (qJD(4) - t172) * t59;
t4 = t223 - t219 / 0.2e1 + t152 * t137;
t2 = -t218 - t224 / 0.2e1 + t152 * t135;
t13 = [0, 0, 0, 0, 0, t121 * qJD(2), t117 * qJD(2), qJD(2) * t108 + t133 * t159, qJD(2) * t109 - t131 * t159, t107 * qJD(3), qJD(2) * t46 - qJD(3) * t47, -t92 * t175, t48 * qJD(4), t92 * t170, t95 * t170, 0, -qJD(2) * t55 - qJD(4) * t17 + t171 * t95, qJD(2) * t56 + qJD(4) * t16 - t171 * t92, t151 * t243, t246 * t245, -t151 * t134, (qJD(5) * t243 + t181) * t134, 0, -qJD(2) * t21 - qJD(4) * t5 - qJD(5) * t8 + t171 * t243, qJD(2) * t22 + qJD(4) * t6 + qJD(5) * t7 - t171 * t59; 0, 0, 0, 0, 0, t165, t166, t168, t167, 0, qJD(3) * t98 + t192, 0, 0, 0, 0, 0, qJD(4) * t71 - t184, qJD(4) * t70 + t183, 0, 0, 0, 0, 0, t139 * t246 - t200, t140 * t246 + t199; 0, 0, 0, 0, 0, 0, 0, t153, -t154, t169, qJD(2) * t98 - t191, 0, 0, 0, 0, 0, t79, -t78, 0, 0, 0, 0, 0, t157, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t190, t57, -t158, 0, qJD(2) * t71 + qJD(4) * t43 - t201, qJD(2) * t70 + qJD(4) * t42 + t202, -t163, t254, t14, -t147, 0, qJD(4) * t11 + qJD(5) * t2 - t214 + t238, -qJD(4) * t12 + qJD(5) * t4 + t213 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t254, t155 * t59 - t250, -t45 * qJD(4) + t155 * t243, 0, qJD(3) * t54 + qJD(4) * t2 + qJD(5) * t10 - t211 + t238, qJD(4) * t4 + qJD(5) * t9 + t212 + t237; 0, 0, 0, 0, 0, -t165, -t166, -t168, -t167, 0, -qJD(3) * t97 - t192, 0, 0, 0, 0, 0, -qJD(4) * t68 + t184, -qJD(4) * t69 - t183, 0, 0, 0, 0, 0, -t142 * t246 + t200, -t141 * t246 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t111 - t179, -qJD(4) * t240 - t178, 0, 0, 0, 0, 0, t248, t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t247; 0, 0, 0, 0, 0, 0, 0, -t153, t154, -t169, qJD(2) * t97 + t191, 0, 0, 0, 0, 0, t158, t57, 0, 0, 0, 0, 0, t147, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, -t177, 0, 0, 0, 0, 0, t182, -t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, -t190, -t78, -t79, 0, qJD(2) * t68 - qJD(3) * t95 + t201, qJD(2) * t69 + qJD(3) * t92 - t202, t163, -t254, -t20, -t157, 0, -qJD(3) * t243 + qJD(5) * t1 + t214 + t235, qJD(5) * t3 - t213 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t178, 0, 0, 0, 0, 0, t233, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t177, 0, 0, 0, 0, 0, -t182, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135 * t226, -t137 * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, 0, -t135 * t156 + t225, -t137 * t156 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t254, -t20, qJD(4) * t54 - t51, 0, -qJD(3) * t45 - qJD(4) * t1 + t211 + t235, -qJD(4) * t3 - t212 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, 0, t135 * t227 - t225, t137 * t227 - t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;
