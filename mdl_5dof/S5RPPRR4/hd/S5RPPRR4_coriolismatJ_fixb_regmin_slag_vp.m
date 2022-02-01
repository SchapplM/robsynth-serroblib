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
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:17:00
% EndTime: 2022-01-23 09:17:04
% DurationCPUTime: 1.75s
% Computational Cost: add. (2101->180), mult. (5009->267), div. (0->0), fcn. (5485->8), ass. (0->167)
t244 = qJD(4) + qJD(5);
t134 = sin(qJ(5));
t130 = sin(pkin(9));
t137 = cos(qJ(4));
t203 = t137 * t130;
t132 = cos(pkin(9));
t135 = sin(qJ(4));
t204 = t135 * t132;
t110 = t203 + t204;
t131 = sin(pkin(8));
t92 = t110 * t131;
t220 = t134 * t92;
t136 = cos(qJ(5));
t202 = t137 * t132;
t205 = t135 * t130;
t237 = t202 - t205;
t95 = t237 * t131;
t91 = t136 * t95;
t240 = t91 - t220;
t59 = t134 * t95 + t136 * t92;
t243 = -t240 ^ 2 + t59 ^ 2;
t252 = t243 * qJD(1);
t206 = t130 * t131;
t171 = pkin(3) * t206 + qJ(2) * t131;
t75 = pkin(4) * t92 + t171;
t251 = t75 * t59;
t250 = t240 * t75;
t133 = cos(pkin(8));
t170 = qJD(1) * t133;
t20 = t59 * t170;
t249 = qJD(1) * t59;
t248 = qJD(4) * t59;
t51 = t240 * t170;
t247 = t59 * qJD(5);
t50 = -t110 * t136 - t134 * t237;
t142 = t50 * t133;
t94 = t110 * t133;
t215 = t136 * t94;
t96 = t237 * t133;
t218 = t134 * t96;
t147 = -t218 / 0.2e1 - t215 / 0.2e1;
t141 = t142 / 0.2e1 + t147;
t231 = t141 * qJD(1);
t246 = t244 * t50 - t231;
t49 = t110 * t134 - t136 * t237;
t143 = t49 * t133;
t214 = t136 * t96;
t219 = t134 * t94;
t148 = t219 / 0.2e1 - t214 / 0.2e1;
t140 = t143 / 0.2e1 + t148;
t232 = t140 * qJD(1);
t245 = t244 * t49 - t232;
t180 = t240 * qJD(1);
t242 = -t132 ^ 2 / 0.2e1 - t130 ^ 2 / 0.2e1;
t241 = qJD(2) * t140 + qJD(3) * t59;
t138 = -t142 / 0.2e1 + t147;
t236 = qJD(2) * t138;
t139 = -t143 / 0.2e1 + t148;
t235 = qJD(2) * t139;
t233 = qJD(2) * t141;
t160 = t91 / 0.2e1;
t228 = pkin(4) * t95;
t226 = t133 * pkin(4);
t225 = pkin(4) * qJD(4);
t224 = pkin(4) * qJD(5);
t111 = -pkin(2) * t133 - qJ(3) * t131 - pkin(1);
t106 = t132 * t111;
t72 = -t132 * t131 * pkin(6) + t106 + (-qJ(2) * t130 - pkin(3)) * t133;
t208 = qJ(2) * t133;
t88 = t111 * t130 + t132 * t208;
t77 = -pkin(6) * t206 + t88;
t42 = t135 * t77 - t137 * t72;
t30 = -pkin(7) * t95 - t42;
t29 = t30 - t226;
t151 = t226 / 0.2e1 - t29 / 0.2e1;
t149 = t30 / 0.2e1 + t151;
t1 = t149 * t134;
t223 = t1 * qJD(1);
t222 = t134 * t30;
t43 = -t135 * t72 - t137 * t77;
t31 = -pkin(7) * t92 - t43;
t221 = t134 * t31;
t217 = t136 * t30;
t216 = t136 * t31;
t3 = t149 * t136;
t213 = t3 * qJD(1);
t11 = -t216 - t222;
t5 = t11 * t133 - t228 * t59 - t250;
t212 = t5 * qJD(1);
t12 = t217 - t221;
t6 = t12 * t133 + t228 * t240 - t251;
t211 = t6 * qJD(1);
t9 = -t136 * t29 + t221;
t7 = -t133 * t9 - t251;
t210 = t7 * qJD(1);
t10 = -t134 * t29 - t216;
t8 = t10 * t133 - t250;
t209 = t8 * qJD(1);
t16 = -t133 * t42 - t171 * t92;
t200 = t16 * qJD(1);
t17 = t133 * t43 - t171 * t95;
t199 = t17 * qJD(1);
t21 = (-t215 - t218) * t133 - t131 * t59;
t198 = t21 * qJD(1);
t22 = (t214 - t219) * t133 + t131 * t240;
t197 = t22 * qJD(1);
t45 = 0.2e1 * t160 - t220;
t191 = t45 * qJD(1);
t128 = t131 ^ 2;
t87 = -t130 * t208 + t106;
t46 = t128 * qJ(2) + (-t130 * t87 + t132 * t88) * t133;
t190 = t46 * qJD(1);
t47 = (t130 * t88 + t132 * t87) * t131;
t189 = t47 * qJD(1);
t48 = t92 ^ 2 - t95 ^ 2;
t188 = t48 * qJD(1);
t54 = t160 - t91 / 0.2e1;
t183 = t54 * qJD(1);
t55 = -t131 * t92 - t133 * t94;
t182 = t55 * qJD(1);
t56 = t131 * t95 + t133 * t96;
t181 = t56 * qJD(1);
t179 = t240 * qJD(4);
t145 = -t204 / 0.2e1 - t203 / 0.2e1;
t68 = (-t110 / 0.2e1 + t145) * t133;
t177 = t68 * qJD(1);
t144 = -t202 / 0.2e1 + t205 / 0.2e1;
t69 = (-t237 / 0.2e1 + t144) * t133;
t176 = t69 * qJD(1);
t175 = t92 * qJD(1);
t174 = t95 * qJD(1);
t173 = t95 * qJD(4);
t97 = (0.1e1 / 0.2e1 - t242) * t131;
t172 = t97 * qJD(1);
t120 = t133 ^ 2 + t128;
t169 = qJD(3) * t133;
t168 = qJD(4) * t133;
t107 = t120 * t130;
t167 = t107 * qJD(1);
t108 = t120 * t132;
t166 = t108 * qJD(1);
t116 = t120 * qJ(2);
t165 = t116 * qJD(1);
t164 = t120 * qJD(1);
t163 = t59 * t180;
t162 = t240 * t249;
t161 = t92 * t174;
t159 = t131 * t170;
t158 = t131 * t169;
t79 = t95 * t170;
t157 = -t79 + t173;
t156 = qJD(5) * t54 + t51;
t155 = pkin(4) * t244;
t154 = -qJD(5) + t170;
t153 = t130 * t159;
t152 = t132 * t159;
t150 = -t247 - t248;
t146 = qJD(5) * t45 + t179 - t51;
t98 = (0.1e1 / 0.2e1 + t242) * t131;
t78 = t92 * t170;
t71 = (t110 / 0.2e1 + t145) * t133;
t70 = (t237 / 0.2e1 + t144) * t133;
t57 = -qJD(4) * t92 + t78;
t14 = -t247 - (qJD(4) - t170) * t59;
t4 = t221 - t217 / 0.2e1 + t151 * t136;
t2 = -t216 - t222 / 0.2e1 + t151 * t134;
t13 = [0, 0, 0, 0, t120 * qJD(2), t116 * qJD(2), qJD(2) * t107 + t132 * t158, qJD(2) * t108 - t130 * t158, qJD(2) * t46 - qJD(3) * t47, -t92 * t173, t48 * qJD(4), t92 * t168, t95 * t168, 0, -qJD(2) * t55 - qJD(4) * t17 + t169 * t95, qJD(2) * t56 + qJD(4) * t16 - t169 * t92, t150 * t240, t244 * t243, -t150 * t133, (qJD(5) * t240 + t179) * t133, 0, -qJD(2) * t21 - qJD(4) * t5 - qJD(5) * t8 + t169 * t240, qJD(2) * t22 + qJD(4) * t6 + qJD(5) * t7 - t169 * t59; 0, 0, 0, 0, t164, t165, t167, t166, qJD(3) * t98 + t190, 0, 0, 0, 0, 0, qJD(4) * t71 - t182, qJD(4) * t70 + t181, 0, 0, 0, 0, 0, t138 * t244 - t198, t139 * t244 + t197; 0, 0, 0, 0, 0, 0, t152, -t153, qJD(2) * t98 - t189, 0, 0, 0, 0, 0, t79, -t78, 0, 0, 0, 0, 0, t156, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t188, t57, -t157, 0, qJD(2) * t71 + qJD(4) * t43 - t199, qJD(2) * t70 + qJD(4) * t42 + t200, -t162, t252, t14, -t146, 0, qJD(4) * t11 + qJD(5) * t2 - t212 + t236, -qJD(4) * t12 + qJD(5) * t4 + t211 + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t252, t154 * t59 - t248, -t45 * qJD(4) + t154 * t240, 0, qJD(3) * t54 + qJD(4) * t2 + qJD(5) * t10 - t209 + t236, qJD(4) * t4 + qJD(5) * t9 + t210 + t235; 0, 0, 0, 0, -t164, -t165, -t167, -t166, -qJD(3) * t97 - t190, 0, 0, 0, 0, 0, -qJD(4) * t68 + t182, -qJD(4) * t69 - t181, 0, 0, 0, 0, 0, -t141 * t244 + t198, -t140 * t244 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t172, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t110 - t177, -qJD(4) * t237 - t176, 0, 0, 0, 0, 0, t246, t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, t245; 0, 0, 0, 0, 0, 0, -t152, t153, qJD(2) * t97 + t189, 0, 0, 0, 0, 0, t157, t57, 0, 0, 0, 0, 0, t146, t14; 0, 0, 0, 0, 0, 0, 0, 0, t172, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, -t175, 0, 0, 0, 0, 0, t180, -t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, -t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t188, -t78, -t79, 0, qJD(2) * t68 - qJD(3) * t95 + t199, qJD(2) * t69 + qJD(3) * t92 - t200, t162, -t252, -t20, -t156, 0, -qJD(3) * t240 + qJD(5) * t1 + t212 + t233, qJD(5) * t3 - t211 + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t176, 0, 0, 0, 0, 0, t231, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t175, 0, 0, 0, 0, 0, -t180, t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t224, -t136 * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, 0, -t134 * t155 + t223, -t136 * t155 + t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t252, -t20, qJD(4) * t54 - t51, 0, -qJD(3) * t45 - qJD(4) * t1 + t209 + t233, -qJD(4) * t3 - t210 + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, 0, t134 * t225 - t223, t136 * t225 - t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;
