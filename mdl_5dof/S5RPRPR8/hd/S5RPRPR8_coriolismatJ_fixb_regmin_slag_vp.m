% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:13
% EndTime: 2019-12-31 18:22:19
% DurationCPUTime: 2.13s
% Computational Cost: add. (1614->211), mult. (3645->340), div. (0->0), fcn. (3618->8), ass. (0->178)
t142 = sin(pkin(9));
t143 = cos(pkin(9));
t147 = cos(qJ(5));
t148 = cos(qJ(3));
t235 = t148 / 0.2e1;
t176 = t147 * t235;
t209 = t147 * t142;
t145 = sin(qJ(5));
t212 = t145 * t143;
t110 = t209 + t212;
t236 = -t148 / 0.2e1;
t177 = t110 * t236;
t211 = t145 * t148;
t46 = t143 * t211 / 0.2e1 + t142 * t176 + t177;
t186 = t148 * qJD(1);
t146 = sin(qJ(3));
t208 = t147 * t143;
t213 = t145 * t142;
t242 = t208 - t213;
t91 = t242 * t146;
t84 = t91 * t186;
t87 = t91 * qJD(5);
t243 = -qJD(3) * t46 + t84 - t87;
t139 = t142 ^ 2;
t140 = t143 ^ 2;
t127 = t140 + t139;
t241 = t242 / 0.2e1;
t240 = -t110 / 0.2e1;
t239 = t110 / 0.2e1;
t136 = -pkin(4) * t143 - pkin(3);
t238 = -t136 / 0.2e1;
t234 = pkin(7) * t148;
t233 = t146 * pkin(3);
t232 = pkin(7) + qJ(4);
t134 = sin(pkin(8)) * pkin(1) + pkin(6);
t215 = t142 * t146;
t113 = t134 * t215;
t125 = -qJ(4) * t148 + t233;
t75 = t143 * t125 + t113;
t45 = t146 * pkin(4) - t143 * t234 + t75;
t231 = t145 * t45;
t214 = t143 * t146;
t76 = t142 * t125 - t134 * t214;
t64 = -t142 * t234 + t76;
t230 = t145 * t64;
t228 = t147 * t45;
t227 = t147 * t64;
t175 = pkin(4) * t142 + t134;
t101 = t175 * t146;
t102 = t175 * t148;
t135 = -cos(pkin(8)) * pkin(1) - pkin(2);
t170 = -pkin(3) * t148 - t146 * qJ(4);
t105 = t135 + t170;
t94 = t143 * t105;
t34 = -pkin(7) * t214 + t94 + (-t134 * t142 - pkin(4)) * t148;
t207 = t148 * t134;
t114 = t143 * t207;
t70 = t142 * t105 + t114;
t42 = -pkin(7) * t215 + t70;
t19 = t145 * t34 + t147 * t42;
t92 = t242 * t148;
t2 = (t227 + t231) * t148 - t19 * t146 + t102 * t91 + t101 * t92;
t226 = t2 * qJD(1);
t225 = t75 * t142;
t224 = t76 * t143;
t69 = -t142 * t207 + t94;
t9 = (t75 * t146 + t69 * t148) * t143 + (t76 * t146 + t70 * t148) * t142;
t223 = t9 * qJD(1);
t18 = t145 * t42 - t147 * t34;
t88 = t110 * t146;
t12 = -t101 * t88 - t148 * t18;
t222 = qJD(1) * t12;
t13 = -t101 * t91 - t148 * t19;
t221 = qJD(1) * t13;
t25 = (t142 * t70 + t143 * t69) * t146;
t219 = qJD(1) * t25;
t90 = t110 * t148;
t32 = t146 * t88 - t148 * t90;
t218 = qJD(1) * t32;
t33 = -t91 * t146 + t148 * t92;
t217 = qJD(1) * t33;
t216 = qJD(1) * t91;
t210 = t146 * t148;
t23 = -t69 * t146 + (t75 - 0.2e1 * t113) * t148;
t206 = t23 * qJD(1);
t24 = t76 * t148 + (-t70 + 0.2e1 * t114) * t146;
t205 = t24 * qJD(1);
t26 = -t88 * t92 - t90 * t91;
t204 = t26 * qJD(1);
t203 = t46 * qJD(1);
t155 = -t212 / 0.2e1 - t209 / 0.2e1;
t48 = (t240 + t155) * t148;
t202 = t48 * qJD(1);
t49 = -t242 * t236 + t143 * t176 - t142 * t211 / 0.2e1;
t201 = t49 * qJD(1);
t200 = t49 * qJD(3);
t50 = (t241 - t208 / 0.2e1 + t213 / 0.2e1) * t148;
t199 = t50 * qJD(1);
t198 = t50 * qJD(3);
t197 = t88 * qJD(5);
t196 = qJD(3) * t110;
t195 = qJD(3) * t136;
t194 = qJD(4) * t148;
t193 = qJD(5) * t110;
t192 = qJD(5) * t148;
t141 = t146 ^ 2;
t107 = t127 * t141;
t191 = t107 * qJD(1);
t104 = t242 * qJD(5);
t190 = t127 * qJD(3);
t130 = t148 ^ 2 - t141;
t189 = t130 * qJD(1);
t188 = t146 * qJD(1);
t187 = t146 * qJD(3);
t185 = t148 * qJD(3);
t184 = t242 * t187;
t183 = t134 * t185;
t182 = t146 * t194;
t181 = t135 * t188;
t180 = t135 * t186;
t179 = t146 * t185;
t178 = t146 * t186;
t173 = qJD(3) * t48 - t84;
t172 = t127 * t148;
t171 = qJD(4) + t195;
t169 = t224 - t225;
t1 = (t228 - t230) * t148 + t18 * t146 - t102 * t88 - t101 * t90;
t168 = t1 * qJD(1);
t11 = t134 ^ 2 * t210 + t69 * t75 + t70 * t76;
t157 = t69 * t142 / 0.2e1 - t70 * t143 / 0.2e1;
t150 = -t207 / 0.2e1 - t157;
t8 = t141 * t134 / 0.2e1 + (t224 / 0.2e1 - t225 / 0.2e1) * t146 + t150 * t148;
t167 = t11 * qJD(1) + t8 * qJD(2);
t123 = t232 * t142;
t124 = t232 * t143;
t72 = -t123 * t145 + t124 * t147;
t151 = t101 * t240 + t72 * t236 + t91 * t238;
t158 = -t230 / 0.2e1 + t228 / 0.2e1;
t3 = t151 + t158;
t47 = (t239 + t155) * t148;
t166 = qJD(1) * t3 + qJD(2) * t47;
t71 = t123 * t147 + t124 * t145;
t152 = -t101 * t242 / 0.2e1 - t88 * t238 + t71 * t235;
t159 = -t231 / 0.2e1 - t227 / 0.2e1;
t4 = t152 + t159;
t165 = qJD(1) * t4 + qJD(2) * t50;
t77 = (-0.1e1 + t127) * t210;
t164 = t8 * qJD(1) + t77 * qJD(2);
t15 = -t110 * t91 - t242 * t88;
t30 = t88 ^ 2 - t91 ^ 2;
t163 = qJD(1) * t30 + qJD(3) * t15;
t37 = -t110 ^ 2 + t242 ^ 2;
t162 = qJD(1) * t15 + qJD(3) * t37;
t161 = -qJD(1) * t88 + qJD(3) * t242;
t160 = t196 + t216;
t28 = -t88 * t239 + t91 * t241;
t156 = qJD(3) * t28 - t88 * t216;
t154 = -qJD(1) * t28 - t196 * t242;
t118 = t127 * qJ(4);
t95 = (0.1e1 / 0.2e1 - t140 / 0.2e1 - t139 / 0.2e1) * t146;
t153 = -qJD(1) * t150 + qJD(2) * t95 - qJD(3) * t118;
t149 = t170 * qJD(3) + t194;
t137 = t187 / 0.2e1;
t106 = (t186 - qJD(5) / 0.2e1) * t146;
t100 = t110 * t187;
t96 = (0.1e1 + t127) * t146 / 0.2e1;
t51 = t155 * t148 + t177;
t41 = t49 * qJD(5);
t40 = t50 * qJD(5);
t29 = -t88 * t186 + t200;
t27 = t28 * qJD(5);
t22 = t207 / 0.2e1 - t157;
t20 = -t198 + (-qJD(5) + t186) * t88;
t14 = t15 * qJD(5);
t7 = t8 * qJD(3);
t6 = -t151 + t158;
t5 = -t152 + t159;
t10 = [0, 0, 0, 0, t179, t130 * qJD(3), 0, 0, 0, t135 * t187, t135 * t185, -t23 * qJD(3) + t143 * t182, t24 * qJD(3) - t142 * t182, -qJD(3) * t9 + qJD(4) * t107, qJD(3) * t11 - qJD(4) * t25, (qJD(3) * t92 - t197) * t91, qJD(3) * t26 + qJD(5) * t30, -t33 * qJD(3) + t192 * t88, -t32 * qJD(3) + t192 * t91, -t179, -t1 * qJD(3) - t13 * qJD(5) + t194 * t91, t2 * qJD(3) + t12 * qJD(5) - t194 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t178, t189, t185, -t187, 0, t181 - t183, t134 * t187 + t180, -t114 * qJD(3) + t149 * t142 - t206, t142 * t183 + t149 * t143 + t205, t169 * qJD(3) - t223, (-pkin(3) * t207 + qJ(4) * t169) * qJD(3) + t22 * qJD(4) + t167, t160 * t92 + t27, t204 + (-t110 * t90 + t242 * t92) * qJD(3) + t14, t100 - t40 - t217, -qJD(5) * t46 + t184 - t218, -t106, (-t102 * t242 + t136 * t90 - t146 * t71) * qJD(3) - t48 * qJD(4) + t6 * qJD(5) - t168, t226 + (t102 * t110 + t136 * t92 - t146 * t72) * qJD(3) + t49 * qJD(4) + t5 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJD(3) * t142 + t143 * t188) * t148, (qJD(3) * t143 - t142 * t188) * t148, t191, qJD(3) * t22 - t219, 0, 0, 0, 0, 0, -t173, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t163, t20, t243, t137, qJD(3) * t6 - qJD(5) * t19 - t221, qJD(3) * t5 + qJD(5) * t18 + t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, -t185, -t143 * t187, t142 * t187, qJD(3) * t172, (qJ(4) * t172 - t233) * qJD(3) + t96 * qJD(4) + t164, 0, 0, 0, 0, 0, qJD(5) * t51 - t184, t100 - t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t51 - t87, t197 - t200; 0, 0, 0, 0, -t178, -t189, 0, 0, 0, -t181, -t180, t206, -t205, t223, qJD(4) * t150 - t167, -t216 * t92 + t27, t14 - t204, -t41 + t217, -qJD(5) * t48 + t218, t106, -qJD(4) * t46 - qJD(5) * t3 + t168, qJD(4) * t50 - qJD(5) * t4 - t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t95 - t164, 0, 0, 0, 0, 0, -qJD(5) * t47, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * qJD(4), t118 * qJD(4), t110 * t104, t37 * qJD(5), 0, 0, 0, t136 * t193, t136 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, -t153, 0, 0, 0, 0, 0, -t203, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t162, t104 - t201, -t193 - t202, -t188 / 0.2e1, -qJD(5) * t72 + t110 * t195 - t166, qJD(5) * t71 + t195 * t242 - t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 * t178, t142 * t178, -t191, -qJD(3) * t150 + t219, 0, 0, 0, 0, 0, -t243, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t153, 0, 0, 0, 0, 0, t193 + t203, t104 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t163, t29, t173, t137, qJD(3) * t3 - qJD(4) * t91 + t221, qJD(3) * t4 + qJD(4) * t88 - t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t47, t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t162, t201, t202, t188 / 0.2e1, -t110 * t171 + t166, -t171 * t242 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, -t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
