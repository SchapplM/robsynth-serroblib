% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:30
% EndTime: 2019-12-31 20:58:35
% DurationCPUTime: 1.75s
% Computational Cost: add. (2421->319), mult. (5322->352), div. (0->0), fcn. (3513->8), ass. (0->178)
t146 = sin(qJ(3));
t147 = sin(qJ(2));
t149 = cos(qJ(2));
t237 = cos(qJ(3));
t81 = t146 * t149 + t237 * t147;
t70 = t81 * qJD(1);
t135 = t149 * pkin(2);
t250 = t135 + pkin(1);
t90 = t250 * qJD(1);
t252 = -t70 * qJ(4) - t90;
t140 = qJDD(2) + qJDD(3);
t141 = qJD(2) + qJD(3);
t245 = t140 * qJ(4) + t141 * qJD(4);
t204 = t149 * qJDD(1);
t206 = qJD(1) * qJD(2);
t251 = t147 * t206 - t204;
t239 = pkin(7) + pkin(6);
t190 = qJDD(1) * t237;
t205 = t147 * qJDD(1);
t177 = t146 * t205 - t149 * t190;
t47 = t141 * t81;
t27 = qJD(1) * t47 + t177;
t197 = t237 * t149;
t184 = qJD(1) * t197;
t208 = qJD(1) * t147;
t195 = t146 * t208;
t68 = -t184 + t195;
t249 = t27 * qJ(5) + t68 * qJD(5);
t222 = t70 * t141;
t67 = t70 ^ 2;
t248 = -t141 ^ 2 - t67;
t92 = t239 * t149;
t86 = qJD(1) * t92;
t228 = t146 * t86;
t230 = qJD(2) * pkin(2);
t91 = t239 * t147;
t84 = qJD(1) * t91;
t78 = -t84 + t230;
t40 = t237 * t78 - t228;
t247 = qJD(4) - t40;
t246 = 0.2e1 * t245;
t144 = qJ(2) + qJ(3);
t133 = sin(t144);
t134 = cos(t144);
t211 = t134 * pkin(3) + t133 * qJ(4);
t148 = sin(qJ(1));
t150 = cos(qJ(1));
t191 = g(1) * t148 - g(2) * t150;
t138 = g(1) * t150;
t210 = g(2) * t148 + t138;
t131 = t140 * pkin(3);
t244 = qJDD(4) - t131;
t215 = t146 * t147;
t178 = t141 * t215;
t187 = -t141 * t184 - t146 * t204 - t147 * t190;
t26 = qJD(1) * t178 + t187;
t243 = -t140 * pkin(4) + t26 * qJ(5) - t70 * qJD(5);
t194 = t237 * qJD(3);
t107 = pkin(2) * t194 + qJD(4);
t119 = t146 * pkin(2) + qJ(4);
t242 = t107 * t141 + t119 * t140 + t245;
t241 = t68 ^ 2;
t240 = pkin(3) + pkin(4);
t238 = t70 * pkin(4);
t236 = pkin(2) * t147;
t120 = t134 * pkin(4);
t32 = t68 * pkin(3) + t252;
t234 = t32 * t68;
t233 = t70 * t68;
t232 = t90 * t68;
t38 = t70 * pkin(3) + t68 * qJ(4);
t45 = -t237 * t84 - t228;
t77 = t237 * t86;
t41 = t146 * t78 + t77;
t53 = -t146 * t91 + t237 * t92;
t231 = qJ(4) * t27;
t229 = t119 * t27;
t227 = t41 * t141;
t44 = -t146 * t84 + t77;
t226 = t44 * t141;
t225 = t68 * qJ(5);
t224 = t68 * t141;
t62 = t70 * qJ(5);
t30 = t62 + t45;
t221 = t107 - t30;
t220 = t107 - t45;
t145 = qJDD(1) * pkin(1);
t219 = t133 * t148;
t218 = t133 * t150;
t217 = t134 * t148;
t216 = t134 * t150;
t214 = qJ(5) - t239;
t24 = t62 + t40;
t213 = qJD(4) - t24;
t142 = t147 ^ 2;
t209 = -t149 ^ 2 + t142;
t207 = qJD(3) * t146;
t203 = t147 * t230;
t202 = pkin(2) * t207;
t64 = t251 * pkin(2) - t145;
t200 = t135 + t211;
t198 = qJD(2) * t239;
t196 = t141 * t207;
t192 = t149 * t206;
t50 = qJDD(2) * pkin(2) + t239 * (-t192 - t205);
t51 = t239 * t251;
t189 = t146 * t50 + t78 * t194 - t86 * t207 - t237 * t51;
t188 = -t146 * t51 + t86 * t194 + t78 * t207 - t237 * t50;
t186 = g(2) * (pkin(3) * t216 + qJ(4) * t218 + t150 * t250);
t125 = -t237 * pkin(2) - pkin(3);
t185 = t81 * qJ(4) + t250;
t25 = t225 + t41;
t29 = t44 + t225;
t183 = -t29 + t202;
t182 = -t44 + t202;
t181 = g(1) * t219 - g(2) * t218;
t180 = g(1) * t217 - g(2) * t216;
t179 = -pkin(3) * t133 - t236;
t176 = -t140 + t233;
t34 = pkin(2) * t208 + t38;
t6 = t189 + t245;
t4 = t27 * pkin(3) + t26 * qJ(4) - t70 * qJD(4) + t64;
t7 = t188 + t244;
t175 = -t250 - t211;
t52 = t146 * t92 + t237 * t91;
t174 = -0.2e1 * pkin(1) * t206 - pkin(6) * qJDD(2);
t85 = t147 * t198;
t87 = t149 * t198;
t15 = -t146 * t87 - t91 * t194 - t92 * t207 - t237 * t85;
t173 = g(1) * t218 + g(2) * t219 - g(3) * t134 - t188;
t172 = -g(1) * t216 - g(2) * t217 - g(3) * t133 + t189;
t46 = -qJD(2) * t197 - t149 * t194 + t178;
t169 = -t46 * qJ(4) + t81 * qJD(4) - t203;
t168 = t173 - t244;
t167 = t90 * t70 + t173;
t166 = t40 * t141 - t172;
t165 = t45 * t141 - t172;
t1 = -t27 * pkin(4) + qJDD(5) - t4;
t164 = t53 * t140 + t15 * t141 + t181;
t16 = t53 * qJD(3) - t146 * t85 + t237 * t87;
t163 = -t52 * t140 - t16 * t141 + t180;
t153 = qJD(2) ^ 2;
t162 = -pkin(6) * t153 + 0.2e1 * t145 + t191;
t154 = qJD(1) ^ 2;
t161 = pkin(1) * t154 - pkin(6) * qJDD(1) + t210;
t160 = t210 * t240 * t133;
t159 = -t141 * t195 - t187;
t14 = -t240 * t68 + qJD(5) - t252;
t158 = t14 * t68 + t172 + t249;
t157 = -t32 * t70 + t168;
t156 = -t14 * t70 - t168 + t243;
t155 = -t177 - t222;
t132 = t141 * qJ(4);
t115 = -pkin(4) + t125;
t96 = qJ(4) * t216;
t94 = qJ(4) * t217;
t93 = pkin(2) * t196;
t80 = -t197 + t215;
t39 = t80 * pkin(3) - t185;
t37 = t132 + t41;
t36 = t80 * qJ(5) + t53;
t35 = -t81 * qJ(5) + t52;
t33 = -t141 * pkin(3) + t247;
t31 = t67 - t241;
t28 = -t240 * t80 + t185;
t23 = -t38 - t238;
t18 = t132 + t25;
t17 = -t34 - t238;
t13 = -t240 * t141 + t247 - t62;
t12 = t155 + t222;
t11 = t159 + t224;
t10 = t47 * pkin(3) - t169;
t9 = t46 * qJ(5) - t81 * qJD(5) + t16;
t8 = t47 * qJ(5) + t80 * qJD(5) + t15;
t5 = -t240 * t47 + t169;
t3 = t6 + t249;
t2 = t7 + t243;
t19 = [qJDD(1), t191, t210, t142 * qJDD(1) + 0.2e1 * t147 * t192, 0.2e1 * t147 * t204 - 0.2e1 * t209 * t206, qJDD(2) * t147 + t153 * t149, qJDD(2) * t149 - t153 * t147, 0, t147 * t174 + t149 * t162, -t147 * t162 + t149 * t174, -t26 * t81 - t70 * t46, t26 * t80 - t81 * t27 + t46 * t68 - t70 * t47, t81 * t140 - t46 * t141, -t80 * t140 - t47 * t141, 0, t68 * t203 - t250 * t27 - t90 * t47 + t64 * t80 + t163, t70 * t203 + t250 * t26 + t90 * t46 + t64 * t81 - t164, t10 * t68 + t39 * t27 + t32 * t47 + t4 * t80 + t163, -t15 * t68 + t16 * t70 - t52 * t26 - t53 * t27 - t33 * t46 - t37 * t47 - t6 * t80 + t7 * t81 - t210, -t10 * t70 + t39 * t26 + t32 * t46 - t4 * t81 + t164, t6 * t53 + t37 * t15 + t4 * t39 + t32 * t10 + t7 * t52 + t33 * t16 - t239 * t138 - t186 + (-g(1) * t175 - g(2) * t239) * t148, -t1 * t80 - t14 * t47 - t35 * t140 - t9 * t141 - t28 * t27 - t5 * t68 + t180, t1 * t81 - t14 * t46 + t36 * t140 + t8 * t141 - t28 * t26 + t5 * t70 + t181, t13 * t46 + t18 * t47 - t2 * t81 + t35 * t26 + t36 * t27 + t3 * t80 + t8 * t68 - t9 * t70 + t210, t3 * t36 + t18 * t8 + t2 * t35 + t13 * t9 + t1 * t28 + t14 * t5 - t186 + (g(1) * t214 - g(2) * t120) * t150 + (-g(1) * (t175 - t120) + g(2) * t214) * t148; 0, 0, 0, -t147 * t154 * t149, t209 * t154, t205, t204, qJDD(2), -g(3) * t149 + t147 * t161, g(3) * t147 + t149 * t161, t233, t31, t11, t12, t140, t226 + (t237 * t140 - t68 * t208 - t196) * pkin(2) + t167, -t232 + (-t140 * t146 - t141 * t194 - t70 * t208) * pkin(2) + t165, -t125 * t140 - t34 * t68 + t157 + t226 - t93, -t229 - t125 * t26 + (t182 + t37) * t70 + (t33 - t220) * t68, t34 * t70 - t165 - t234 + t242, t6 * t119 + t7 * t125 - t32 * t34 - g(1) * (t150 * t179 + t96) - g(2) * (t148 * t179 + t94) - g(3) * t200 + t220 * t37 + t182 * t33, -t115 * t140 + t29 * t141 + t17 * t68 - t156 - t93, -t30 * t141 - t17 * t70 + t158 + t242, t115 * t26 + t229 + (-t18 - t183) * t70 + (-t13 + t221) * t68, t3 * t119 + t2 * t115 - t14 * t17 - g(1) * (-t150 * t236 + t96) - g(2) * (-t148 * t236 + t94) - g(3) * (t120 + t200) + t221 * t18 + t160 + t183 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t31, t11, t12, t140, t167 + t227, t166 - t232, -t38 * t68 + t131 + t157 + t227, pkin(3) * t26 - t231 + (t37 - t41) * t70 + (t33 - t247) * t68, t38 * t70 - t166 - t234 + t246, t6 * qJ(4) - t7 * pkin(3) - t32 * t38 - t33 * t41 - g(1) * (-pkin(3) * t218 + t96) - g(2) * (-pkin(3) * t219 + t94) - g(3) * t211 + t247 * t37, t140 * t240 + t25 * t141 + t23 * t68 - t156, -t24 * t141 - t23 * t70 + t158 + t246, t231 - t240 * t26 + (-t18 + t25) * t70 + (-t13 + t213) * t68, t3 * qJ(4) - t2 * t240 - t13 * t25 - t14 * t23 - g(1) * t96 - g(2) * t94 - g(3) * (t120 + t211) + t213 * t18 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t11, t248, -t37 * t141 - t157, t176, t248, t26 - t224, -t18 * t141 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155 - t222, t159 - t224, -t67 - t241, t13 * t70 - t18 * t68 + t1 + t191;];
tau_reg = t19;
