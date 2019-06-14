% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:01:23
% EndTime: 2019-05-05 15:01:28
% DurationCPUTime: 1.80s
% Computational Cost: add. (3920->235), mult. (7417->255), div. (0->0), fcn. (4146->6), ass. (0->167)
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t134 = cos(qJ(4));
t175 = qJD(1) * qJD(4);
t167 = t134 * t175;
t131 = sin(qJ(4));
t174 = t131 * qJDD(1);
t108 = -t167 - t174;
t102 = qJDD(5) - t108;
t177 = qJD(1) * t134;
t104 = -t133 * qJD(4) + t130 * t177;
t106 = qJD(4) * t130 + t133 * t177;
t184 = t106 * t104;
t210 = t102 + t184;
t190 = t133 * t210;
t101 = t106 ^ 2;
t117 = qJD(1) * t131 + qJD(5);
t207 = t117 ^ 2;
t214 = -t101 - t207;
t29 = t130 * t214 + t190;
t244 = pkin(8) * t29;
t243 = t131 * t29;
t157 = pkin(4) * t131 - pkin(8) * t134;
t202 = pkin(1) + qJ(3);
t147 = t157 + t202;
t196 = t130 * t210;
t27 = -t133 * t214 + t196;
t242 = t147 * t27;
t168 = t131 * t175;
t173 = t134 * qJDD(1);
t109 = -t168 + t173;
t161 = t133 * qJDD(4) - t109 * t130;
t148 = qJD(5) * t106 - t161;
t92 = t117 * t106;
t46 = t148 - t92;
t208 = t104 ^ 2;
t85 = t208 - t207;
t241 = t131 * t46 + t134 * (-t133 * t85 + t196);
t240 = t130 * t85 + t190;
t185 = t104 * t117;
t152 = -qJDD(4) * t130 - t109 * t133;
t63 = -qJD(5) * t104 - t152;
t216 = t63 - t185;
t198 = t130 * t216;
t213 = t101 - t208;
t217 = t148 + t92;
t238 = -t131 * t213 + t134 * (t133 * t217 + t198);
t201 = qJ(2) - pkin(7);
t211 = t102 - t184;
t189 = t133 * t211;
t209 = -t207 - t208;
t220 = t130 * t209 + t189;
t195 = t130 * t211;
t219 = t133 * t209 - t195;
t231 = t131 * t219 - t134 * t217;
t237 = t147 * t220 + t201 * t231;
t236 = pkin(8) * t219;
t215 = t63 + t185;
t86 = -t101 + t207;
t230 = t134 * (-t130 * t86 + t189) + t131 * t215;
t212 = t101 + t208;
t229 = pkin(4) * t212;
t228 = t133 * t86 + t195;
t227 = qJ(6) * t216;
t222 = t134 * t212;
t218 = -t130 * t217 + t133 * t216;
t137 = qJD(1) ^ 2;
t206 = 2 * qJD(3);
t205 = pkin(5) * t148;
t204 = pkin(5) * t207;
t203 = pkin(5) * t133;
t129 = t137 * pkin(7);
t132 = sin(qJ(1));
t135 = cos(qJ(1));
t165 = g(1) * t132 - t135 * g(2);
t153 = -qJDD(2) + t165;
t144 = qJ(2) * t137 + t153;
t162 = t202 * qJDD(1);
t142 = t162 + t144;
t38 = -pkin(4) * t108 - pkin(8) * t109 - t129 + (t206 + (pkin(4) * t134 + pkin(8) * t131) * qJD(4)) * qJD(1) + t142;
t136 = qJD(4) ^ 2;
t149 = t137 * t157;
t126 = qJDD(1) * qJ(2);
t154 = g(1) * t135 + g(2) * t132;
t151 = 0.2e1 * qJD(2) * qJD(1) - t154;
t145 = qJDD(3) + t151;
t83 = -t202 * t137 + t126 + t145;
t75 = -qJDD(1) * pkin(7) + t83;
t65 = t134 * g(3) - t131 * t75;
t43 = -t136 * pkin(4) + qJDD(4) * pkin(8) - t131 * t149 - t65;
t22 = t130 * t38 + t133 * t43;
t64 = g(3) * t131 + t134 * t75;
t42 = qJDD(4) * pkin(4) + pkin(8) * t136 - t134 * t149 + t64;
t139 = -pkin(5) * t92 + 0.2e1 * qJD(6) * t106 + t42;
t138 = t139 + t227;
t13 = t138 - t205;
t200 = t13 * t134;
t197 = t130 * t215;
t191 = t133 * t215;
t188 = t134 * t42;
t187 = qJ(6) * t133;
t186 = qJDD(1) * pkin(1);
t183 = t117 * t130;
t182 = t117 * t133;
t127 = t131 ^ 2;
t181 = t127 * t137;
t128 = t134 ^ 2;
t180 = t128 * t137;
t169 = t131 * t137 * t134;
t179 = t131 * (qJDD(4) + t169);
t178 = t127 + t128;
t176 = qJD(6) * t117;
t172 = qJD(1) * t206;
t171 = t104 * t182;
t170 = t131 * t184;
t166 = -qJ(6) * t130 - pkin(4);
t21 = t130 * t43 - t133 * t38;
t6 = t130 * t21 + t133 * t22;
t69 = pkin(5) * t104 - qJ(6) * t106;
t164 = t102 * qJ(6) - t104 * t69 + t22;
t114 = 0.2e1 * t176;
t160 = t114 + t164;
t84 = t106 * t183;
t159 = t134 * (t133 * t63 - t84) + t170;
t158 = t104 * t183 - t133 * t148;
t11 = t160 - t204;
t12 = -t102 * pkin(5) - qJ(6) * t207 + t106 * t69 + qJDD(6) + t21;
t156 = pkin(5) * t12 - qJ(6) * t11;
t155 = pkin(5) * t215 + qJ(6) * t46;
t5 = t130 * t22 - t133 * t21;
t31 = -t131 * t65 + t134 * t64;
t150 = qJ(6) * t210 + t164;
t146 = (-t104 * t130 - t106 * t133) * t117;
t143 = t134 * (t84 - t171) + t131 * t102;
t141 = t134 * (t130 * t148 + t171) - t170;
t140 = pkin(5) * t211 + qJ(6) * t209 - t12;
t82 = t142 + t172;
t122 = 0.2e1 * t126;
t112 = t178 * t137;
t111 = t178 * qJDD(1);
t110 = -0.2e1 * t168 + t173;
t107 = 0.2e1 * t167 + t174;
t103 = t134 * (qJDD(4) - t169);
t93 = t144 + t186;
t81 = -t179 + t134 * (-t136 - t180);
t80 = t131 * (-t136 - t181) + t103;
t74 = -t129 + t82;
t52 = (qJD(5) + t117) * t104 + t152;
t47 = (-qJD(5) + t117) * t106 + t161;
t40 = t106 * t182 + t130 * t63;
t26 = t133 * t47 + t197;
t25 = -t133 * t46 + t197;
t24 = t130 * t47 - t191;
t23 = -t130 * t46 - t191;
t18 = t134 * t52 - t243;
t16 = t134 * t216 + t243;
t15 = t131 * t26 + t222;
t14 = t131 * t25 + t222;
t10 = qJ(6) * t212 + t12;
t9 = (-t217 - t148) * pkin(5) + t138;
t8 = t139 - t205 + 0.2e1 * t227;
t7 = (-t207 + t212) * pkin(5) + t160;
t4 = t131 * t6 + t188;
t3 = t11 * t133 + t12 * t130;
t2 = t11 * t130 - t12 * t133;
t1 = t131 * t3 + t200;
t17 = [0, 0, 0, 0, 0, qJDD(1), t165, t154, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t153 - 0.2e1 * t186, t122 + t151, qJ(2) * (-pkin(1) * t137 + t126 + t151) + pkin(1) * t93, qJDD(1), 0, 0, 0, 0, 0, 0, t122 + t145, t153 + 0.2e1 * t162 + t172, qJ(2) * t83 + t202 * t82, (t109 - t168) * t134, -t107 * t134 - t110 * t131, t103 - t131 * (t136 - t180), (-t108 + t167) * t131, t134 * (-t136 + t181) - t179, 0, t202 * t107 + t131 * t74 + t201 * t80, t202 * t110 + t134 * t74 + t201 * t81, -t201 * t111 - t202 * t112 - t31, t201 * t31 + t202 * t74, t159, -t238, t230, t141, -t241, t143, -t130 * t188 - t131 * t21 + t237, -t131 * t22 - t133 * t188 + t201 * t18 - t242, -t134 * t5 + t147 * t24 + t201 * t15, t147 * t5 + t201 * t4, t159, t230, t238, t143, t241, t141, t134 * (-t130 * t9 - t187 * t217) + t131 * t140 + t237, t134 * (t10 * t133 - t130 * t7) - t131 * t155 + t147 * t23 + t201 * t14, t134 * (-pkin(5) * t198 + t133 * t8) - t131 * (pkin(5) * t214 - t150 - 0.2e1 * t176 + t204) + t242 + t201 * t16, -t131 * t156 + (-pkin(5) * t130 + t187) * t200 + t147 * t2 + t201 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t137, -t93, 0, 0, 0, 0, 0, 0, 0, -t137, -qJDD(1), -t82, 0, 0, 0, 0, 0, 0, -t107, -t110, t112, -t74, 0, 0, 0, 0, 0, 0, -t220, t27, -t24, -t5, 0, 0, 0, 0, 0, 0, -t220, -t23, -t27, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t137, t83, 0, 0, 0, 0, 0, 0, t80, t81, -t111, t31, 0, 0, 0, 0, 0, 0, t231, t18, t15, t4, 0, 0, 0, 0, 0, 0, t231, t14, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, (-t127 + t128) * t137, t173, -t169, -t174, qJDD(4), t64, t65, 0, 0, t40, t218, t228, t158, t240, t146, -pkin(4) * t217 + t133 * t42 + t236, pkin(4) * t52 - t130 * t42 - t244, pkin(8) * t26 + t229 + t6, pkin(4) * t42 + pkin(8) * t6, t40, t228, -t218, t146, -t240, t158, t133 * t9 + t166 * t217 + t236, pkin(8) * t25 + t10 * t130 + t133 * t7 + t229, t244 + t130 * t8 + (pkin(4) + t203) * t216, pkin(8) * t3 + (-t166 + t203) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t213, t215, -t184, -t46, t102, -t21, -t22, 0, 0, t184, t215, -t213, t102, t46, -t184, t140, -t155, t114 + (-t207 - t214) * pkin(5) + t150, -t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, t215, t214, t12;];
tauJ_reg  = t17;
