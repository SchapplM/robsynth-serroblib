% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:17:37
% EndTime: 2019-05-04 20:17:44
% DurationCPUTime: 2.98s
% Computational Cost: add. (7310->268), mult. (13771->382), div. (0->0), fcn. (10496->14), ass. (0->177)
t148 = cos(qJ(4));
t197 = qJD(3) * qJD(4);
t189 = t148 * t197;
t145 = sin(qJ(4));
t195 = t145 * qJDD(3);
t101 = 0.2e1 * t189 + t195;
t134 = sin(pkin(12));
t136 = sin(pkin(7));
t137 = sin(pkin(6));
t138 = cos(pkin(12));
t140 = cos(pkin(7));
t141 = cos(pkin(6));
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t129 = t145 ^ 2;
t151 = qJD(3) ^ 2;
t125 = t129 * t151;
t150 = qJD(4) ^ 2;
t113 = -t125 - t150;
t202 = t148 * t151;
t192 = t145 * t202;
t110 = qJDD(4) - t192;
t203 = t148 * t110;
t80 = t113 * t145 + t203;
t174 = t101 * t149 + t146 * t80;
t76 = t110 * t145 - t113 * t148;
t45 = t136 * t174 + t140 * t76;
t247 = (t138 * (-t136 * t76 + t140 * t174) + t134 * (-t101 * t146 + t149 * t80)) * t137 + t141 * t45;
t190 = t145 * t197;
t194 = t148 * qJDD(3);
t103 = -0.2e1 * t190 + t194;
t130 = t148 ^ 2;
t126 = t130 * t151;
t115 = -t126 - t150;
t109 = qJDD(4) + t192;
t210 = t109 * t145;
t79 = -t115 * t148 + t210;
t171 = -t103 * t149 + t146 * t79;
t74 = t109 * t148 + t115 * t145;
t44 = t136 * t171 - t140 * t74;
t245 = (t138 * (t136 * t74 + t140 * t171) + t134 * (t103 * t146 + t149 * t79)) * t137 + t141 * t44;
t240 = pkin(9) * t80;
t237 = pkin(9) * t79;
t144 = sin(qJ(6));
t147 = cos(qJ(6));
t199 = qJD(3) * t148;
t96 = qJD(4) * t144 + t147 * t199;
t98 = qJD(4) * t147 - t144 * t199;
t71 = t98 * t96;
t102 = t189 + t195;
t91 = qJDD(6) + t102;
t227 = -t71 + t91;
t232 = t144 * t227;
t231 = t147 * t227;
t230 = t102 + t189;
t135 = sin(pkin(11));
t139 = cos(pkin(11));
t108 = g(1) * t135 - g(2) * t139;
t131 = -g(3) + qJDD(1);
t169 = t141 * t108 + t131 * t137;
t184 = -g(1) * t139 - g(2) * t135;
t153 = -t134 * t184 + t138 * t169;
t162 = -t108 * t137 + t131 * t141 + qJDD(2);
t229 = t136 * t162 + t140 * t153;
t204 = t145 * qJ(5);
t182 = -pkin(4) * t148 - t204;
t59 = t134 * t169 + t138 * t184;
t37 = t146 * t229 + t149 * t59;
t35 = -t151 * pkin(3) + qJDD(3) * pkin(9) + t37;
t188 = t151 * t182 + t35;
t226 = -pkin(4) * t150 + t148 * t188;
t225 = t203 + (t126 - t150) * t145;
t90 = t96 ^ 2;
t224 = t98 ^ 2;
t200 = qJD(3) * t145;
t119 = qJD(6) + t200;
t116 = t119 ^ 2;
t223 = 2 * qJD(5);
t222 = pkin(4) + pkin(10);
t221 = t119 * t96;
t161 = t190 - t194;
t196 = qJDD(4) * qJ(5);
t198 = pkin(5) * t200 - qJD(4) * pkin(10) + t223;
t50 = -t136 * t153 + t140 * t162;
t217 = t145 * t50;
t19 = -pkin(5) * t161 - pkin(10) * t126 + qJD(4) * t198 + t196 + t217 + t226;
t220 = t144 * t19;
t62 = t71 + t91;
t219 = t144 * t62;
t187 = qJDD(4) * t144 - t147 * t161;
t159 = (-qJD(6) + t119) * t98 - t187;
t64 = -qJD(6) * t96 + qJDD(4) * t147 + t144 * t161;
t166 = t64 + t221;
t38 = t144 * t159 - t147 * t166;
t218 = t145 * t38;
t216 = t147 * t19;
t215 = t147 * t62;
t49 = t148 * t50;
t207 = t119 * t144;
t206 = t119 * t147;
t205 = t138 * t140;
t105 = (t129 + t130) * qJDD(3);
t106 = t125 + t126;
t201 = pkin(3) * t106 + pkin(9) * t105;
t193 = t145 * t71;
t191 = -t116 - t224;
t25 = t145 * t35 - t49;
t26 = t148 * t35 + t217;
t15 = t145 * t25 + t148 * t26;
t186 = -qJDD(4) * pkin(4) - qJ(5) * t150 + qJDD(5) - t49;
t185 = t146 * t59 - t149 * t229;
t20 = -qJDD(4) * pkin(10) + (t102 - t189) * pkin(5) + (-pkin(10) * t202 + t188) * t145 + t186;
t34 = -qJDD(3) * pkin(3) - pkin(9) * t151 + t185;
t157 = t161 * pkin(4) - qJ(5) * t230 + t34;
t24 = -pkin(5) * t126 - pkin(10) * t194 + (qJD(4) * t222 - t198) * t200 + t157;
t10 = t144 * t20 + t147 * t24;
t9 = t144 * t24 - t147 * t20;
t6 = t10 * t144 - t147 * t9;
t3 = t145 * t6 + t148 * t19;
t7 = t10 * t147 + t144 * t9;
t183 = t146 * t3 - t149 * t7;
t156 = qJD(4) * t223 + t226;
t154 = t156 + t196;
t21 = t154 + t217;
t22 = t188 * t145 + t186;
t13 = t145 * t22 + t148 * t21;
t27 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t200 + t157;
t181 = t13 * t146 - t149 * t27;
t180 = t146 * t15 - t149 * t34;
t60 = -t90 - t224;
t29 = t148 * t60 + t218;
t39 = t144 * t166 + t147 * t159;
t179 = t146 * t29 - t149 * t39;
t65 = -t116 - t90;
t41 = t144 * t65 + t231;
t52 = (qJD(6) + t119) * t98 + t187;
t31 = t145 * t41 + t148 * t52;
t42 = t147 * t65 - t232;
t178 = t146 * t31 - t149 * t42;
t167 = t64 - t221;
t47 = t147 * t191 - t219;
t33 = t145 * t47 + t148 * t167;
t48 = -t144 * t191 - t215;
t177 = t146 * t33 - t149 * t48;
t176 = t146 * t37 - t149 * t185;
t170 = t105 * t146 + t106 * t149;
t168 = (-t125 + t150) * t148 + t210;
t165 = qJDD(3) * t149 - t146 * t151;
t164 = -qJDD(3) * t146 - t149 * t151;
t163 = pkin(3) - t182;
t160 = -t148 * t222 - pkin(3) - t204;
t107 = t125 - t126;
t85 = t165 * t136;
t84 = t164 * t136;
t83 = t116 - t224;
t82 = t90 - t116;
t73 = t230 * t145;
t72 = t103 * t148;
t70 = -t90 + t224;
t68 = t148 * t101 + t103 * t145;
t66 = t170 * t136;
t63 = -qJD(6) * t98 - t187;
t40 = t141 * t66 + (t134 * (t105 * t149 - t106 * t146) + t170 * t205) * t137;
t32 = t145 * t167 - t148 * t47;
t30 = t145 * t52 - t148 * t41;
t28 = t145 * t60 - t148 * t38;
t18 = t136 * t176 + t140 * t50;
t17 = t136 * t177 + t140 * t32;
t16 = t136 * t178 + t140 * t30;
t14 = t145 * t26 - t148 * t25;
t12 = t145 * t21 - t148 * t22;
t11 = t136 * t179 + t140 * t28;
t5 = t136 * t180 + t14 * t140;
t4 = t12 * t140 + t136 * t181;
t2 = t145 * t19 - t148 * t6;
t1 = t136 * t183 + t140 * t2;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t162 + (t134 * t59 + t138 * t153) * t137, 0, 0, 0, 0, 0, 0, t141 * t85 + (t134 * t164 + t165 * t205) * t137, t141 * t84 + (-t134 * t165 + t164 * t205) * t137, 0, t141 * t18 + (t134 * (t146 * t185 + t149 * t37) + t138 * (-t136 * t50 + t140 * t176)) * t137, 0, 0, 0, 0, 0, 0, -t245, -t247, t40, t141 * t5 + (t134 * (t146 * t34 + t149 * t15) + t138 * (-t136 * t14 + t140 * t180)) * t137, 0, 0, 0, 0, 0, 0, t40, t245, t247, t141 * t4 + (t134 * (t13 * t149 + t146 * t27) + t138 * (-t12 * t136 + t140 * t181)) * t137, 0, 0, 0, 0, 0, 0, t141 * t16 + (t134 * (t146 * t42 + t149 * t31) + t138 * (-t136 * t30 + t140 * t178)) * t137, t141 * t17 + (t134 * (t146 * t48 + t149 * t33) + t138 * (-t136 * t32 + t140 * t177)) * t137, t141 * t11 + (t134 * (t146 * t39 + t149 * t29) + t138 * (-t136 * t28 + t140 * t179)) * t137, t141 * t1 + (t134 * (t146 * t7 + t149 * t3) + t138 * (-t136 * t2 + t140 * t183)) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, 0, 0, 0, 0, 0, 0, t85, t84, 0, t18, 0, 0, 0, 0, 0, 0, -t44, -t45, t66, t5, 0, 0, 0, 0, 0, 0, t66, t44, t45, t4, 0, 0, 0, 0, 0, 0, t16, t17, t11, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t185, -t37, 0, 0, t73, t68, t168, t72, t225, 0, pkin(3) * t103 - t148 * t34 - t237, -pkin(3) * t101 + t145 * t34 - t240, t15 + t201, -pkin(3) * t34 + pkin(9) * t15, 0, -t168, -t225, t73, t68, t72, (pkin(4) * t106 + t154) * t148 + (qJ(5) * t106 + t22 + t49) * t145 + t201, -t103 * t163 + t148 * t27 + t237, t145 * (-pkin(4) * t190 + t200 * t223 - t157) + t240 + t163 * t101, pkin(9) * t13 - t163 * t27, t193 + t148 * (-t144 * t64 - t206 * t98), t145 * t70 + t148 * (t144 * t52 - t147 * t167), t145 * t166 + t148 * (-t147 * t83 - t232), -t193 + t148 * (-t147 * t63 - t207 * t96), t145 * t159 + t148 * (-t144 * t82 - t215), t145 * t91 + t148 * (t144 * t96 + t147 * t98) * t119, t145 * (pkin(5) * t41 - t9) + t148 * (pkin(5) * t52 + t216) + pkin(9) * t31 + t160 * t42, t145 * (pkin(5) * t47 - t10) + t148 * (pkin(5) * t167 - t220) + pkin(9) * t33 + t160 * t48, pkin(5) * t218 + t148 * (pkin(5) * t60 - t7) + pkin(9) * t29 + t160 * t39, t160 * t7 + (pkin(5) + pkin(9)) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, t107, t195, t192, t194, qJDD(4), -t25, -t26, 0, 0, qJDD(4), -t195, -t194, -t192, t107, t192, (-pkin(4) * t145 + qJ(5) * t148) * qJDD(3), -pkin(4) * t109 - qJ(5) * t115 + t22, -pkin(4) * t113 + t217 + (qJDD(4) + t110) * qJ(5) + t156, -pkin(4) * t22 + qJ(5) * t21, t147 * t64 - t207 * t98, -t144 * t167 - t147 * t52, -t144 * t83 + t231, -t144 * t63 + t206 * t96, t147 * t82 - t219, (t144 * t98 - t147 * t96) * t119, qJ(5) * t52 - t222 * t41 + t220, qJ(5) * t167 - t222 * t47 + t216, qJ(5) * t60 - t222 * t38 - t6, qJ(5) * t19 - t222 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t109, t113, t22, 0, 0, 0, 0, 0, 0, t41, t47, t38, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, t166, -t71, t159, t91, -t9, -t10, 0, 0;];
tauJ_reg  = t8;
