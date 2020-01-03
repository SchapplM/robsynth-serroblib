% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:16
% EndTime: 2019-12-31 18:49:22
% DurationCPUTime: 2.64s
% Computational Cost: add. (4345->336), mult. (10702->399), div. (0->0), fcn. (8230->12), ass. (0->177)
t146 = sin(qJ(4));
t149 = cos(qJ(3));
t143 = sin(pkin(8));
t144 = cos(pkin(8));
t147 = sin(qJ(3));
t207 = qJD(1) * qJD(3);
t196 = t149 * t207;
t197 = t147 * t207;
t206 = t143 * qJDD(1);
t201 = -t143 * t196 - t144 * t197 - t147 * t206;
t205 = t144 * qJDD(1);
t169 = t149 * t205 + t201;
t240 = cos(qJ(4));
t216 = t149 * t144;
t217 = t147 * t143;
t179 = -t216 + t217;
t89 = t179 * qJD(1);
t97 = t149 * t143 + t147 * t144;
t90 = t97 * qJD(1);
t176 = t146 * t89 - t240 * t90;
t200 = t144 * t196 + t147 * t205 + t149 * t206;
t64 = t143 * t197 - t200;
t157 = qJD(4) * t176 + t146 * t64 + t240 * t169;
t142 = qJD(3) + qJD(4);
t254 = t142 * t176;
t261 = t157 - t254;
t198 = qJD(4) * t240;
t210 = qJD(4) * t146;
t175 = t146 * t169 - t89 * t198 - t90 * t210 - t240 * t64;
t58 = -t146 * t90 - t240 * t89;
t227 = t58 * t142;
t260 = t175 - t227;
t237 = t58 ^ 2;
t259 = t176 ^ 2;
t184 = t259 - t237;
t127 = t144 * pkin(2) + pkin(1);
t102 = -t127 * qJD(1) + qJD(2);
t69 = t89 * pkin(3) + t102;
t27 = -pkin(4) * t58 + qJ(5) * t176 + t69;
t258 = t27 * t58;
t257 = t69 * t58;
t256 = t58 * t176;
t232 = pkin(6) + qJ(2);
t110 = t232 * t143;
t111 = t232 * t144;
t71 = -t147 * t110 + t149 * t111;
t255 = t179 * pkin(7) - t71;
t223 = qJDD(1) * pkin(1);
t148 = sin(qJ(1));
t150 = cos(qJ(1));
t250 = g(1) * t148 - g(2) * t150;
t178 = -qJDD(2) + t223 + t250;
t141 = pkin(8) + qJ(3);
t135 = qJ(4) + t141;
t125 = cos(t135);
t124 = sin(t135);
t221 = t124 * t150;
t222 = t124 * t148;
t208 = qJD(1) * qJD(2);
t245 = t232 * qJDD(1) + t208;
t77 = t245 * t143;
t78 = t245 * t144;
t193 = -t147 * t78 - t149 * t77;
t98 = qJD(1) * t110;
t99 = qJD(1) * t111;
t68 = -t147 * t98 + t149 * t99;
t36 = -t68 * qJD(3) + t193;
t18 = qJDD(3) * pkin(3) + t64 * pkin(7) + t36;
t211 = qJD(3) * t149;
t204 = -t147 * t77 + t149 * t78 - t98 * t211;
t229 = t147 * t99;
t35 = -qJD(3) * t229 + t204;
t21 = pkin(7) * t169 + t35;
t67 = -t149 * t98 - t229;
t47 = -t90 * pkin(7) + t67;
t46 = qJD(3) * pkin(3) + t47;
t48 = -t89 * pkin(7) + t68;
t4 = -t146 * t21 + t240 * t18 - t48 * t198 - t46 * t210;
t167 = g(1) * t221 + g(2) * t222 - g(3) * t125 + t4;
t137 = qJDD(3) + qJDD(4);
t132 = t137 * pkin(4);
t247 = t132 - qJDD(5);
t155 = -t176 * t27 - t167 - t247;
t253 = t69 * t176 + t167;
t34 = -pkin(4) * t176 - t58 * qJ(5);
t252 = t125 * pkin(4) + t124 * qJ(5);
t129 = t137 * qJ(5);
t130 = t142 * qJD(5);
t251 = t129 + t130;
t172 = t97 * qJD(2);
t249 = t255 * qJD(3) - t172;
t248 = qJ(2) * qJDD(1);
t241 = t97 * pkin(7);
t218 = t147 * t111;
t70 = -t149 * t110 - t218;
t53 = t70 - t241;
t33 = t146 * t53 - t240 * t255;
t177 = t146 * t255 + t240 * t53;
t209 = t143 * qJD(2);
t224 = qJD(2) * t216 - t110 * t211;
t40 = -t147 * t209 + (-t218 - t241) * qJD(3) + t224;
t8 = t177 * qJD(4) + t249 * t146 + t240 * t40;
t246 = t124 * t250 + t33 * t137 + t8 * t142;
t243 = t90 ^ 2;
t242 = qJD(3) ^ 2;
t233 = t90 * t89;
t230 = t146 * t48;
t202 = t240 * t48;
t26 = t146 * t46 + t202;
t228 = t26 * t142;
t29 = t240 * t47 - t230;
t225 = pkin(3) * t198 + qJD(5) - t29;
t220 = t125 * t148;
t219 = t125 * t150;
t215 = t90 * qJD(3);
t25 = t240 * t46 - t230;
t214 = qJD(5) - t25;
t139 = t143 ^ 2;
t140 = t144 ^ 2;
t212 = t139 + t140;
t199 = qJD(1) * t217;
t3 = t146 * t18 + t46 * t198 + t240 * t21 - t48 * t210;
t192 = t212 * qJD(1) ^ 2;
t191 = 0.2e1 * t212;
t28 = t146 * t47 + t202;
t190 = pkin(3) * t210 - t28;
t133 = sin(t141);
t189 = -pkin(3) * t133 - pkin(4) * t124;
t188 = g(1) * t150 + g(2) * t148;
t170 = t179 * qJD(3);
t171 = t97 * qJD(3);
t66 = -t146 * t179 + t240 * t97;
t38 = t66 * qJD(4) - t146 * t170 + t240 * t171;
t166 = t240 * t179;
t65 = t146 * t97 + t166;
t186 = -t157 * t65 - t38 * t58;
t185 = -t259 - t237;
t181 = t137 * t65 + t142 * t38;
t173 = -g(1) * t219 - g(2) * t220 - g(3) * t124 + t3;
t101 = -t127 * qJDD(1) + qJDD(2);
t168 = pkin(3) * t171;
t165 = t178 + t223;
t164 = t175 + t227;
t163 = t25 * t142 - t173;
t134 = cos(t141);
t162 = -g(3) * t134 + t133 * t188;
t9 = t33 * qJD(4) + t146 * t40 - t240 * t249;
t161 = g(1) * t220 - g(2) * t219 + t137 * t177 - t9 * t142;
t37 = t142 * t166 + t146 * t171 + t97 * t210;
t160 = t157 * t66 - t175 * t65 + t176 * t38 - t37 * t58;
t75 = pkin(3) * t179 - t127;
t158 = t157 * t33 - t175 * t177 - t176 * t9 + t58 * t8 - t188;
t156 = t191 * t208 - t188;
t153 = -t157 - t254;
t49 = qJDD(2) - t201 * pkin(3) + (-pkin(1) + (-t149 * pkin(3) - pkin(2)) * t144) * qJDD(1);
t5 = -pkin(4) * t157 - qJ(5) * t175 + qJD(5) * t176 + t49;
t138 = -pkin(7) - t232;
t128 = -t240 * pkin(3) - pkin(4);
t126 = t146 * pkin(3) + qJ(5);
t123 = pkin(3) * t134;
t104 = qJ(5) * t219;
t103 = qJ(5) * t220;
t100 = t123 + t127;
t92 = t150 * t100;
t86 = t89 ^ 2;
t51 = -qJD(3) * t71 - t172;
t50 = (-qJD(3) * t111 - t209) * t147 + t224;
t31 = t65 * pkin(4) - t66 * qJ(5) + t75;
t30 = t90 * pkin(3) + t34;
t24 = t142 * qJ(5) + t26;
t20 = -t142 * pkin(4) + t214;
t19 = t66 * t137 - t37 * t142;
t11 = t38 * pkin(4) + t37 * qJ(5) - t66 * qJD(5) + t168;
t6 = t175 * t66 + t176 * t37;
t2 = -t4 - t247;
t1 = t3 + t251;
t7 = [0, 0, 0, 0, 0, qJDD(1), t250, t188, 0, 0, t139 * qJDD(1), 0.2e1 * t143 * t205, 0, t140 * qJDD(1), 0, 0, t165 * t144, -t165 * t143, t191 * t248 + t156, pkin(1) * t178 + (t212 * t248 + t156) * qJ(2), -t170 * t90 - t64 * t97, t97 * t169 + t64 * t179 + (t179 * t89 - t90 * t97) * qJD(3), t97 * qJDD(3) - t242 * t179, -t169 * t179 + t171 * t89, -qJDD(3) * t179 - t242 * t97, 0, t127 * t169 + t101 * t179 + t70 * qJDD(3) + t250 * t134 + (t102 * t97 + t51) * qJD(3), -t71 * qJDD(3) + t101 * t97 + t127 * t64 - t250 * t133 + (-t102 * t179 - t50) * qJD(3), -t50 * t89 + t71 * t169 - t35 * t179 - t51 * t90 + t70 * t64 - t36 * t97 + (t179 * t67 - t68 * t97) * qJD(3) - t188, t35 * t71 + t68 * t50 + t36 * t70 + t67 * t51 - t101 * t127 - g(1) * (-t148 * t127 + t150 * t232) - g(2) * (t150 * t127 + t148 * t232), t6, t160, t19, t186, -t181, 0, -t157 * t75 - t168 * t58 + t69 * t38 + t49 * t65 + t161, -t168 * t176 + t175 * t75 - t69 * t37 + t49 * t66 - t246, t25 * t37 - t26 * t38 - t3 * t65 - t4 * t66 + t158, t3 * t33 + t26 * t8 + t4 * t177 - t25 * t9 + t49 * t75 - g(1) * (-t148 * t100 - t150 * t138) - g(2) * (-t148 * t138 + t92) + t69 * t168, t6, t19, -t160, 0, t181, t186, -t11 * t58 - t157 * t31 + t27 * t38 + t5 * t65 + t161, -t1 * t65 + t2 * t66 - t20 * t37 - t24 * t38 + t158, t11 * t176 - t175 * t31 + t27 * t37 - t5 * t66 + t246, -g(2) * t92 + t1 * t33 + t27 * t11 - t2 * t177 + t20 * t9 + t24 * t8 + t5 * t31 + (g(1) * t138 - g(2) * t252) * t150 + (-g(1) * (-t100 - t252) + g(2) * t138) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t206, -t192, -qJ(2) * t192 - t178, 0, 0, 0, 0, 0, 0, -t169 + t215, (-t89 - t199) * qJD(3) + t200, -t86 - t243, t67 * t90 + t68 * t89 + t101 - t250, 0, 0, 0, 0, 0, 0, t153, t164, t185, -t176 * t25 - t26 * t58 - t250 + t49, 0, 0, 0, 0, 0, 0, t153, t185, -t164, t176 * t20 - t24 * t58 - t250 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, -t86 + t243, (t89 - t199) * qJD(3) + t200, -t233, t169 + t215, qJDD(3), -t102 * t90 + t162 + t193, g(3) * t133 + t102 * t89 + t188 * t134 + (t67 + t229) * qJD(3) - t204, 0, 0, t256, t184, t260, -t256, t261, t137, t28 * t142 + (t240 * t137 - t142 * t210 + t58 * t90) * pkin(3) + t253, t29 * t142 - t257 + (-t137 * t146 - t142 * t198 + t176 * t90) * pkin(3) - t173, t25 * t58 - t26 * t176 + t28 * t176 - t29 * t58 + (-t240 * t175 + t146 * t157 + (-t146 * t176 + t240 * t58) * qJD(4)) * pkin(3), t25 * t28 - t26 * t29 + (t240 * t4 + t146 * t3 - t69 * t90 + (-t146 * t25 + t240 * t26) * qJD(4) + t162) * pkin(3), t256, t260, -t184, t137, -t261, -t256, -t128 * t137 - t142 * t190 + t30 * t58 - t155, t126 * t157 + t128 * t175 + (-t20 + t225) * t58 + (-t190 - t24) * t176, t126 * t137 + t225 * t142 - t176 * t30 + t173 + t251 + t258, t1 * t126 + t2 * t128 - t27 * t30 - g(1) * (t150 * t189 + t104) - g(2) * (t148 * t189 + t103) - g(3) * (t123 + t252) + t225 * t24 + t190 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t184, t260, -t256, t261, t137, t228 + t253, t163 - t257, 0, 0, t256, t260, -t184, t137, -t261, -t256, t34 * t58 + t132 - t155 + t228, -pkin(4) * t175 + t157 * qJ(5) - (t24 - t26) * t176 - (t20 - t214) * t58, -t176 * t34 + 0.2e1 * t129 + 0.2e1 * t130 - t163 + t258, t1 * qJ(5) - t2 * pkin(4) - t27 * t34 - t20 * t26 - g(1) * (-pkin(4) * t221 + t104) - g(2) * (-pkin(4) * t222 + t103) - g(3) * t252 + t214 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 + t256, t260, -t142 ^ 2 - t259, -t24 * t142 + t155;];
tau_reg = t7;
