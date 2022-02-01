% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP4
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:53
% EndTime: 2022-01-23 09:32:59
% DurationCPUTime: 2.79s
% Computational Cost: add. (7259->299), mult. (17799->378), div. (0->0), fcn. (12163->8), ass. (0->183)
t234 = 2 * qJD(2);
t151 = cos(pkin(8));
t195 = t151 * qJDD(1);
t138 = -qJDD(3) + t195;
t131 = -qJDD(4) + t138;
t152 = sin(qJ(4));
t153 = sin(qJ(3));
t155 = cos(qJ(4));
t156 = cos(qJ(3));
t150 = sin(pkin(8));
t203 = qJD(1) * t150;
t109 = (-t156 * t152 - t153 * t155) * t203;
t201 = qJD(1) * t156;
t181 = t150 * t201;
t182 = t153 * t203;
t111 = -t152 * t182 + t155 * t181;
t209 = t111 * t109;
t239 = t131 - t209;
t247 = t239 * pkin(4);
t246 = t152 * t239;
t245 = t155 * t239;
t196 = t150 * qJDD(1);
t135 = t156 * t196;
t119 = -qJD(3) * t182 + t135;
t202 = qJD(1) * t151;
t139 = -qJD(3) + t202;
t146 = t150 ^ 2;
t158 = qJD(1) ^ 2;
t206 = t146 * t158;
t188 = t156 * t206;
t171 = -t151 * pkin(2) - t150 * pkin(6);
t125 = t171 * qJD(1);
t154 = sin(qJ(1));
t157 = cos(qJ(1));
t170 = g(1) * t157 + g(2) * t154;
t123 = -pkin(1) * t158 + qJDD(1) * qJ(2) - t170;
t240 = qJD(1) * t234 + t123;
t174 = -g(3) * t150 + t240 * t151;
t77 = t125 * t202 + t174;
t180 = g(1) * t154 - t157 * g(2);
t214 = qJ(2) * t158;
t163 = qJDD(2) - t180 - t214;
t168 = -pkin(1) + t171;
t99 = t168 * qJDD(1) + t163;
t89 = t156 * t99;
t43 = -pkin(3) * t138 - pkin(7) * t119 + t89 + (pkin(7) * t139 * t203 - pkin(3) * t188 - t77) * t153;
t115 = -pkin(3) * t139 - pkin(7) * t181;
t118 = (qJD(3) * t201 + qJDD(1) * t153) * t150;
t137 = t153 ^ 2 * t206;
t62 = t153 * t99 + t156 * t77;
t44 = -pkin(3) * t137 - pkin(7) * t118 + t115 * t139 + t62;
t24 = t152 * t44 - t155 * t43;
t132 = -qJD(4) + t139;
t210 = t109 * t132;
t169 = -t118 * t152 + t119 * t155;
t68 = qJD(4) * t109 + t169;
t59 = t210 + t68;
t241 = qJ(5) * t59;
t242 = 2 * qJD(5);
t162 = t111 * t242 + t24 + t241 + t247;
t161 = -t162 - t247;
t107 = t109 ^ 2;
t130 = t132 ^ 2;
t71 = -t130 - t107;
t47 = t152 * t71 - t245;
t46 = pkin(3) * t47;
t244 = t161 + t46;
t177 = -t155 * t118 - t119 * t152;
t166 = qJD(4) * t111 - t177;
t25 = t152 * t43 + t155 * t44;
t87 = -pkin(4) * t132 - qJ(5) * t111;
t14 = -t107 * pkin(4) - qJ(5) * t166 + t109 * t242 + t132 * t87 + t25;
t108 = t111 ^ 2;
t85 = -t108 - t130;
t167 = pkin(4) * t85 - t14;
t72 = t131 + t209;
t217 = t152 * t72;
t52 = t155 * t85 + t217;
t51 = pkin(3) * t52;
t243 = t51 + t167;
t124 = t139 * t181;
t92 = t124 - t118;
t238 = t150 * t92;
t211 = qJDD(1) * pkin(1);
t120 = -t163 + t211;
t147 = t151 ^ 2;
t237 = qJ(2) * t206 + t147 * t214 - t120 - t211;
t236 = t182 * (qJD(3) - t139) - t135;
t235 = pkin(4) * t166 - t107 * qJ(5) + t111 * t87 + qJDD(5);
t136 = t139 ^ 2;
t7 = t152 * t25 - t155 * t24;
t233 = pkin(3) * t7;
t48 = t155 * t71 + t246;
t30 = t153 * t48 + t156 * t47;
t232 = pkin(2) * t30;
t215 = t155 * t72;
t53 = -t152 * t85 + t215;
t32 = t153 * t53 + t156 * t52;
t231 = pkin(2) * t32;
t57 = (-qJD(4) - t132) * t111 + t177;
t35 = t152 * t57 - t155 * t59;
t36 = t152 * t59 + t155 * t57;
t16 = t153 * t36 + t156 * t35;
t230 = pkin(6) * t16;
t229 = pkin(6) * t30;
t228 = pkin(6) * t32;
t227 = pkin(7) * t35;
t226 = pkin(7) * t47;
t225 = pkin(7) * t52;
t224 = g(3) * t151;
t223 = t156 * t7;
t69 = -t107 - t108;
t222 = qJ(2) * (t151 * (-t153 * t35 + t156 * t36) + t150 * t69) - pkin(1) * t16;
t98 = t132 * t111;
t55 = t166 - t98;
t221 = qJ(2) * (t151 * (-t153 * t47 + t156 * t48) + t150 * t55) - pkin(1) * t30;
t60 = -(-qJD(4) + t132) * t109 + t169;
t220 = qJ(2) * (t151 * (-t153 * t52 + t156 * t53) + t150 * t60) - pkin(1) * t32;
t173 = t118 * pkin(3) - pkin(7) * t137 + t224;
t176 = -t115 * t156 - t125;
t49 = (t123 + (t234 - t176) * qJD(1)) * t150 + t173;
t218 = t152 * t49;
t216 = t155 * t49;
t213 = qJ(5) * t152;
t212 = qJ(5) * t155;
t208 = t132 * t152;
t207 = t132 * t155;
t129 = t153 * t188;
t116 = -t129 + t138;
t205 = t153 * t116;
t117 = -t129 - t138;
t204 = t156 * t117;
t194 = t51 - t25;
t11 = pkin(4) * t162;
t3 = t14 * t152 - t155 * t162;
t192 = pkin(3) * t3 - t11;
t189 = t156 ^ 2 * t206;
t187 = t151 * t209;
t34 = pkin(3) * t35;
t186 = -pkin(2) * t16 - t34;
t185 = -pkin(3) * t69 + pkin(7) * t36;
t184 = -pkin(3) * t55 + pkin(7) * t48;
t183 = -pkin(3) * t60 + pkin(7) * t53;
t179 = t150 * (t240 * t150 + t224) + t151 * t174;
t8 = t152 * t24 + t155 * t25;
t61 = t153 * t77 - t89;
t175 = t24 - t46;
t37 = t153 * t62 - t156 * t61;
t28 = t49 + t235;
t144 = t147 * qJDD(1);
t143 = t146 * qJDD(1);
t126 = (t146 + t147) * t158;
t122 = -t137 + t189;
t121 = -t137 - t136;
t104 = -t189 - t136;
t94 = -t108 + t130;
t93 = t107 - t130;
t91 = t124 + t118;
t90 = -t135 + (qJD(3) + t139) * t182;
t84 = t121 * t153 + t204;
t79 = t108 - t107;
t78 = t104 * t156 + t205;
t76 = t224 + (t123 + (t234 + t125) * qJD(1)) * t150;
t63 = -t153 * t91 + t156 * t90;
t58 = t68 - t210;
t56 = t166 + t98;
t54 = pkin(4) * t59;
t39 = -pkin(4) * t60 + qJ(5) * t72;
t38 = t151 * t131 + t150 * (t156 * (-t109 * t155 - t111 * t152) - t153 * (-t109 * t152 + t111 * t155)) * t132;
t27 = t150 * (t156 * (t111 * t208 + t155 * t68) - t153 * (-t111 * t207 + t152 * t68)) + t187;
t26 = t150 * (t156 * (t109 * t207 + t152 * t166) - t153 * (t109 * t208 - t155 * t166)) - t187;
t22 = -qJ(5) * t85 + t28;
t21 = t150 * (t156 * (t155 * t93 + t217) - t153 * (t152 * t93 - t215)) + t151 * t56;
t20 = t150 * (t156 * (-t152 * t94 - t245) - t153 * (t155 * t94 - t246)) - t151 * t59;
t17 = -pkin(4) * t55 + qJ(5) * t71 - t173 + (t176 * qJD(1) - t240) * t150 - t235;
t13 = t150 * (t156 * (-t152 * t58 - t155 * t55) - t153 * (-t152 * t55 + t155 * t58)) - t151 * t79;
t9 = t162 + t241;
t6 = -pkin(4) * t69 + qJ(5) * t57 + t14;
t5 = -pkin(4) * t28 + qJ(5) * t14;
t4 = t14 * t155 + t152 * t162;
t2 = t153 * t8 + t223;
t1 = t153 * t4 + t156 * t3;
t10 = [0, 0, 0, 0, 0, qJDD(1), t180, t170, 0, 0, t143, 0.2e1 * t150 * t195, 0, t144, 0, 0, -t237 * t151, t237 * t150, pkin(1) * t126 + qJ(2) * (t144 + t143) + t179, pkin(1) * t120 + qJ(2) * t179, (t150 * (t139 * t182 + t119) - t151 * t153 * t206) * t156, t150 * (t153 * t236 + t156 * t92) - t151 * t122, t150 * (t204 - t153 * (t136 - t189)) + t151 * t90, (t151 * t188 - t238) * t153, t150 * (t156 * (t137 - t136) + t205) + t151 * t91, t151 * t138, t150 * (-pkin(6) * t84 + t153 * t76) + t151 * (-pkin(2) * t84 + t61) - pkin(1) * t84 + qJ(2) * (t151 * (-t117 * t153 + t121 * t156) - t238), t150 * (-pkin(6) * t78 + t156 * t76) + t151 * (-pkin(2) * t78 + t62) - pkin(1) * t78 + qJ(2) * (t151 * (-t104 * t153 + t116 * t156) - t236 * t150), -t150 * t37 + qJ(2) * (t151 * (-t153 * t90 - t156 * t91) - t150 * (t137 + t189)) + t168 * t63, qJ(2) * (t151 * (t153 * t61 + t156 * t62) + t150 * t76) + t168 * t37, t27, t13, t20, t26, t21, t38, t150 * (t156 * (t218 - t226) - t153 * (t184 - t216) - t229) + t151 * (t175 - t232) + t221, t150 * (t156 * (t216 - t225) - t153 * (t183 + t218) - t228) + t151 * (-t194 - t231) + t220, t150 * (t156 * (-t7 - t227) - t153 * (t185 + t8) - t230) + t151 * t186 + t222, t150 * (-pkin(7) * t223 - t153 * (-pkin(3) * t49 + pkin(7) * t8) - pkin(6) * t2) + t151 * (-pkin(2) * t2 - t233) - pkin(1) * t2 + qJ(2) * (t151 * (-t153 * t7 + t156 * t8) + t150 * t49), t27, t13, t20, t26, t21, t38, t150 * (t156 * (-t152 * t17 + t212 * t239 - t226) - t153 * (t155 * t17 + t213 * t239 + t184) - t229) + t151 * (-t232 - t244) + t221, t150 * (t156 * (-t152 * t39 + t155 * t22 - t225) - t153 * (t152 * t22 + t155 * t39 + t183) - t228) + t151 * (-t231 - t243) + t220, t150 * (t156 * (-t152 * t6 + t155 * t9 - t227) - t153 * (t152 * t9 + t155 * t6 + t185) - t230) + t151 * (t186 + t54) + t222, t150 * (t156 * (-pkin(7) * t3 - t152 * t5 + t162 * t212) - t153 * (-pkin(3) * t28 + pkin(7) * t4 + t155 * t5 + t162 * t213) - pkin(6) * t1) + t151 * (-pkin(2) * t1 - t192) - pkin(1) * t1 + qJ(2) * (t151 * (-t153 * t3 + t156 * t4) + t150 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, t196, -t126, -t120, 0, 0, 0, 0, 0, 0, t84, t78, t63, t37, 0, 0, 0, 0, 0, 0, t30, t32, t16, t2, 0, 0, 0, 0, 0, 0, t30, t32, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t122, -t90, -t129, -t91, -t138, -t61, -t62, 0, 0, -t209, t79, t59, t209, -t56, -t131, -t175, t194, t34, t233, -t209, t79, t59, t209, -t56, -t131, t244, t243, -t54 + t34, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t79, t59, t209, -t56, -t131, -t24, -t25, 0, 0, -t209, t79, t59, t209, -t56, -t131, t161, t167, -t54, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t60, t69, t28;];
tauJ_reg = t10;
