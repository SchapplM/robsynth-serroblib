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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:06:54
% EndTime: 2019-12-05 18:07:05
% DurationCPUTime: 2.85s
% Computational Cost: add. (7259->299), mult. (17799->378), div. (0->0), fcn. (12163->8), ass. (0->183)
t233 = 2 * qJD(2);
t150 = cos(pkin(8));
t194 = t150 * qJDD(1);
t138 = -qJDD(3) + t194;
t131 = -qJDD(4) + t138;
t151 = sin(qJ(4));
t152 = sin(qJ(3));
t154 = cos(qJ(4));
t155 = cos(qJ(3));
t149 = sin(pkin(8));
t202 = qJD(1) * t149;
t109 = (-t151 * t155 - t152 * t154) * t202;
t200 = qJD(1) * t155;
t180 = t149 * t200;
t181 = t152 * t202;
t111 = -t151 * t181 + t154 * t180;
t208 = t111 * t109;
t238 = t131 - t208;
t246 = t238 * pkin(4);
t245 = t151 * t238;
t244 = t154 * t238;
t195 = t149 * qJDD(1);
t135 = t155 * t195;
t119 = -qJD(3) * t181 + t135;
t201 = qJD(1) * t150;
t139 = -qJD(3) + t201;
t145 = t149 ^ 2;
t157 = qJD(1) ^ 2;
t205 = t145 * t157;
t187 = t155 * t205;
t171 = -pkin(2) * t150 - pkin(6) * t149;
t125 = t171 * qJD(1);
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t169 = -g(2) * t153 + g(3) * t156;
t123 = -pkin(1) * t157 + qJDD(1) * qJ(2) - t169;
t239 = qJD(1) * t233 + t123;
t174 = -g(1) * t149 + t150 * t239;
t77 = t125 * t201 + t174;
t170 = g(2) * t156 + g(3) * t153;
t213 = qJ(2) * t157;
t161 = qJDD(2) - t170 - t213;
t167 = -pkin(1) + t171;
t99 = qJDD(1) * t167 + t161;
t89 = t155 * t99;
t43 = -pkin(3) * t138 - pkin(7) * t119 + t89 + (pkin(7) * t139 * t202 - pkin(3) * t187 - t77) * t152;
t115 = -pkin(3) * t139 - pkin(7) * t180;
t118 = (qJD(3) * t200 + qJDD(1) * t152) * t149;
t137 = t152 ^ 2 * t205;
t62 = t152 * t99 + t155 * t77;
t44 = -pkin(3) * t137 - pkin(7) * t118 + t115 * t139 + t62;
t24 = t151 * t44 - t154 * t43;
t132 = -qJD(4) + t139;
t209 = t109 * t132;
t168 = -t118 * t151 + t119 * t154;
t68 = qJD(4) * t109 + t168;
t59 = t209 + t68;
t240 = qJ(5) * t59;
t241 = 2 * qJD(5);
t162 = t111 * t241 + t24 + t240 + t246;
t160 = -t162 - t246;
t107 = t109 ^ 2;
t130 = t132 ^ 2;
t71 = -t130 - t107;
t47 = t151 * t71 - t244;
t46 = pkin(3) * t47;
t243 = t160 + t46;
t177 = -t118 * t154 - t119 * t151;
t165 = qJD(4) * t111 - t177;
t25 = t151 * t43 + t154 * t44;
t87 = -pkin(4) * t132 - qJ(5) * t111;
t14 = -pkin(4) * t107 - qJ(5) * t165 + t109 * t241 + t132 * t87 + t25;
t108 = t111 ^ 2;
t85 = -t108 - t130;
t166 = pkin(4) * t85 - t14;
t72 = t131 + t208;
t216 = t151 * t72;
t52 = t154 * t85 + t216;
t51 = pkin(3) * t52;
t242 = t51 + t166;
t124 = t139 * t180;
t92 = t124 - t118;
t237 = t149 * t92;
t210 = qJDD(1) * pkin(1);
t120 = -t161 + t210;
t146 = t150 ^ 2;
t236 = qJ(2) * t205 + t146 * t213 - t120 - t210;
t235 = (qJD(3) - t139) * t181 - t135;
t234 = pkin(4) * t165 - qJ(5) * t107 + t111 * t87 + qJDD(5);
t136 = t139 ^ 2;
t7 = t151 * t25 - t154 * t24;
t232 = pkin(3) * t7;
t48 = t154 * t71 + t245;
t30 = t152 * t48 + t155 * t47;
t231 = pkin(2) * t30;
t214 = t154 * t72;
t53 = -t151 * t85 + t214;
t32 = t152 * t53 + t155 * t52;
t230 = pkin(2) * t32;
t57 = (-qJD(4) - t132) * t111 + t177;
t35 = t151 * t57 - t154 * t59;
t36 = t151 * t59 + t154 * t57;
t16 = t152 * t36 + t155 * t35;
t229 = pkin(6) * t16;
t228 = pkin(6) * t30;
t227 = pkin(6) * t32;
t226 = pkin(7) * t35;
t225 = pkin(7) * t47;
t224 = pkin(7) * t52;
t223 = g(1) * t150;
t222 = t155 * t7;
t69 = -t107 - t108;
t221 = qJ(2) * (t150 * (-t152 * t35 + t155 * t36) + t149 * t69) - pkin(1) * t16;
t98 = t132 * t111;
t55 = t165 - t98;
t220 = qJ(2) * (t150 * (-t152 * t47 + t155 * t48) + t149 * t55) - pkin(1) * t30;
t60 = -(-qJD(4) + t132) * t109 + t168;
t219 = qJ(2) * (t150 * (-t152 * t52 + t155 * t53) + t149 * t60) - pkin(1) * t32;
t173 = pkin(3) * t118 - pkin(7) * t137 + t223;
t176 = -t115 * t155 - t125;
t49 = (t123 + (t233 - t176) * qJD(1)) * t149 + t173;
t217 = t151 * t49;
t215 = t154 * t49;
t212 = qJ(5) * t151;
t211 = qJ(5) * t154;
t207 = t132 * t151;
t206 = t132 * t154;
t129 = t152 * t187;
t116 = -t129 + t138;
t204 = t152 * t116;
t117 = -t129 - t138;
t203 = t155 * t117;
t193 = t51 - t25;
t11 = pkin(4) * t162;
t3 = t14 * t151 - t154 * t162;
t191 = pkin(3) * t3 - t11;
t188 = t155 ^ 2 * t205;
t186 = t150 * t208;
t34 = pkin(3) * t35;
t185 = -pkin(2) * t16 - t34;
t184 = -pkin(3) * t69 + pkin(7) * t36;
t183 = -pkin(3) * t55 + pkin(7) * t48;
t182 = -pkin(3) * t60 + pkin(7) * t53;
t179 = t149 * (t149 * t239 + t223) + t150 * t174;
t8 = t151 * t24 + t154 * t25;
t61 = t152 * t77 - t89;
t175 = t24 - t46;
t37 = t152 * t62 - t155 * t61;
t28 = t49 + t234;
t144 = t146 * qJDD(1);
t143 = t145 * qJDD(1);
t126 = (t145 + t146) * t157;
t122 = -t137 + t188;
t121 = -t137 - t136;
t104 = -t188 - t136;
t94 = -t108 + t130;
t93 = t107 - t130;
t91 = t124 + t118;
t90 = -t135 + (qJD(3) + t139) * t181;
t84 = t121 * t152 + t203;
t79 = t108 - t107;
t78 = t104 * t155 + t204;
t76 = t223 + (t123 + (t233 + t125) * qJD(1)) * t149;
t63 = -t152 * t91 + t155 * t90;
t58 = t68 - t209;
t56 = t165 + t98;
t54 = pkin(4) * t59;
t39 = -pkin(4) * t60 + qJ(5) * t72;
t38 = t150 * t131 + t149 * (t155 * (-t109 * t154 - t111 * t151) - t152 * (-t109 * t151 + t111 * t154)) * t132;
t27 = t149 * (t155 * (t111 * t207 + t154 * t68) - t152 * (-t111 * t206 + t151 * t68)) + t186;
t26 = t149 * (t155 * (t109 * t206 + t151 * t165) - t152 * (t109 * t207 - t154 * t165)) - t186;
t22 = -qJ(5) * t85 + t28;
t21 = t149 * (t155 * (t154 * t93 + t216) - t152 * (t151 * t93 - t214)) + t150 * t56;
t20 = t149 * (t155 * (-t151 * t94 - t244) - t152 * (t154 * t94 - t245)) - t150 * t59;
t17 = -pkin(4) * t55 + qJ(5) * t71 - t173 + (qJD(1) * t176 - t239) * t149 - t234;
t13 = t149 * (t155 * (-t151 * t58 - t154 * t55) - t152 * (-t151 * t55 + t154 * t58)) - t150 * t79;
t9 = t162 + t240;
t6 = -pkin(4) * t69 + qJ(5) * t57 + t14;
t5 = -pkin(4) * t28 + qJ(5) * t14;
t4 = t14 * t154 + t151 * t162;
t2 = t152 * t8 + t222;
t1 = t152 * t4 + t155 * t3;
t10 = [0, 0, 0, 0, 0, qJDD(1), t170, t169, 0, 0, t143, 0.2e1 * t149 * t194, 0, t144, 0, 0, -t236 * t150, t236 * t149, pkin(1) * t126 + qJ(2) * (t144 + t143) + t179, pkin(1) * t120 + qJ(2) * t179, (t149 * (t139 * t181 + t119) - t150 * t152 * t205) * t155, t149 * (t152 * t235 + t155 * t92) - t150 * t122, t149 * (t203 - t152 * (t136 - t188)) + t150 * t90, (t150 * t187 - t237) * t152, t149 * (t155 * (t137 - t136) + t204) + t150 * t91, t150 * t138, t149 * (-pkin(6) * t84 + t152 * t76) + t150 * (-pkin(2) * t84 + t61) - pkin(1) * t84 + qJ(2) * (t150 * (-t117 * t152 + t121 * t155) - t237), t149 * (-pkin(6) * t78 + t155 * t76) + t150 * (-pkin(2) * t78 + t62) - pkin(1) * t78 + qJ(2) * (t150 * (-t104 * t152 + t116 * t155) - t235 * t149), -t149 * t37 + qJ(2) * (t150 * (-t152 * t90 - t155 * t91) - t149 * (t137 + t188)) + t167 * t63, qJ(2) * (t150 * (t152 * t61 + t155 * t62) + t149 * t76) + t167 * t37, t27, t13, t20, t26, t21, t38, t149 * (t155 * (t217 - t225) - t152 * (t183 - t215) - t228) + t150 * (t175 - t231) + t220, t149 * (t155 * (t215 - t224) - t152 * (t182 + t217) - t227) + t150 * (-t193 - t230) + t219, t149 * (t155 * (-t7 - t226) - t152 * (t184 + t8) - t229) + t150 * t185 + t221, t149 * (-pkin(7) * t222 - t152 * (-pkin(3) * t49 + pkin(7) * t8) - pkin(6) * t2) + t150 * (-pkin(2) * t2 - t232) - pkin(1) * t2 + qJ(2) * (t150 * (-t152 * t7 + t155 * t8) + t149 * t49), t27, t13, t20, t26, t21, t38, t149 * (t155 * (-t151 * t17 + t211 * t238 - t225) - t152 * (t154 * t17 + t212 * t238 + t183) - t228) + t150 * (-t231 - t243) + t220, t149 * (t155 * (-t151 * t39 + t154 * t22 - t224) - t152 * (t151 * t22 + t154 * t39 + t182) - t227) + t150 * (-t230 - t242) + t219, t149 * (t155 * (-t151 * t6 + t154 * t9 - t226) - t152 * (t151 * t9 + t154 * t6 + t184) - t229) + t150 * (t185 + t54) + t221, t149 * (t155 * (-pkin(7) * t3 - t151 * t5 + t162 * t211) - t152 * (-pkin(3) * t28 + pkin(7) * t4 + t154 * t5 + t162 * t212) - pkin(6) * t1) + t150 * (-pkin(2) * t1 - t191) - pkin(1) * t1 + qJ(2) * (t150 * (-t152 * t3 + t155 * t4) + t149 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, t195, -t126, -t120, 0, 0, 0, 0, 0, 0, t84, t78, t63, t37, 0, 0, 0, 0, 0, 0, t30, t32, t16, t2, 0, 0, 0, 0, 0, 0, t30, t32, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t122, -t90, -t129, -t91, -t138, -t61, -t62, 0, 0, -t208, t79, t59, t208, -t56, -t131, -t175, t193, t34, t232, -t208, t79, t59, t208, -t56, -t131, t243, t242, -t54 + t34, t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t79, t59, t208, -t56, -t131, -t24, -t25, 0, 0, -t208, t79, t59, t208, -t56, -t131, t160, t166, -t54, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t60, t69, t28;];
tauJ_reg = t10;
