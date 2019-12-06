% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:08
% DurationCPUTime: 2.55s
% Computational Cost: add. (5576->332), mult. (14791->416), div. (0->0), fcn. (10620->6), ass. (0->191)
t174 = qJD(2) + qJD(3);
t173 = qJD(4) + t174;
t178 = sin(qJ(3));
t179 = sin(qJ(2));
t225 = qJD(1) * t179;
t213 = t178 * t225;
t180 = cos(qJ(3));
t181 = cos(qJ(2));
t231 = t180 * t181;
t137 = qJD(1) * t231 - t213;
t149 = t178 * t181 + t179 * t180;
t138 = qJD(1) * t149;
t177 = sin(qJ(4));
t253 = cos(qJ(4));
t98 = -t253 * t137 + t138 * t177;
t241 = t173 * t98;
t232 = t178 * t179;
t199 = t174 * t232;
t220 = qJD(1) * qJD(2);
t211 = t181 * t220;
t223 = qJD(3) * t180;
t214 = t181 * t223;
t227 = -qJD(1) * t214 - t180 * t211;
t106 = qJD(1) * t199 + t227;
t192 = t149 * qJD(3);
t193 = t149 * qJD(2);
t186 = -t193 - t192;
t185 = t186 * qJD(1);
t212 = qJD(4) * t253;
t222 = qJD(4) * t177;
t30 = t253 * t106 - t137 * t212 + t138 * t222 - t177 * t185;
t264 = -t30 + t241;
t195 = t177 * t137 + t253 * t138;
t255 = t195 ^ 2;
t96 = t98 ^ 2;
t263 = -t96 + t255;
t170 = -t181 * pkin(2) - pkin(1);
t158 = qJD(1) * t170;
t117 = -pkin(3) * t137 + t158;
t210 = pkin(4) * t98 + qJD(5);
t59 = t117 + t210;
t248 = t59 * t195;
t240 = t30 * qJ(5);
t254 = -pkin(7) - pkin(6);
t159 = t254 * t179;
t153 = qJD(1) * t159;
t244 = qJD(2) * pkin(2);
t145 = t153 + t244;
t216 = qJD(2) * t254;
t203 = qJD(1) * t216;
t146 = t179 * t203;
t147 = t181 * t203;
t160 = t254 * t181;
t155 = qJD(1) * t160;
t224 = qJD(3) * t178;
t57 = t145 * t223 + t180 * t146 + t178 * t147 + t155 * t224;
t38 = pkin(8) * t185 + t57;
t143 = t180 * t155;
t108 = t145 * t178 - t143;
t205 = -t146 * t178 + t180 * t147;
t58 = -t108 * qJD(3) + t205;
t39 = pkin(8) * t106 + t58;
t209 = -t177 * t38 + t253 * t39;
t139 = t178 * t155;
t107 = t180 * t145 + t139;
t132 = t138 * pkin(8);
t77 = t107 - t132;
t66 = pkin(3) * t174 + t77;
t250 = pkin(8) * t137;
t78 = t108 + t250;
t72 = t253 * t78;
t34 = t177 * t66 + t72;
t7 = -t34 * qJD(4) + t209;
t188 = t7 + t240;
t221 = t195 * qJD(5);
t3 = t188 - t221;
t262 = t3 - t248;
t242 = t117 * t195;
t261 = t7 - t242;
t239 = t98 * qJ(5);
t260 = t98 * t195;
t237 = t195 * t173;
t184 = t253 * t186;
t31 = -qJD(1) * t184 + t195 * qJD(4) - t177 * t106;
t259 = -t31 + t237;
t208 = -t177 * t39 - t66 * t212 + t78 * t222 - t253 * t38;
t196 = t117 * t98 + t208;
t172 = t179 * t244;
t258 = 0.2e1 * t172;
t257 = -0.2e1 * t220;
t90 = t195 * qJ(5);
t119 = t178 * t159 - t180 * t160;
t252 = pkin(2) * t179;
t251 = pkin(3) * t138;
t249 = pkin(8) * t149;
t70 = t177 * t78;
t33 = t253 * t66 - t70;
t18 = t33 - t90;
t17 = pkin(4) * t173 + t18;
t246 = t17 - t18;
t204 = pkin(3) * t212;
t245 = -t177 * pkin(3) * t31 - t98 * t204;
t41 = t253 * t77 - t70;
t114 = -t153 * t178 + t143;
t80 = t114 - t250;
t115 = t180 * t153 + t139;
t81 = -t132 + t115;
t45 = t177 * t80 + t253 * t81;
t233 = t178 * t160;
t118 = t180 * t159 + t233;
t91 = t118 - t249;
t198 = -t231 + t232;
t92 = -t198 * pkin(8) + t119;
t54 = t177 * t91 + t253 * t92;
t236 = t138 * t137;
t235 = t138 * t158;
t234 = t177 * t178;
t183 = qJD(1) ^ 2;
t230 = t181 * t183;
t182 = qJD(2) ^ 2;
t229 = t182 * t179;
t228 = t182 * t181;
t226 = t179 ^ 2 - t181 ^ 2;
t169 = pkin(2) * t180 + pkin(3);
t104 = t169 * t212 + (-t178 * t222 + (t253 * t180 - t234) * qJD(3)) * pkin(2);
t215 = t253 * t178;
t105 = -t169 * t222 + (-t178 * t212 + (-t177 * t180 - t215) * qJD(3)) * pkin(2);
t134 = pkin(2) * t215 + t177 * t169;
t219 = -t104 * t98 - t105 * t195 - t134 * t31;
t171 = pkin(2) * t225;
t218 = t179 * t230;
t154 = t179 * t216;
t156 = t181 * t216;
t217 = t180 * t154 + t178 * t156 + t159 * t223;
t40 = -t177 * t77 - t72;
t44 = -t177 * t81 + t253 * t80;
t53 = -t177 * t92 + t253 * t91;
t206 = pkin(1) * t257;
t202 = t179 * t211;
t133 = -pkin(2) * t234 + t253 * t169;
t65 = pkin(4) * t195 + t251;
t19 = t34 - t239;
t201 = -t17 * t98 + t19 * t195;
t200 = t195 * t34 - t33 * t98;
t197 = qJ(5) * t31 + t208;
t51 = -pkin(8) * t193 + (t233 - t249) * qJD(3) + t217;
t116 = -qJD(2) * t231 + t199 - t214;
t64 = -t119 * qJD(3) - t154 * t178 + t180 * t156;
t52 = pkin(8) * t116 + t64;
t10 = t177 * t52 + t91 * t212 - t92 * t222 + t253 * t51;
t194 = -t137 * t158 - t57;
t191 = t253 * t198;
t2 = -qJD(5) * t98 - t197;
t122 = t198 * pkin(3) + t170;
t190 = t59 * t98 - t2;
t113 = t253 * t149 - t177 * t198;
t11 = -t54 * qJD(4) - t177 * t51 + t253 * t52;
t187 = (-t72 + (-pkin(3) * t173 - t66) * t177) * qJD(4) + t209;
t103 = -t186 * pkin(3) + t172;
t82 = -pkin(3) * t185 + qJD(2) * t171;
t12 = t31 * pkin(4) + t82;
t168 = t253 * pkin(3) + pkin(4);
t157 = t173 * t204;
t130 = pkin(4) + t133;
t120 = t171 + t251;
t112 = t149 * t177 + t191;
t89 = t105 * t173;
t88 = t104 * t173;
t79 = -t137 ^ 2 + t138 ^ 2;
t73 = t112 * pkin(4) + t122;
t69 = t138 * t174 + t185;
t68 = -t227 + (-t137 - t213) * t174;
t63 = t160 * t224 + t217;
t60 = t171 + t65;
t47 = t113 * qJD(4) - t177 * t116 - t184;
t46 = qJD(4) * t191 + t253 * t116 + t149 * t222 - t177 * t186;
t43 = t47 * t173;
t42 = t46 * t173;
t29 = -qJ(5) * t112 + t54;
t28 = -qJ(5) * t113 + t53;
t25 = t47 * pkin(4) + t103;
t24 = -t90 + t45;
t23 = t44 + t239;
t22 = -t90 + t41;
t21 = t40 + t239;
t9 = t112 * t31 + t47 * t98;
t8 = -t113 * t30 - t195 * t46;
t5 = t46 * qJ(5) - t113 * qJD(5) + t11;
t4 = -qJ(5) * t47 - qJD(5) * t112 + t10;
t1 = t112 * t30 - t113 * t31 - t195 * t47 + t46 * t98;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t202, t226 * t257, t228, -0.2e1 * t202, -t229, 0, -pkin(6) * t228 + t179 * t206, pkin(6) * t229 + t181 * t206, 0, 0, -t106 * t149 - t116 * t138, t106 * t198 - t116 * t137 + t138 * t186 + t149 * t185, -t116 * t174, t137 * t186 - t198 * t185, t186 * t174, 0, t64 * t174 + t158 * t192 + (-t137 * t252 + t158 * t149) * qJD(2) + (t170 * t192 + (t170 * t149 + t198 * t252) * qJD(2)) * qJD(1), -t106 * t170 - t116 * t158 + t138 * t258 - t174 * t63, t118 * t106 + t107 * t116 + t108 * t186 + t119 * t185 + t63 * t137 - t64 * t138 - t58 * t149 - t57 * t198, t107 * t64 + t108 * t63 + t118 * t58 + t119 * t57 + t158 * t258, t8, t1, -t42, t9, -t43, 0, t103 * t98 + t11 * t173 + t112 * t82 + t117 * t47 + t122 * t31, -t10 * t173 + t103 * t195 + t113 * t82 - t117 * t46 - t122 * t30, -t10 * t98 - t11 * t195 + t112 * t208 - t113 * t7 + t30 * t53 - t31 * t54 + t33 * t46 - t34 * t47, t10 * t34 + t103 * t117 + t11 * t33 + t122 * t82 - t208 * t54 + t53 * t7, t8, t1, -t42, t9, -t43, 0, t112 * t12 + t173 * t5 + t25 * t98 + t31 * t73 + t47 * t59, t113 * t12 - t173 * t4 + t195 * t25 - t30 * t73 - t46 * t59, -t112 * t2 - t113 * t3 + t17 * t46 - t19 * t47 - t195 * t5 + t28 * t30 - t29 * t31 - t4 * t98, t12 * t73 + t17 * t5 + t19 * t4 + t2 * t29 + t25 * t59 + t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, t226 * t183, 0, t218, 0, 0, t183 * pkin(1) * t179, pkin(1) * t230, 0, 0, -t236, t79, t68, t236, t69, 0, t137 * t171 - t114 * t174 - t235 + (t143 + (-pkin(2) * t174 - t145) * t178) * qJD(3) + t205, t115 * t174 + (-t138 * t225 - t174 * t223) * pkin(2) + t194, (t108 + t114) * t138 + (t107 - t115) * t137 + (t106 * t180 + (t180 * t137 + t178 * t138) * qJD(3) + t178 * t185) * pkin(2), -t107 * t114 - t108 * t115 + (-t158 * t225 + t57 * t178 + t180 * t58 + (-t107 * t178 + t108 * t180) * qJD(3)) * pkin(2), t260, t263, t264, -t260, t259, 0, -t120 * t98 - t44 * t173 + t261 + t89, -t120 * t195 + t173 * t45 + t196 - t88, t133 * t30 + t195 * t44 + t45 * t98 + t200 + t219, -t117 * t120 + t133 * t7 - t134 * t208 + (t104 - t45) * t34 + (t105 - t44) * t33, t260, t263, t264, -t260, t259, 0, -t23 * t173 - t60 * t98 + t262 + t89, t173 * t24 - t195 * t60 + t190 - t88, t130 * t30 + t195 * t23 + t24 * t98 + t201 + t219, t130 * t3 + t134 * t2 - t59 * t60 + (t104 - t24) * t19 + (t105 - t23) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t79, t68, t236, t69, 0, t108 * t174 - t235 + t58, t107 * t174 + t194, 0, 0, t260, t263, t264, -t260, t259, 0, -t40 * t173 - t98 * t251 + t187 - t242, t173 * t41 - t195 * t251 - t157 + t196, t40 * t195 + t41 * t98 + (t195 * t222 + t253 * t30) * pkin(3) + t200 + t245, -t33 * t40 - t34 * t41 + (t253 * t7 - t117 * t138 - t177 * t208 + (-t177 * t33 + t253 * t34) * qJD(4)) * pkin(3), t260, t263, t264, -t260, t259, 0, -t21 * t173 - t65 * t98 + t187 - t221 + t240 - t248, t173 * t22 - t195 * t65 - t157 + t190, t168 * t30 + t22 * t98 + (pkin(3) * t222 + t21) * t195 + t201 + t245, t3 * t168 - t17 * t21 - t19 * t22 - t59 * t65 + (t177 * t2 + (-t17 * t177 + t253 * t19) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t263, t264, -t260, t259, 0, t34 * t173 + t261, t173 * t33 + t196, 0, 0, t260, t263, t264, -t260, t259, 0, t19 * t173 + (-t210 - t59) * t195 + t188, -pkin(4) * t255 + t173 * t18 + (qJD(5) + t59) * t98 + t197, pkin(4) * t30 - t246 * t98, t262 * pkin(4) + t246 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t237, -t30 - t241, -t96 - t255, t17 * t195 + t19 * t98 + t12;];
tauc_reg = t6;
