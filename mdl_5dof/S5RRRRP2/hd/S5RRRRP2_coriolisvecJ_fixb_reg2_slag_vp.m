% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP2
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:13
% EndTime: 2019-12-05 18:48:20
% DurationCPUTime: 1.79s
% Computational Cost: add. (3550->267), mult. (6151->330), div. (0->0), fcn. (3931->6), ass. (0->176)
t152 = sin(qJ(4));
t153 = sin(qJ(3));
t155 = cos(qJ(3));
t239 = cos(qJ(4));
t113 = t152 * t155 + t239 * t153;
t149 = qJD(1) + qJD(2);
t94 = t113 * t149;
t219 = t94 * qJ(5);
t154 = sin(qJ(2));
t224 = pkin(1) * qJD(1);
t197 = t154 * t224;
t122 = pkin(7) * t149 + t197;
t181 = pkin(8) * t149 + t122;
t85 = t181 * t155;
t74 = t152 * t85;
t84 = t181 * t153;
t77 = qJD(3) * pkin(3) - t84;
t37 = t239 * t77 - t74;
t23 = t37 - t219;
t150 = t153 ^ 2;
t151 = t155 ^ 2;
t206 = t150 + t151;
t183 = qJD(4) * t239;
t189 = t239 * t155;
t242 = -qJD(3) * t189 - t155 * t183;
t240 = -pkin(8) - pkin(7);
t190 = qJD(3) * t240;
t117 = t153 * t190;
t118 = t155 * t190;
t131 = t240 * t153;
t146 = t155 * pkin(8);
t132 = pkin(7) * t155 + t146;
t214 = t152 * t153;
t165 = t189 - t214;
t156 = cos(qJ(2));
t196 = t156 * t224;
t202 = qJD(4) * t152;
t229 = t239 * t117 + t152 * t118 + t131 * t183 - t132 * t202 - t165 * t196;
t79 = t152 * t131 + t239 * t132;
t228 = -t79 * qJD(4) + t113 * t196 - t152 * t117 + t239 * t118;
t148 = qJD(3) + qJD(4);
t241 = t94 ^ 2;
t238 = pkin(1) * t156;
t192 = t149 * t214;
t92 = -t149 * t189 + t192;
t237 = t94 * t92;
t142 = -pkin(3) * t155 - pkin(2);
t96 = t142 * t149 - t196;
t236 = t96 * t94;
t139 = pkin(1) * t154 + pkin(7);
t235 = -pkin(8) - t139;
t65 = t148 * t113;
t217 = -t65 * qJ(5) + qJD(5) * t165;
t234 = t217 + t229;
t170 = t148 * t214;
t64 = t170 + t242;
t168 = t64 * qJ(5) - t113 * qJD(5);
t233 = t168 + t228;
t22 = pkin(4) * t148 + t23;
t232 = t22 - t23;
t223 = pkin(1) * qJD(2);
t193 = qJD(1) * t223;
t137 = t154 * t193;
t205 = qJD(3) * t153;
t188 = t149 * t205;
t101 = pkin(3) * t188 + t137;
t47 = t65 * t149;
t33 = pkin(4) * t47 + t101;
t182 = pkin(4) * t92 + qJD(5);
t49 = t182 + t96;
t231 = -t165 * t33 + t49 * t65;
t230 = t33 * t113 - t49 * t64;
t227 = -t101 * t165 + t96 * t65;
t176 = pkin(3) * t183;
t226 = -t152 * pkin(3) * t47 - t92 * t176;
t225 = t101 * t113 - t96 * t64;
t40 = -t239 * t84 - t74;
t222 = t148 * t92;
t209 = t242 * t149;
t46 = t149 * t170 + t209;
t221 = t46 * qJ(5);
t220 = t92 * qJ(5);
t109 = t235 * t153;
t110 = t139 * t155 + t146;
t59 = t152 * t109 + t239 * t110;
t216 = qJ(5) * t113;
t215 = t149 * t153;
t213 = t154 * t155;
t157 = qJD(3) ^ 2;
t212 = t157 * t153;
t145 = t157 * t155;
t211 = qJD(5) + t49;
t123 = -pkin(2) * t149 - t196;
t204 = qJD(3) * t155;
t210 = t123 * t204 + t153 * t137;
t175 = t156 * t193;
t208 = t206 * t175;
t207 = t150 - t151;
t203 = qJD(3) * t156;
t201 = -qJD(1) - t149;
t200 = -qJD(2) + t149;
t199 = pkin(3) * t215;
t198 = t239 * pkin(3);
t195 = t156 * t223;
t194 = pkin(3) * t202;
t144 = t154 * t223;
t143 = pkin(3) * t205;
t76 = t239 * t85;
t147 = t149 ^ 2;
t191 = t153 * t147 * t155;
t187 = t149 * t204;
t186 = t153 * t203;
t38 = t152 * t77 + t76;
t24 = t38 - t220;
t169 = qJD(3) * t181;
t56 = -t153 * t169 + t155 * t175;
t57 = -t153 * t175 - t155 * t169;
t178 = -t152 * t57 - t77 * t183 + t85 * t202 - t239 * t56;
t166 = qJ(5) * t47 + t178;
t4 = -qJD(5) * t92 - t166;
t180 = -t152 * t56 + t239 * t57;
t11 = -qJD(4) * t38 + t180;
t159 = t11 + t221;
t5 = -t94 * qJD(5) + t159;
t185 = -t5 * t113 + t165 * t4 + t22 * t64 - t24 * t65;
t184 = -t11 * t113 - t165 * t178 + t37 * t64 - t38 * t65;
t55 = pkin(4) * t65 + t143;
t39 = t152 * t84 - t76;
t179 = qJD(3) * t235;
t58 = t239 * t109 - t110 * t152;
t78 = t239 * t131 - t132 * t152;
t177 = t206 * qJD(2);
t174 = t153 * t187;
t171 = t55 - t197;
t90 = -pkin(4) * t165 + t142;
t167 = t92 * t96 + t178;
t81 = t153 * t179 + t155 * t195;
t82 = -t153 * t195 + t155 * t179;
t18 = t109 * t183 - t110 * t202 + t152 * t82 + t239 * t81;
t164 = t143 - t197;
t163 = -t123 * t149 - t175;
t162 = -t154 * t215 + t155 * t203;
t161 = -t148 * t192 - t209;
t160 = t211 * t92 + t166;
t19 = -t59 * qJD(4) - t152 * t81 + t239 * t82;
t158 = (-t76 + (-pkin(3) * t148 - t77) * t152) * qJD(4) + t180;
t141 = -pkin(2) - t238;
t140 = t198 + pkin(4);
t133 = t148 * t176;
t126 = t142 - t238;
t121 = -0.2e1 * t174;
t120 = 0.2e1 * t174;
t119 = t144 + t143;
t108 = t165 * qJ(5);
t102 = t123 * t205;
t91 = -0.2e1 * t207 * t149 * qJD(3);
t89 = t92 ^ 2;
t83 = t90 - t238;
t68 = pkin(4) * t94 + t199;
t61 = t65 * t148;
t60 = t64 * t148;
t51 = t108 + t79;
t50 = t78 - t216;
t48 = t144 + t55;
t42 = t108 + t59;
t41 = t58 - t216;
t34 = -t89 + t241;
t31 = t161 + t222;
t28 = -t219 + t40;
t27 = t39 + t220;
t15 = -t165 * t47 + t65 * t92;
t14 = -t113 * t46 - t64 * t94;
t9 = t168 + t19;
t8 = t18 + t217;
t3 = -t113 * t47 - t165 * t46 + t64 * t92 - t65 * t94;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149 * t144 - t137, t201 * t195, 0, 0, t120, t91, t145, t121, -t212, 0, t141 * t188 - t139 * t145 + t102 + (t201 * t213 - t186) * t223, t139 * t212 + t141 * t187 - t162 * t223 + t210, t149 * t177 * t238 + t208, ((qJD(1) * t141 + t123) * t154 + (qJD(1) * t139 + t122) * t156 * t206) * t223, t14, t3, -t60, t15, -t61, 0, t119 * t92 + t126 * t47 + t148 * t19 + t227, t119 * t94 - t126 * t46 - t148 * t18 + t225, -t18 * t92 - t19 * t94 + t46 * t58 - t47 * t59 + t184, t101 * t126 + t11 * t58 + t119 * t96 - t178 * t59 + t18 * t38 + t19 * t37, t14, t3, -t60, t15, -t61, 0, t148 * t9 + t47 * t83 + t48 * t92 + t231, -t148 * t8 - t46 * t83 + t48 * t94 + t230, t41 * t46 - t42 * t47 - t8 * t92 - t9 * t94 + t185, t22 * t9 + t24 * t8 + t33 * t83 + t4 * t42 + t41 * t5 + t48 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149 * t197 - t137, t200 * t196, 0, 0, t120, t91, t145, t121, -t212, 0, -pkin(2) * t188 - pkin(7) * t145 + t102 + (t200 * t213 + t186) * t224, -pkin(2) * t187 + pkin(7) * t212 + t162 * t224 + t210, -t206 * t149 * t196 + t208, ((-pkin(2) * qJD(2) - t123) * t154 + (pkin(7) * t177 - t206 * t122) * t156) * t224, t14, t3, -t60, t15, -t61, 0, t142 * t47 + t228 * t148 + t164 * t92 + t227, -t142 * t46 - t229 * t148 + t164 * t94 + t225, -t228 * t94 - t229 * t92 + t46 * t78 - t47 * t79 + t184, t101 * t142 + t11 * t78 + t164 * t96 - t178 * t79 + t228 * t37 + t229 * t38, t14, t3, -t60, t15, -t61, 0, t233 * t148 + t171 * t92 + t47 * t90 + t231, -t234 * t148 + t171 * t94 - t46 * t90 + t230, -t233 * t94 - t234 * t92 + t46 * t50 - t47 * t51 + t185, t171 * t49 + t233 * t22 + t234 * t24 + t33 * t90 + t4 * t51 + t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t207 * t147, 0, t191, 0, 0, t163 * t153, t163 * t155, 0, 0, t237, t34, t31, -t237, 0, 0, -t39 * t148 - t92 * t199 + t158 - t236, t148 * t40 - t94 * t199 - t133 + t167, t46 * t198 + (-t37 + t40) * t92 + (t38 + t39 + t194) * t94 + t226, -t37 * t39 - t38 * t40 + (-t96 * t215 + t239 * t11 - t178 * t152 + (-t152 * t37 + t239 * t38) * qJD(4)) * pkin(3), t237, t34, t31, -t237, 0, 0, -t27 * t148 - t211 * t94 - t68 * t92 + t158 + t221, t148 * t28 - t68 * t94 - t133 + t160, t140 * t46 + (-t22 + t28) * t92 + (t24 + t27 + t194) * t94 + t226, t5 * t140 - t22 * t27 - t24 * t28 - t49 * t68 + (t152 * t4 + (-t152 * t22 + t239 * t24) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t34, t31, -t237, 0, 0, t38 * t148 + t11 - t236, t148 * t37 + t167, 0, 0, t237, t34, t31, -t237, 0, 0, t24 * t148 + (-t182 - t49) * t94 + t159, -pkin(4) * t241 + t148 * t23 + t160, pkin(4) * t46 - t232 * t92, t232 * t24 + (-t49 * t94 + t5) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t148 + t47, t161 - t222, -t89 - t241, t22 * t94 + t24 * t92 + t33;];
tauc_reg = t1;
