% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:17
% EndTime: 2019-03-09 02:01:24
% DurationCPUTime: 2.29s
% Computational Cost: add. (5288->339), mult. (12658->416), div. (0->0), fcn. (9151->8), ass. (0->174)
t133 = cos(qJ(5));
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t132 = sin(qJ(4));
t222 = cos(qJ(4));
t112 = t222 * t127 + t132 * t129;
t104 = t112 * qJD(1);
t187 = t132 * t127;
t173 = qJD(4) * t187;
t116 = qJD(1) * t173;
t174 = t222 * t129;
t117 = qJD(1) * t174;
t146 = -qJD(4) * t117 + t116;
t181 = t133 * qJD(4);
t131 = sin(qJ(5));
t184 = qJD(5) * t131;
t55 = -qJD(5) * t181 + t104 * t184 + t133 * t146;
t87 = qJD(4) * t131 + t104 * t133;
t80 = t87 * t184;
t237 = -t133 * t55 - t80;
t183 = qJD(5) * t133;
t149 = t174 - t187;
t142 = t149 * qJD(3);
t232 = qJD(1) * t142;
t119 = sin(pkin(9)) * pkin(1) + qJ(3);
t115 = t119 * qJD(1);
t123 = t129 * qJD(2);
t208 = pkin(7) * qJD(1);
t83 = t123 + (-t115 - t208) * t127;
t94 = t127 * qJD(2) + t129 * t115;
t84 = t129 * t208 + t94;
t235 = -t132 * t84 + t222 * t83;
t32 = qJD(4) * t235 + t232;
t46 = t132 * t83 + t222 * t84;
t43 = qJD(4) * pkin(8) + t46;
t102 = qJD(1) * t187 - t117;
t113 = -cos(pkin(9)) * pkin(1) - pkin(3) * t129 - pkin(2);
t99 = t113 * qJD(1) + qJD(3);
t54 = t102 * pkin(4) - t104 * pkin(8) + t99;
t107 = t112 * qJD(4);
t95 = qJD(1) * t107;
t64 = t95 * pkin(4) + pkin(8) * t146;
t171 = t131 * t32 - t133 * t64 + t43 * t183 + t54 * t184;
t18 = t131 * t54 + t133 * t43;
t98 = qJD(5) + t102;
t219 = t18 * t98;
t240 = -t171 + t219;
t225 = pkin(5) * t95;
t2 = t171 - t225;
t14 = qJ(6) * t98 + t18;
t221 = t14 * t98;
t242 = -t2 + t221;
t17 = -t131 * t43 + t133 * t54;
t3 = t131 * t64 + t133 * t32 + t54 * t183 - t43 * t184;
t241 = -t17 * t98 + t3;
t90 = t131 * t95;
t151 = t98 * t183 + t90;
t92 = t133 * t95;
t239 = t98 * t184 - t92;
t56 = qJD(5) * t87 - t131 * t146;
t238 = -t87 * t98 + t56;
t85 = t104 * t131 - t181;
t203 = t133 * t85;
t205 = t131 * t87;
t213 = -t131 * t56 - t85 * t183;
t236 = (-t203 + t205) * t102 + t213 - t237;
t215 = pkin(7) + t119;
t108 = t215 * t127;
t109 = t215 * t129;
t234 = -t222 * t108 - t132 * t109;
t106 = -qJD(4) * t174 + t173;
t191 = t106 * t131;
t139 = -t112 * t151 + t98 * t191;
t214 = t107 * t85 - t149 * t56;
t233 = t139 - t214;
t231 = t104 * qJD(4);
t42 = -qJD(4) * pkin(4) - t235;
t21 = t85 * pkin(5) - t87 * qJ(6) + t42;
t223 = pkin(8) * t95;
t230 = t21 * t98 - t223;
t229 = -t104 * t107 - t146 * t149;
t68 = -pkin(4) * t149 - pkin(8) * t112 + t113;
t72 = -t132 * t108 + t222 * t109;
t210 = t131 * t68 + t133 * t72;
t49 = t234 * qJD(4) + t142;
t76 = pkin(4) * t107 + pkin(8) * t106;
t10 = -qJD(5) * t210 - t131 * t49 + t133 * t76;
t228 = t87 ^ 2;
t227 = t98 ^ 2;
t226 = t104 ^ 2;
t224 = pkin(8) * t87;
t143 = t112 * qJD(3);
t33 = qJD(1) * t143 + t46 * qJD(4);
t218 = t33 * t234;
t169 = t85 * t98;
t217 = t87 * t85;
t75 = pkin(4) * t104 + pkin(8) * t102;
t23 = t131 * t75 + t133 * t235;
t189 = t112 * t133;
t190 = t106 * t133;
t60 = t85 * t190;
t212 = -t56 * t189 + t60;
t211 = t87 * t107 + t149 * t55;
t209 = t106 * t102 - t112 * t95;
t207 = t104 * t85;
t206 = t149 * t95;
t202 = t133 * t98;
t201 = t33 * t149;
t200 = t55 * t131;
t199 = t56 * t133;
t198 = t87 * t104;
t197 = t95 * qJ(6);
t196 = t98 * t104;
t163 = pkin(5) * t131 - qJ(6) * t133;
t195 = -t131 * qJD(6) + t98 * t163 - t46;
t193 = t104 * t102;
t186 = qJD(6) - t17;
t185 = t127 ^ 2 + t129 ^ 2;
t182 = t106 * qJD(4);
t61 = t87 * t191;
t180 = -t61 + (t183 * t87 - t200) * t112;
t179 = pkin(8) * qJD(5) * t98;
t74 = t98 * t190;
t177 = t85 ^ 2 - t228;
t172 = t112 * t184;
t168 = qJD(1) * t185;
t7 = t56 * pkin(5) + t55 * qJ(6) - t87 * qJD(6) + t33;
t167 = -t7 - t179;
t166 = t33 + t179;
t165 = (qJD(5) * t85 - t55) * pkin(8);
t164 = pkin(5) * t133 + qJ(6) * t131;
t162 = t127 * (-t115 * t127 + t123) - t129 * t94;
t13 = -pkin(5) * t98 + t186;
t161 = t13 * t133 - t131 * t14;
t160 = -t13 * t131 - t133 * t14;
t159 = t131 * t18 + t133 * t17;
t158 = t131 * t17 - t133 * t18;
t22 = -t131 * t235 + t133 * t75;
t28 = -t131 * t72 + t133 * t68;
t155 = t102 * t202 + t151;
t154 = -t131 * t102 * t98 - t239;
t153 = t21 * t87 + t171;
t9 = t131 * t76 + t133 * t49 + t68 * t183 - t72 * t184;
t150 = t98 * t42 - t223;
t148 = t85 * t172 + t212;
t147 = -t98 * t172 + t95 * t189 - t74;
t145 = t131 * t169 - t199;
t144 = -t154 - t207;
t138 = -t213 * t112 - t85 * t191;
t1 = qJD(6) * t98 + t197 + t3;
t137 = t161 * qJD(5) + t1 * t133 + t2 * t131;
t136 = -t159 * qJD(5) + t131 * t171 + t3 * t133;
t135 = t139 + t214;
t50 = qJD(4) * t72 + t143;
t114 = -pkin(4) - t164;
t101 = t102 ^ 2;
t97 = t107 * qJD(4);
t52 = pkin(5) * t87 + qJ(6) * t85;
t51 = pkin(8) * t199;
t35 = t107 * t98 - t206;
t34 = t163 * t112 - t234;
t27 = -t55 + t169;
t25 = pkin(5) * t149 - t28;
t24 = -qJ(6) * t149 + t210;
t20 = -pkin(5) * t104 - t22;
t19 = qJ(6) * t104 + t23;
t16 = t155 - t198;
t15 = t87 * t202 - t200;
t12 = -t163 * t106 + (t164 * qJD(5) - qJD(6) * t133) * t112 + t50;
t11 = t237 * t112 - t87 * t190;
t8 = -t107 * pkin(5) - t10;
t6 = t147 + t211;
t5 = qJ(6) * t107 - qJD(6) * t149 + t9;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t168 (t119 * t168 - t162) * qJD(3), -t104 * t106 - t112 * t146, t209 + t229, -t182, t102 * t107 - t206, -t97, 0, -qJD(4) * t50 + t107 * t99 + t113 * t95, -t49 * qJD(4) - t99 * t106 - t113 * t146, -t49 * t102 + t50 * t104 + t106 * t235 - t46 * t107 + t33 * t112 + t146 * t234 + t149 * t32 - t72 * t95, -t235 * t50 + t32 * t72 + t46 * t49 - t218, t11, t148 - t180, t6, t138, t233, t35, -t42 * t191 + t10 * t98 + t17 * t107 + t171 * t149 + t28 * t95 + t50 * t85 - t234 * t56 + (t33 * t131 + t42 * t183) * t112, -t42 * t190 - t18 * t107 + t3 * t149 - t210 * t95 + t50 * t87 + t234 * t55 - t9 * t98 + (t33 * t133 - t42 * t184) * t112, -t10 * t87 + t28 * t55 - t210 * t56 - t9 * t85 + t159 * t106 + (t158 * qJD(5) - t3 * t131 + t133 * t171) * t112, t10 * t17 - t171 * t28 + t18 * t9 + t210 * t3 + t42 * t50 - t218, t11, t6, -t60 + (-t184 * t85 + t199) * t112 + t180, t35, -t233, t138, -t21 * t191 - t13 * t107 + t2 * t149 + t12 * t85 - t25 * t95 + t34 * t56 - t8 * t98 + (t7 * t131 + t183 * t21) * t112, -t24 * t56 - t25 * t55 - t5 * t85 + t8 * t87 - t161 * t106 + (qJD(5) * t160 - t1 * t131 + t2 * t133) * t112, t21 * t190 - t1 * t149 + t14 * t107 - t12 * t87 + t24 * t95 + t34 * t55 + t5 * t98 + (-t7 * t133 + t184 * t21) * t112, t1 * t24 + t12 * t21 + t13 * t8 + t14 * t5 + t2 * t25 + t34 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t182, t209 - t229, -t106 * t46 - t107 * t235 + t112 * t32 - t201, 0, 0, 0, 0, 0, 0, t135, t239 * t112 + t211 + t74, -t61 + (-t200 + (t131 * t85 + t133 * t87) * qJD(5)) * t112 + t212, t106 * t158 + t42 * t107 + t112 * t136 - t201, 0, 0, 0, 0, 0, 0, t135, t148 + t180, t147 - t211, t106 * t160 + t21 * t107 + t112 * t137 - t149 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 * qJD(1) ^ 2, t162 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t231, -t116 + (t117 - t102) * qJD(4), -t101 - t226, t102 * t46 + t104 * t235, 0, 0, 0, 0, 0, 0, t154 - t207, -t133 * t227 - t198 - t90, t236, -t42 * t104 + t241 * t131 + t240 * t133, 0, 0, 0, 0, 0, 0, -t131 * t227 - t207 + t92, t236, t155 + t198, -t21 * t104 + t242 * t133 + (t13 * t98 + t1) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t101 + t226, -t116 + (t117 + t102) * qJD(4), -t193, 0, 0 -(qJD(3) + t99) * t104, t99 * t102 - t232, 0, 0, t15 (-t203 - t205) * t102 + t213 + t237, t16, t145, -t144, -t196, -pkin(4) * t56 - t17 * t104 + t150 * t131 - t166 * t133 - t22 * t98 - t46 * t85, pkin(4) * t55 + t18 * t104 + t166 * t131 + t150 * t133 + t23 * t98 - t46 * t87, t22 * t87 + t23 * t85 - t51 + (-t102 * t17 + t3 + (-t17 + t224) * qJD(5)) * t133 + (t165 - t240) * t131, -t33 * pkin(4) + pkin(8) * t136 - t17 * t22 - t18 * t23 - t42 * t46, t15, t16, t80 + (t102 * t87 + t56) * t131 + (t55 + t169) * t133, -t196, t144, t145, t13 * t104 + t114 * t56 + t230 * t131 + t167 * t133 + t195 * t85 + t20 * t98, t19 * t85 - t20 * t87 - t51 + (t102 * t13 + t1 + (t13 + t224) * qJD(5)) * t133 + (t165 - t242) * t131, -t14 * t104 + t114 * t55 + t167 * t131 - t230 * t133 - t19 * t98 - t195 * t87, pkin(8) * t137 + t7 * t114 - t13 * t20 - t14 * t19 + t195 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, -t177, t27, -t217, -t238, t95, -t42 * t87 + t240, t42 * t85 - t241, 0, 0, t217, t27, t177, t95, t238, -t217, -t52 * t85 - t153 + t219 + 0.2e1 * t225, pkin(5) * t55 - t56 * qJ(6) + (t14 - t18) * t87 + (t13 - t186) * t85, 0.2e1 * t197 - t21 * t85 + t52 * t87 + (0.2e1 * qJD(6) - t17) * t98 + t3, -t2 * pkin(5) + t1 * qJ(6) - t13 * t18 + t14 * t186 - t21 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217 - t231, t27, -t227 - t228, t153 - t221 - t225;];
tauc_reg  = t4;
