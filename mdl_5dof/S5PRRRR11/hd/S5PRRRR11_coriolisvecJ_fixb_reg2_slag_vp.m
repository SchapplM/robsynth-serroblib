% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:46:00
% DurationCPUTime: 1.98s
% Computational Cost: add. (4558->293), mult. (14684->421), div. (0->0), fcn. (11294->8), ass. (0->185)
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t207 = qJD(4) * t155;
t208 = qJD(4) * t152;
t150 = cos(pkin(5));
t204 = t150 * qJD(2);
t140 = qJD(3) + t204;
t153 = sin(qJ(3));
t149 = sin(pkin(5));
t245 = pkin(7) + pkin(8);
t200 = t149 * t245;
t181 = t153 * t200;
t156 = cos(qJ(3));
t244 = pkin(2) * t150;
t143 = t156 * t244;
t213 = qJD(1) * t149;
t216 = qJD(2) * t143 + t156 * t213;
t87 = -qJD(2) * t181 + t216;
t77 = t140 * pkin(3) + t87;
t172 = qJD(3) * t181;
t202 = qJD(2) * qJD(3);
t188 = t156 * t202;
t175 = t150 * t188;
t210 = qJD(3) * t149;
t189 = qJD(1) * t210;
t217 = pkin(2) * t175 + t156 * t189;
t79 = -qJD(2) * t172 + t217;
t142 = t153 * t244;
t162 = -t156 * t200 - t142;
t192 = t153 * t213;
t80 = (t162 * qJD(2) - t192) * qJD(3);
t221 = t149 * t156;
t88 = t192 + (t245 * t221 + t142) * qJD(2);
t184 = -t152 * t80 - t155 * t79 - t77 * t207 + t88 * t208;
t211 = qJD(2) * t149;
t170 = t152 * t156 + t153 * t155;
t246 = t170 * qJD(4);
t248 = t170 * qJD(3) + t246;
t65 = t248 * t211;
t8 = -t65 * pkin(9) - t184;
t185 = -t152 * t79 + t155 * t80;
t84 = t155 * t88;
t40 = t152 * t77 + t84;
t20 = -t40 * qJD(4) + t185;
t222 = t149 * t153;
t199 = t152 * t222;
t201 = qJD(3) + qJD(4);
t168 = t201 * t199;
t195 = t156 * t211;
t177 = t155 * t195;
t218 = t201 * t177;
t64 = qJD(2) * t168 - t218;
t9 = t64 * pkin(9) + t20;
t197 = -t151 * t8 + t154 * t9;
t196 = t153 * t211;
t178 = t152 * t196;
t112 = -t177 + t178;
t243 = t112 * pkin(9);
t36 = t40 - t243;
t228 = t154 * t36;
t136 = qJD(4) + t140;
t117 = t170 * t149;
t114 = qJD(2) * t117;
t109 = t114 * pkin(9);
t82 = t152 * t88;
t39 = t155 * t77 - t82;
t35 = -t109 + t39;
t33 = t136 * pkin(4) + t35;
t6 = t151 * t33 + t228;
t2 = -t6 * qJD(5) + t197;
t171 = t151 * t112 - t154 * t114;
t145 = t150 * qJD(1);
t124 = (-pkin(3) * t156 - pkin(2)) * t149;
t212 = qJD(2) * t124;
t115 = t145 + t212;
t73 = t112 * pkin(4) + t115;
t237 = t73 * t171;
t251 = t2 + t237;
t206 = qJD(5) * t151;
t190 = t151 * t9 - t36 * t206;
t1 = (qJD(5) * t33 + t8) * t154 + t190;
t60 = t154 * t112 + t151 * t114;
t236 = t73 * t60;
t250 = t236 - t1;
t249 = t60 * t171;
t16 = t171 ^ 2 - t60 ^ 2;
t132 = qJD(5) + t136;
t205 = qJD(5) * t154;
t26 = t112 * t205 + t114 * t206 + t151 * t65 + t154 * t64;
t12 = t60 * t132 - t26;
t160 = t171 * qJD(5) + t151 * t64 - t154 * t65;
t13 = -t132 * t171 + t160;
t215 = pkin(7) * t221 + t142;
t110 = pkin(8) * t221 + t215;
t99 = t150 * pkin(3) + t143 - t181;
t52 = t155 * t110 + t152 * t99;
t104 = qJD(2) * t215 + t192;
t242 = t114 * pkin(4);
t241 = t153 * pkin(3);
t144 = t155 * pkin(3) + pkin(4);
t219 = t152 * t154;
t45 = -t152 * t87 - t84;
t37 = t45 + t243;
t46 = t155 * t87 - t82;
t38 = -t109 + t46;
t235 = -t151 * t38 + t154 * t37 + t144 * t206 - (-t152 * t205 + (-t151 * t155 - t219) * qJD(4)) * pkin(3);
t220 = t151 * t152;
t234 = t151 * t37 + t154 * t38 - t144 * t205 - (-t152 * t206 + (t154 * t155 - t220) * qJD(4)) * pkin(3);
t198 = t155 * t221;
t116 = -t198 + t199;
t159 = t149 * t248;
t193 = t156 * t210;
t78 = -qJD(4) * t198 - t155 * t193 + t168;
t29 = t116 * t205 + t117 * t206 + t151 * t159 + t154 * t78;
t72 = -t151 * t116 + t154 * t117;
t233 = t160 * t72 + t29 * t60;
t232 = t78 * t112 - t117 * t65;
t231 = t150 * t160;
t230 = t150 * t65;
t229 = t151 * t36;
t227 = t29 * t132;
t226 = t78 * t136;
t225 = t114 * t112;
t224 = t115 * t114;
t146 = t149 ^ 2;
t223 = t146 * qJD(2) ^ 2;
t214 = t153 ^ 2 - t156 ^ 2;
t209 = qJD(3) * t153;
t203 = qJD(3) - t140;
t194 = t149 * t209;
t191 = -pkin(4) * t132 - t33;
t51 = -t152 * t110 + t155 * t99;
t183 = pkin(3) * t194;
t182 = pkin(3) * t196;
t180 = t153 * t156 * t223;
t179 = t140 * t193;
t176 = qJD(2) * t194;
t5 = t154 * t33 - t229;
t174 = -t171 * t6 - t5 * t60;
t30 = t72 * qJD(5) - t151 * t78 + t154 * t159;
t71 = t154 * t116 + t151 * t117;
t173 = -t171 * t30 - t71 * t26;
t41 = t150 * pkin(4) - t117 * pkin(9) + t51;
t44 = -t116 * pkin(9) + t52;
t21 = -t151 * t44 + t154 * t41;
t22 = t151 * t41 + t154 * t44;
t169 = t146 * t153 * t188;
t167 = t115 * t112 + t184;
t138 = qJD(3) * t143;
t100 = t138 - t172;
t101 = t162 * qJD(3);
t31 = t155 * t100 + t152 * t101 - t110 * t208 + t99 * t207;
t126 = -pkin(2) * t211 + t145;
t163 = qJD(3) * (-pkin(2) * qJD(2) * t146 + t126 * t149);
t119 = t215 * qJD(3);
t32 = -t52 * qJD(4) - t152 * t100 + t155 * t101;
t158 = t114 * t159 - t116 * t64;
t125 = t149 * t175;
t122 = -pkin(7) * t222 + t143;
t121 = pkin(3) * t219 + t151 * t144;
t120 = -pkin(3) * t220 + t154 * t144;
t118 = -pkin(7) * t194 + t138;
t103 = -pkin(7) * t196 + t216;
t95 = t104 * qJD(3);
t94 = -pkin(7) * t176 + t217;
t90 = t182 + t242;
t89 = t116 * pkin(4) + t124;
t55 = t136 * t159;
t54 = t64 * t150;
t53 = (pkin(4) * t246 + (t170 * pkin(4) + t241) * qJD(3)) * t149;
t50 = pkin(3) * t176 + t65 * pkin(4);
t48 = -t112 ^ 2 + t114 ^ 2;
t43 = -qJD(2) * t159 + t114 * t136;
t42 = t112 * t136 - t201 * t178 + t218;
t28 = t30 * t132;
t25 = t26 * t150;
t24 = t78 * pkin(9) + t32;
t23 = -pkin(9) * t159 + t31;
t11 = t154 * t35 - t229;
t10 = -t151 * t35 - t228;
t4 = -t22 * qJD(5) - t151 * t23 + t154 * t24;
t3 = t21 * qJD(5) + t151 * t24 + t154 * t23;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t140 + t204) * t194, t125 - t179, 0, (t153 * t94 - t156 * t95 + (-t103 * t153 + t104 * t156) * qJD(3)) * t149, 0, 0, 0, 0, 0, 0, -t55 + t230, -t54 + t226, t158 + t232, -t20 * t116 - t184 * t117 - t40 * t78 + (-t39 * t246 + (-t39 * t170 + t204 * t241) * qJD(3)) * t149, 0, 0, 0, 0, 0, 0, -t28 - t231, -t25 + t227, t173 + t233, t1 * t72 + t50 * t150 - t2 * t71 - t6 * t29 - t5 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t169, -0.2e1 * t214 * t146 * t202, t125 + t179, -0.2e1 * t169, (-t140 - t204) * t194, 0, -t119 * t140 - t95 * t150 + t153 * t163, -t118 * t140 - t94 * t150 + t156 * t163, (t153 * t95 + t156 * t94 + (-t103 * t156 - t104 * t153) * qJD(3) + (t118 * t156 + t119 * t153 + (-t122 * t156 - t153 * t215) * qJD(3)) * qJD(2)) * t149, -t103 * t119 + t104 * t118 - t95 * t122 + t215 * t94, -t114 * t78 - t64 * t117, -t158 + t232, -t54 - t226, t112 * t159 + t65 * t116, -t55 - t230, 0, t124 * t65 + t32 * t136 + t20 * t150 + (t115 * t246 + (t115 * t170 + (qJD(2) * t116 + t112) * t241) * qJD(3)) * t149, 0.2e1 * t114 * t183 - t115 * t78 - t124 * t64 - t31 * t136 + t150 * t184, -t31 * t112 - t32 * t114 + t116 * t184 - t20 * t117 - t40 * t159 + t39 * t78 + t51 * t64 - t52 * t65, -t184 * t52 + t20 * t51 + t40 * t31 + t39 * t32 + (t115 + t212) * t183, t171 * t29 - t26 * t72, -t173 + t233, -t25 - t227, -t160 * t71 + t60 * t30, -t28 + t231, 0, t4 * t132 + t2 * t150 - t160 * t89 + t73 * t30 + t50 * t71 + t53 * t60, -t1 * t150 - t3 * t132 - t171 * t53 - t89 * t26 - t73 * t29 + t50 * t72, -t1 * t71 + t160 * t22 + t171 * t4 - t2 * t72 + t21 * t26 + t5 * t29 - t3 * t60 - t6 * t30, t1 * t22 + t2 * t21 + t6 * t3 + t5 * t4 + t50 * t89 + t73 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t214 * t223, t203 * t195, t180, -t203 * t196, 0, -t153 * t189 + t104 * t140 + (-t126 * t222 - t119) * qJD(2), t103 * t140 + (pkin(7) * t209 - t126 * t156) * t211 - t217, 0, 0, t225, t48, t42, -t225, t43, 0, -t112 * t182 - t224 - t45 * t136 + (-t84 + (-pkin(3) * t136 - t77) * t152) * qJD(4) + t185, t46 * t136 + (-t114 * t196 - t136 * t207) * pkin(3) + t167, (t40 + t45) * t114 + (-t39 + t46) * t112 + (-t152 * t65 + t155 * t64 + (-t112 * t155 + t114 * t152) * qJD(4)) * pkin(3), -t39 * t45 - t40 * t46 + (-t115 * t196 - t152 * t184 + t155 * t20 + (-t152 * t39 + t155 * t40) * qJD(4)) * pkin(3), -t249, t16, t12, t249, t13, 0, -t235 * t132 - t90 * t60 + t251, t234 * t132 + t171 * t90 + t250, t120 * t26 + t121 * t160 - t171 * t235 + t234 * t60 + t174, t1 * t121 + t2 * t120 - t234 * t6 - t235 * t5 - t73 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, t48, t42, -t225, t43, 0, t40 * t136 + t20 - t224, t39 * t136 + t167, 0, 0, -t249, t16, t12, t249, t13, 0, -t60 * t242 - t10 * t132 + t237 + (t191 * t151 - t228) * qJD(5) + t197, t171 * t242 + t11 * t132 + t236 + (t191 * qJD(5) - t8) * t154 - t190, -t10 * t171 + t11 * t60 + (t151 * t160 + t154 * t26 + (-t151 * t171 - t154 * t60) * qJD(5)) * pkin(4) + t174, -t5 * t10 - t6 * t11 + (t1 * t151 - t114 * t73 + t154 * t2 + (-t151 * t5 + t154 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249, t16, t12, t249, t13, 0, t6 * t132 + t251, t5 * t132 + t250, 0, 0;];
tauc_reg = t7;
