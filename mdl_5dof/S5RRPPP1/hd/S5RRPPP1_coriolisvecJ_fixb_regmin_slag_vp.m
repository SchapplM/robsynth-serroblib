% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:30
% EndTime: 2019-12-31 19:24:38
% DurationCPUTime: 2.08s
% Computational Cost: add. (2730->330), mult. (8465->462), div. (0->0), fcn. (6005->6), ass. (0->163)
t176 = (qJD(1) * qJD(2));
t216 = -2 * t176;
t136 = sin(pkin(5));
t138 = cos(pkin(5));
t137 = cos(pkin(8));
t179 = qJD(3) * t137;
t139 = sin(qJ(2));
t186 = qJD(1) * t139;
t140 = cos(qJ(2));
t195 = t136 * t140;
t155 = pkin(2) * t139 - qJ(3) * t195;
t106 = t155 * qJD(1);
t168 = qJ(3) * t138 + pkin(7);
t158 = qJD(1) * t168;
t107 = t139 * t158;
t108 = t140 * t158;
t135 = sin(pkin(8));
t198 = t135 * t138;
t199 = t135 * t136;
t37 = t106 * t199 - t137 * t107 - t108 * t198;
t203 = -qJD(4) * t138 + t37 + (qJ(4) * t186 - t179) * t136;
t177 = t138 * qJD(2);
t185 = qJD(1) * t140;
t112 = t136 * t185 - t177;
t109 = t112 ^ 2;
t192 = t138 * t140;
t193 = t137 * t139;
t100 = t135 * t192 + t193;
t184 = qJD(2) * t136;
t64 = t100 * qJD(1) + t135 * t184;
t215 = -t64 ^ 2 - t109;
t115 = -qJ(3) * t136 * t139 - pkin(2) * t140 - pkin(1);
t116 = t168 * t139;
t214 = (t115 * t136 - t116 * t138) * t137;
t213 = t136 ^ 2 + t138 ^ 2;
t182 = qJD(2) * t139;
t178 = qJD(3) * t139;
t74 = t155 * qJD(2) - t136 * t178;
t157 = qJD(2) * t168;
t75 = qJD(3) * t192 - t139 * t157;
t76 = -t138 * t178 - t140 * t157;
t27 = t137 * t75 + t76 * t198 + t74 * t199;
t17 = -(qJ(4) * t182 - qJD(4) * t140) * t136 - t27;
t173 = t138 * t185;
t114 = t173 + t184;
t82 = pkin(7) * t185 + t114 * qJ(3);
t95 = qJD(2) * pkin(2) - t107;
t96 = t115 * qJD(1);
t30 = -t135 * t82 + (t136 * t96 + t138 * t95) * t137;
t174 = t135 * t186;
t183 = qJD(2) * t137;
t63 = -t136 * t183 - t137 * t173 + t174;
t210 = t63 * t64;
t209 = pkin(3) + qJ(5);
t46 = t138 * t106 + t108 * t136;
t165 = t138 * t174;
t88 = t137 * t185 - t165;
t154 = -qJ(4) * t88 + t46;
t152 = t135 * t140 + t138 * t193;
t87 = t152 * qJD(1);
t208 = t209 * t87 + t154 - (-qJD(4) * t135 - qJD(5) * t137) * t136;
t61 = t74 * qJD(1);
t207 = t137 * t61;
t206 = t137 * t74;
t44 = t213 * t137 * t186 + t135 * t114;
t205 = t44 * t112;
t194 = t137 * t138;
t92 = t135 * t107;
t161 = t108 * t194 - t92;
t172 = t209 * t139;
t180 = qJD(3) * t136;
t200 = t106 * t137;
t204 = -qJD(5) * t138 + t135 * t180 - pkin(4) * t88 - (-qJD(1) * t172 - t200) * t136 - t161;
t202 = pkin(4) * t87 - t203;
t201 = qJD(4) * t64;
t197 = t135 * t139;
t196 = t136 * t137;
t142 = qJD(1) ^ 2;
t191 = t140 * t142;
t141 = qJD(2) ^ 2;
t190 = t141 * t139;
t189 = t141 * t140;
t117 = t168 * t140;
t103 = t135 * t117;
t188 = pkin(3) * t195 + t103;
t102 = pkin(2) * t198 + qJ(3) * t196;
t187 = t139 ^ 2 - t140 ^ 2;
t181 = qJD(3) * t112;
t52 = t75 * qJD(1) + qJD(2) * t180;
t62 = t76 * qJD(1);
t15 = t137 * t52 + t62 * t198 + t61 * t199;
t31 = t137 * t82 + t95 * t198 + t96 * t199;
t41 = t115 * t199 - t116 * t198 + t137 * t117;
t175 = -pkin(2) * t137 - pkin(3);
t171 = t139 * t176;
t170 = t140 * t176;
t169 = -qJ(4) * t135 - pkin(2);
t126 = t136 * t171;
t167 = -t126 + t210;
t39 = -t136 * t62 + t138 * t61;
t42 = -t136 * t76 + t138 * t74;
t49 = t138 * t115 + t116 * t136;
t166 = pkin(1) * t216;
t47 = t135 * t52;
t164 = -t62 * t194 + t47;
t66 = t135 * t75;
t163 = -t76 * t194 + t66;
t162 = qJD(2) * t172;
t43 = -t136 * t95 + t138 * t96 + qJD(3);
t9 = -qJ(4) * t126 + t112 * qJD(4) - t15;
t45 = t114 * t137 - t213 * t174;
t160 = -t44 * t64 + t45 * t63;
t83 = -qJ(4) * t138 - t102;
t23 = qJ(4) * t112 - t31;
t78 = -qJD(2) * t165 + t137 * t170;
t153 = -qJ(4) * t100 + t49;
t35 = qJ(4) * t195 - t41;
t89 = t152 * qJD(2);
t77 = qJD(1) * t89;
t3 = -pkin(4) * t77 - t9;
t150 = -qJ(4) * t64 + t43;
t149 = -t112 * t45 + t78;
t148 = -t112 * t63 + t78;
t90 = t140 * t183 - t177 * t197;
t147 = -qJ(4) * t90 - qJD(4) * t100 + t42;
t146 = t152 * t176;
t145 = qJD(4) - t30;
t5 = pkin(3) * t77 - qJ(4) * t78 - t201 + t39;
t144 = t77 - t205;
t1 = t77 * qJ(5) + t63 * qJD(5) + t5;
t12 = (-pkin(3) * t171 - t207) * t136 + t164;
t143 = pkin(4) * t78 + (-qJD(1) * t162 - t207) * t136 + t164;
t128 = qJ(3) * t199;
t101 = pkin(2) * t194 - t128;
t99 = -t137 * t192 + t197;
t85 = (-pkin(3) * t137 + t169) * t136;
t84 = t175 * t138 + t128;
t59 = (-t209 * t137 + t169) * t136;
t58 = pkin(4) * t196 - t83;
t53 = pkin(4) * t199 + t128 + (-qJ(5) + t175) * t138;
t40 = -t103 + t214;
t38 = t188 - t214;
t36 = t92 + (t106 * t136 - t108 * t138) * t137;
t34 = pkin(3) * t99 + t153;
t33 = (-pkin(3) * t186 - t200) * t136 + t161;
t29 = pkin(3) * t87 + t154;
t28 = -pkin(4) * t99 - t35;
t26 = -t66 + (t136 * t74 + t138 * t76) * t137;
t25 = t209 * t99 + t153;
t24 = t116 * t194 + pkin(4) * t100 + (qJ(5) * t140 - t115 * t137) * t136 + t188;
t22 = pkin(3) * t112 + t145;
t20 = (-pkin(3) * t182 - t206) * t136 + t163;
t16 = pkin(3) * t63 + t150;
t14 = -t47 + (t136 * t61 + t138 * t62) * t137;
t13 = pkin(3) * t89 + t147;
t11 = -pkin(4) * t63 + qJD(5) - t23;
t10 = -pkin(4) * t89 - t17;
t8 = pkin(4) * t90 + (qJD(5) * t140 - t162 - t206) * t136 + t163;
t7 = t209 * t63 + t150;
t6 = pkin(4) * t64 + t209 * t112 + t145;
t4 = qJD(5) * t99 + t209 * t89 + t147;
t2 = qJD(5) * t112 + t143;
t18 = [0, 0, 0, 0.2e1 * t139 * t170, t187 * t216, t189, -t190, 0, -pkin(7) * t189 + t139 * t166, pkin(7) * t190 + t140 * t166, -t112 * t26 + t39 * t99 + t42 * t63 + t43 * t89 + t49 * t77 + (-t14 * t140 + (qJD(1) * t40 + t30) * t182) * t136, t100 * t39 + t112 * t27 + t42 * t64 + t43 * t90 + t49 * t78 + (t140 * t15 + (-qJD(1) * t41 - t31) * t182) * t136, -t100 * t14 - t15 * t99 - t26 * t64 - t27 * t63 - t30 * t90 - t31 * t89 - t40 * t78 - t41 * t77, t14 * t40 + t15 * t41 + t26 * t30 + t27 * t31 + t39 * t49 + t42 * t43, t100 * t12 + t17 * t63 + t20 * t64 + t22 * t90 + t23 * t89 + t35 * t77 + t38 * t78 + t9 * t99, -t112 * t20 - t13 * t63 - t16 * t89 - t34 * t77 - t5 * t99 + (-t12 * t140 + (qJD(1) * t38 + t22) * t182) * t136, -t100 * t5 + t112 * t17 - t13 * t64 - t16 * t90 - t34 * t78 + (t140 * t9 + (-qJD(1) * t35 - t23) * t182) * t136, t12 * t38 + t13 * t16 + t17 * t23 + t20 * t22 + t34 * t5 + t35 * t9, -t10 * t63 + t100 * t2 - t11 * t89 + t24 * t78 - t28 * t77 - t3 * t99 + t6 * t90 + t64 * t8, -t1 * t100 - t10 * t112 - t25 * t78 - t4 * t64 - t7 * t90 + (-t140 * t3 + (qJD(1) * t28 + t11) * t182) * t136, t1 * t99 + t112 * t8 + t25 * t77 + t4 * t63 + t7 * t89 + (t140 * t2 + (-qJD(1) * t24 - t6) * t182) * t136, t1 * t25 + t10 * t11 + t2 * t24 + t28 * t3 + t4 * t7 + t6 * t8; 0, 0, 0, -t139 * t191, t187 * t142, 0, 0, 0, t142 * pkin(1) * t139, pkin(1) * t191, t112 * t36 + t138 * t14 - t43 * t87 - t46 * t63 + (t135 * t181 - pkin(2) * t77 - t137 * t39 + (qJD(2) * t101 - t30) * t186) * t136, -t112 * t37 - t138 * t15 - t43 * t88 - t46 * t64 + (t112 * t179 - pkin(2) * t78 + t135 * t39 + (-qJD(2) * t102 + t31) * t186) * t136, -t101 * t78 - t102 * t77 + t30 * t88 + t31 * t87 + t36 * t64 + t37 * t63 + (-t135 * t14 + t137 * t15 + (t135 * t64 - t137 * t63) * qJD(3)) * t136, t101 * t14 + t102 * t15 - t30 * t36 - t31 * t37 - t43 * t46 + (-pkin(2) * t39 + (-t135 * t30 + t137 * t31) * qJD(3)) * t136, -t22 * t88 - t23 * t87 - t33 * t64 + t77 * t83 + t78 * t84 + t203 * t63 + (-t137 * t9 + (qJD(3) * t64 + t12) * t135) * t136, t112 * t33 + t12 * t138 + t16 * t87 + t29 * t63 - t77 * t85 + (t137 * t5 + (qJD(4) * t63 - t181) * t135 + (qJD(2) * t84 - t22) * t186) * t136, -t138 * t9 + t16 * t88 + t29 * t64 - t78 * t85 + t203 * t112 + ((-t5 + t201) * t135 + (-qJD(2) * t83 + t23) * t186) * t136, t12 * t84 - t16 * t29 - t22 * t33 + t5 * t85 + t83 * t9 + t203 * t23 + (qJD(3) * t22 - qJD(4) * t16) * t199, t11 * t87 + t53 * t78 - t58 * t77 - t6 * t88 + t204 * t64 - t202 * t63 + (t135 * t2 + t137 * t3) * t136, t138 * t3 - t59 * t78 + t7 * t88 + t208 * t64 - t202 * t112 + (-t1 * t135 + (qJD(2) * t58 - t11) * t186) * t136, -t138 * t2 + t59 * t77 - t7 * t87 - t208 * t63 + t204 * t112 + (-t1 * t137 + (-qJD(2) * t53 + t6) * t186) * t136, t1 * t59 + t202 * t11 + t2 * t53 + t204 * t6 - t208 * t7 + t3 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t149, t160, t30 * t44 - t31 * t45 + t39, t160, -t146 + t205, -t149, -t22 * t44 + t23 * t45 + t5, t160, -t149, t144, -t11 * t45 - t44 * t6 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t167, t215, -t112 * t23 + t16 * t64 + t12, t148, t215, t167, t64 * t7 + (qJD(5) + t11) * t112 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 * t64 - t146, t126 + t210, -t63 ^ 2 - t109, -t112 * t6 - t63 * t7 + t3;];
tauc_reg = t18;
