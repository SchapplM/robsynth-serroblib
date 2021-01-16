% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:35
% EndTime: 2021-01-15 21:23:44
% DurationCPUTime: 2.06s
% Computational Cost: add. (2927->230), mult. (7960->328), div. (0->0), fcn. (6191->8), ass. (0->157)
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t145 = sin(pkin(9));
t146 = cos(pkin(9));
t152 = cos(qJ(2));
t193 = t146 * t152;
t180 = qJD(1) * t193;
t149 = sin(qJ(2));
t186 = qJD(1) * t149;
t112 = t145 * t186 - t180;
t148 = sin(qJ(4));
t151 = cos(qJ(4));
t124 = t145 * t152 + t146 * t149;
t188 = qJD(1) * t124;
t163 = t112 * t148 - t151 * t188;
t71 = -t151 * t112 - t148 * t188;
t26 = t147 * t71 - t150 * t163;
t62 = t150 * t71;
t27 = -t147 * t163 - t62;
t223 = t26 * t27;
t142 = qJD(2) + qJD(4);
t141 = qJD(5) + t142;
t197 = t141 * t27;
t114 = t124 * qJD(2);
t103 = qJD(1) * t114;
t182 = qJD(1) * qJD(2);
t178 = t149 * t182;
t133 = t145 * t178;
t177 = t152 * t182;
t104 = t146 * t177 - t133;
t155 = qJD(4) * t163 - t151 * t103 - t104 * t148;
t183 = qJD(5) * t147;
t184 = qJD(4) * t151;
t185 = qJD(4) * t148;
t23 = -t148 * t103 + t151 * t104 - t112 * t184 - t185 * t188;
t6 = qJD(5) * t62 + t147 * t155 + t150 * t23 + t163 * t183;
t219 = t6 + t197;
t156 = -qJD(5) * t26 - t147 * t23 + t150 * t155;
t198 = t141 * t26;
t214 = t156 + t198;
t220 = t26 ^ 2 - t27 ^ 2;
t205 = pkin(7) * t188;
t204 = -qJ(3) - pkin(6);
t132 = t204 * t152;
t129 = qJD(1) * t132;
t118 = t145 * t129;
t131 = t204 * t149;
t128 = qJD(1) * t131;
t200 = qJD(2) * pkin(2);
t122 = t128 + t200;
t73 = t146 * t122 + t118;
t49 = qJD(2) * pkin(3) - t205 + t73;
t206 = pkin(7) * t112;
t194 = t146 * t129;
t74 = t145 * t122 - t194;
t54 = t74 - t206;
t165 = -t148 * t49 - t151 * t54;
t226 = pkin(8) * t71;
t13 = -t165 + t226;
t138 = -pkin(2) * t152 - pkin(1);
t187 = qJD(1) * t138;
t130 = qJD(3) + t187;
t79 = pkin(3) * t112 + t130;
t40 = -pkin(4) * t71 + t79;
t230 = t13 * t183 + t40 * t27;
t172 = qJD(2) * t204;
t109 = qJD(3) * t152 + t149 * t172;
t91 = t109 * qJD(1);
t110 = -qJD(3) * t149 + t152 * t172;
t92 = t110 * qJD(1);
t55 = -t145 * t91 + t146 * t92;
t38 = -pkin(7) * t104 + t55;
t56 = t145 * t92 + t146 * t91;
t39 = -pkin(7) * t103 + t56;
t159 = -(qJD(4) * t49 + t39) * t151 - t148 * t38 + t54 * t185;
t2 = pkin(8) * t155 - t159;
t158 = qJD(4) * t165 - t148 * t39 + t151 * t38;
t3 = -pkin(8) * t23 + t158;
t215 = -t147 * t2 + t150 * t3 - t40 * t26;
t195 = t142 * t71;
t228 = t23 - t195;
t196 = t142 * t163;
t227 = t155 - t196;
t209 = pkin(4) * t163;
t225 = pkin(8) * t163;
t224 = t163 * t71;
t221 = t163 ^ 2 - t71 ^ 2;
t218 = -t79 * t71 + t159;
t217 = (-t13 * t141 - t3) * t147 + t230;
t216 = t163 * t79 + t158;
t213 = -0.2e1 * t182;
t137 = pkin(2) * t146 + pkin(3);
t208 = pkin(2) * t145;
t108 = t137 * t148 + t151 * t208;
t77 = -t128 * t145 + t194;
t57 = t77 + t206;
t78 = t146 * t128 + t118;
t58 = t78 - t205;
t212 = -t108 * qJD(4) + t148 * t58 - t151 * t57;
t162 = t137 * t151 - t148 * t208;
t211 = -t162 * qJD(4) + t148 * t57 + t151 * t58;
t210 = qJD(5) - t141;
t207 = pkin(2) * t149;
t203 = t212 + t226;
t202 = t211 + t225;
t65 = t146 * t109 + t145 * t110;
t199 = t13 * t150;
t154 = qJD(1) ^ 2;
t192 = t152 * t154;
t153 = qJD(2) ^ 2;
t191 = t153 * t149;
t190 = t153 * t152;
t82 = t145 * t131 - t146 * t132;
t189 = t149 ^ 2 - t152 ^ 2;
t140 = t149 * t200;
t139 = pkin(2) * t186;
t174 = -t148 * t54 + t151 * t49;
t12 = t174 + t225;
t10 = pkin(4) * t142 + t12;
t179 = -pkin(4) * t141 - t10;
t135 = pkin(2) * t178;
t80 = pkin(3) * t103 + t135;
t87 = pkin(3) * t114 + t140;
t86 = pkin(3) * t188 + t139;
t64 = -t109 * t145 + t146 * t110;
t169 = pkin(1) * t213;
t81 = t146 * t131 + t132 * t145;
t168 = 0.2e1 * t188;
t167 = -t10 * t147 - t199;
t123 = t145 * t149 - t193;
t75 = t151 * t123 + t124 * t148;
t76 = -t123 * t148 + t124 * t151;
t35 = t147 * t76 + t150 * t75;
t36 = -t147 * t75 + t150 * t76;
t66 = -pkin(7) * t124 + t81;
t67 = -pkin(7) * t123 + t82;
t164 = -t148 * t66 - t151 * t67;
t94 = pkin(3) * t123 + t138;
t117 = t123 * qJD(2);
t45 = pkin(7) * t117 + t64;
t46 = -pkin(7) * t114 + t65;
t160 = t148 * t45 + t151 * t46 + t66 * t184 - t67 * t185;
t157 = qJD(4) * t164 - t148 * t46 + t151 * t45;
t107 = pkin(4) + t162;
t50 = pkin(4) * t75 + t94;
t41 = t86 - t209;
t34 = qJD(4) * t76 + t151 * t114 - t117 * t148;
t33 = -qJD(4) * t75 - t114 * t148 - t117 * t151;
t19 = pkin(4) * t34 + t87;
t18 = -pkin(4) * t155 + t80;
t17 = -pkin(8) * t75 - t164;
t16 = -pkin(8) * t76 - t148 * t67 + t151 * t66;
t9 = qJD(5) * t36 + t147 * t33 + t150 * t34;
t8 = -qJD(5) * t35 - t147 * t34 + t150 * t33;
t5 = -pkin(8) * t33 + t157;
t4 = -pkin(8) * t34 + t160;
t1 = [0, 0, 0, 0.2e1 * t149 * t177, t189 * t213, t190, -t191, 0, -pkin(6) * t190 + t149 * t169, pkin(6) * t191 + t152 * t169, t103 * t138 + t114 * t130 + (t64 + (qJD(1) * t123 + t112) * t207) * qJD(2), t104 * t138 - t117 * t130 + (t168 * t207 - t65) * qJD(2), -t103 * t82 - t104 * t81 - t112 * t65 - t114 * t74 + t117 * t73 - t123 * t56 - t124 * t55 - t188 * t64, t55 * t81 + t56 * t82 + t64 * t73 + t65 * t74 + (t130 + t187) * t140, -t163 * t33 + t23 * t76, t155 * t76 + t163 * t34 - t23 * t75 + t33 * t71, t33 * t142, -t34 * t142, 0, t142 * t157 - t155 * t94 + t79 * t34 - t71 * t87 + t80 * t75, -t142 * t160 - t163 * t87 + t94 * t23 + t79 * t33 + t80 * t76, t26 * t8 + t36 * t6, t156 * t36 - t26 * t9 - t27 * t8 - t35 * t6, t8 * t141, -t9 * t141, 0, t19 * t27 - t50 * t156 + t18 * t35 + t40 * t9 + (-t147 * t4 + t150 * t5 + (-t147 * t16 - t150 * t17) * qJD(5)) * t141, t19 * t26 + t50 * t6 + t18 * t36 + t40 * t8 - (t147 * t5 + t150 * t4 + (-t147 * t17 + t150 * t16) * qJD(5)) * t141; 0, 0, 0, -t149 * t192, t189 * t154, 0, 0, 0, t154 * pkin(1) * t149, pkin(1) * t192, -qJD(2) * t77 - t112 * t139 - t130 * t188 + t55, qJD(2) * t78 + t112 * t130 - t139 * t188 - t56, (t74 + t77) * t188 + (-t73 + t78) * t112 + (-t103 * t145 - t104 * t146) * pkin(2), -t73 * t77 - t74 * t78 + (-t130 * t186 + t145 * t56 + t146 * t55) * pkin(2), t224, t221, t228, t227, 0, t212 * t142 + t71 * t86 + t216, t211 * t142 + t163 * t86 + t218, t223, t220, t219, t214, 0, -t41 * t27 + (t202 * t147 + t203 * t150) * t141 + ((-t107 * t147 - t108 * t150) * t141 + t167) * qJD(5) + t215, -t41 * t26 + (-t3 + (qJD(5) * t108 - t203) * t141) * t147 + (-qJD(5) * t10 - t2 + (-qJD(5) * t107 + t202) * t141) * t150 + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168 * qJD(2), -t133 + (-t112 + t180) * qJD(2), -t112 ^ 2 - t188 ^ 2, t112 * t74 + t188 * t73 + t135, 0, 0, 0, 0, 0, -t155 - t196, t23 + t195, 0, 0, 0, 0, 0, -t156 + t198, t6 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t221, t228, t227, 0, -t142 * t165 + t216, t142 * t174 + t218, t223, t220, t219, t214, 0, t27 * t209 - (-t12 * t147 - t199) * t141 + (t147 * t179 - t199) * qJD(5) + t215, t26 * t209 + (qJD(5) * t179 + t12 * t141 - t2) * t150 + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t220, t219, t214, 0, t210 * t167 + t215, (-t210 * t10 - t2) * t150 + t217;];
tauc_reg = t1;
