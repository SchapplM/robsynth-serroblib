% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:08
% EndTime: 2021-01-15 12:17:14
% DurationCPUTime: 1.65s
% Computational Cost: add. (1799->180), mult. (3438->272), div. (0->0), fcn. (3712->6), ass. (0->153)
t134 = sin(qJ(5));
t137 = cos(qJ(3));
t138 = -pkin(1) - pkin(6);
t234 = qJ(4) - t138;
t118 = t234 * t137;
t133 = cos(pkin(8));
t135 = sin(qJ(3));
t160 = t234 * t135;
t214 = sin(pkin(8));
t231 = -t214 * t118 - t133 * t160;
t216 = t231 * t134;
t136 = cos(qJ(5));
t215 = t231 * t136;
t233 = 0.2e1 * t134;
t123 = t214 * t135;
t211 = t133 * t137;
t114 = -t123 + t211;
t115 = t133 * t135 + t214 * t137;
t223 = qJD(3) * pkin(3);
t232 = (t214 * t114 - t115 * t133) * t223;
t112 = t114 ^ 2;
t113 = t115 ^ 2;
t78 = t112 + t113;
t158 = t112 / 0.2e1 + t113 / 0.2e1;
t77 = -t133 * t118 + t214 * t160;
t131 = t134 ^ 2;
t132 = t136 ^ 2;
t122 = t132 - t131;
t57 = t136 * t114;
t161 = t57 * t233;
t142 = qJD(1) * t161 - t122 * qJD(3);
t230 = -t77 / 0.2e1;
t229 = -t114 / 0.2e1;
t228 = -t134 / 0.2e1;
t227 = -t136 / 0.2e1;
t226 = t136 / 0.2e1;
t225 = t137 * pkin(3);
t125 = t135 * pkin(3) + qJ(2);
t139 = t115 * pkin(4) - t114 * pkin(7) + t125;
t18 = -t136 * t139 + t216;
t62 = t114 * pkin(4) + t115 * pkin(7) + t225;
t220 = t136 * t62;
t1 = t220 * t115 + (-t18 + t216) * t114;
t222 = t1 * qJD(1);
t221 = t134 * t62;
t19 = t134 * t139 + t215;
t2 = -t221 * t115 + (-t19 + t215) * t114;
t219 = t2 * qJD(1);
t175 = t77 / 0.2e1 + t230;
t4 = t175 * t115;
t218 = t4 * qJD(1);
t51 = t134 * t114;
t10 = t18 * t115 + t51 * t77;
t213 = t10 * qJD(1);
t11 = -t19 * t115 - t57 * t77;
t212 = t11 * qJD(1);
t17 = -t114 * t77 - t115 * t231;
t210 = t17 * qJD(1);
t25 = -0.1e1 / 0.2e1 - t158;
t20 = t25 * t134;
t209 = t20 * qJD(1);
t208 = t25 * qJD(1);
t27 = t25 * t136;
t207 = t27 * qJD(1);
t176 = t112 - t113;
t29 = t176 * t134;
t206 = t29 * qJD(1);
t30 = t78 * t134;
t205 = t30 * qJD(1);
t31 = t176 * t136;
t204 = t31 * qJD(1);
t141 = -t214 * t115 / 0.2e1 + t133 * t229;
t38 = (-t137 / 0.2e1 + t141) * pkin(3);
t203 = t38 * qJD(1);
t42 = t114 * t225 - t125 * t115;
t202 = t42 * qJD(1);
t43 = t125 * t114 + t115 * t225;
t201 = t43 * qJD(1);
t200 = t51 * qJD(1);
t52 = t134 * t115;
t199 = t52 * qJD(1);
t56 = t136 * t115;
t198 = t56 * qJD(1);
t64 = t78 * t136;
t197 = t64 * qJD(1);
t196 = t78 * qJD(1);
t195 = qJD(1) * qJ(2);
t194 = qJD(3) * t136;
t193 = qJD(5) * t134;
t192 = qJD(5) * t136;
t109 = t123 / 0.2e1 - t211 / 0.2e1;
t191 = t109 * qJD(1);
t190 = t114 * qJD(1);
t189 = t114 * qJD(3);
t188 = t114 * qJD(4);
t187 = t115 * qJD(1);
t186 = t115 * qJD(2);
t106 = t115 * qJD(3);
t185 = t115 * qJD(4);
t121 = t135 ^ 2 - t137 ^ 2;
t184 = t121 * qJD(1);
t183 = t125 * qJD(1);
t182 = t135 * qJD(1);
t181 = t135 * qJD(3);
t180 = t137 * qJD(1);
t179 = t137 * qJD(3);
t174 = qJ(2) * t182;
t173 = qJ(2) * t180;
t172 = t115 * t193;
t171 = t115 * t192;
t170 = t114 * t187;
t169 = t114 * t106;
t168 = t134 * t192;
t167 = t134 * t187;
t166 = t134 * t194;
t165 = t136 * t190;
t164 = t136 * t187;
t163 = t135 * t180;
t159 = -qJD(5) - t187;
t157 = -t51 * qJD(3) - t171;
t155 = qJD(3) * t161;
t12 = t125 * t225;
t154 = t12 * qJD(1) + t4 * qJD(2);
t124 = t214 * pkin(3) + pkin(7);
t126 = -t133 * pkin(3) - pkin(4);
t153 = -t114 * t124 - t115 * t126;
t152 = t159 * t136;
t149 = t124 * t115 / 0.2e1 + t126 * t229;
t140 = t62 / 0.2e1 + t149;
t8 = t175 * t134 - t140 * t136;
t151 = -t126 * t134 * qJD(3) - t8 * qJD(1);
t6 = t140 * t134 + t175 * t136;
t150 = -t6 * qJD(1) - t126 * t194;
t49 = (t131 / 0.2e1 - t132 / 0.2e1) * t114;
t148 = -t49 * qJD(1) + t166;
t147 = t114 * t152;
t146 = -t109 * qJD(5) + t170;
t144 = t136 * t112 * t134 * qJD(1) + t49 * qJD(3);
t63 = t122 * t112;
t143 = t63 * qJD(1) + t155;
t103 = t109 * qJD(3);
t102 = t136 * t189;
t45 = t52 * qJD(5);
t44 = t49 * qJD(5);
t37 = t225 / 0.2e1 + t141 * pkin(3);
t36 = -t193 - t199;
t28 = t78 * t227 + t226;
t24 = 0.1e1 / 0.2e1 - t158;
t21 = t158 * t134 + t228;
t9 = t134 * t230 + t77 * t228 + t220 / 0.2e1 - t149 * t136;
t7 = -t221 / 0.2e1 + t149 * t134 + (-t226 + t227) * t77;
t3 = t4 * qJD(3);
t5 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t135 * t179, t121 * qJD(3), 0, 0, 0, qJ(2) * t179 + qJD(2) * t135, -qJ(2) * t181 + qJD(2) * t137, t43 * qJD(3) + t186, t114 * qJD(2) + t42 * qJD(3), t78 * qJD(4), t125 * qJD(2) + t12 * qJD(3) + t17 * qJD(4), -t112 * t168 - t132 * t169, -t63 * qJD(5) + t115 * t155, t31 * qJD(3) - t114 * t172, -t29 * qJD(3) - t114 * t171, t169, t1 * qJD(3) + t30 * qJD(4) + t11 * qJD(5) + t136 * t186, t2 * qJD(3) + t64 * qJD(4) + t10 * qJD(5) - t134 * t186; 0, 0, 0, 0, qJD(1), t195, 0, 0, 0, 0, 0, t182, t180, t187, t190, 0, t24 * qJD(4) + t183 + t3, 0, 0, 0, 0, 0, t28 * qJD(5) + t164, t21 * qJD(5) - t167; 0, 0, 0, 0, 0, 0, -t163, t184, -t181, -t179, 0, -t138 * t181 + t173, -t138 * t179 - t174, -qJD(3) * t231 + t201, -t77 * qJD(3) + t202, -t232, (-t133 * t231 + t214 * t77) * t223 + t37 * qJD(4) + t154, -t44 + (-t132 * t190 - t166) * t115, -0.2e1 * t114 * t168 + t142 * t115, t134 * t189 + t204, t102 - t206, t146, t222 + (t134 * t153 - t215) * qJD(3) + t9 * qJD(5), t219 + (t136 * t153 + t216) * qJD(3) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t24 * qJD(2) + t37 * qJD(3) + t210, 0, 0, 0, 0, 0, t205, t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -t143, t159 * t51, t147, -t103, t28 * qJD(2) + t9 * qJD(3) - t19 * qJD(5) + t212, t21 * qJD(2) + t7 * qJD(3) + t18 * qJD(5) + t213; 0, 0, 0, 0, -qJD(1), -t195, 0, 0, 0, 0, 0, -t182, -t180, -t187, -t190, 0, t25 * qJD(4) - t183 + t3, 0, 0, 0, 0, 0, t27 * qJD(5) - t164, -t20 * qJD(5) + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t179, -t106, -t189, 0, t218 + t232, 0, 0, 0, 0, 0, -t51 * qJD(5) - t106 * t136, -t57 * qJD(5) + t106 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157 + t207, -t57 * qJD(3) + t172 - t209; 0, 0, 0, 0, 0, 0, t163, -t184, 0, 0, 0, -t173, t174, -t188 - t201, t185 - t202, 0, t38 * qJD(4) - t154, t132 * t170 - t44, t147 * t233, t56 * qJD(5) - t204, -t45 + t206, -t146, t8 * qJD(5) - t136 * t188 - t222, t51 * qJD(4) + t6 * qJD(5) - t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t122 * qJD(5), 0, 0, 0, t126 * t193, t126 * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t187, 0, t203, 0, 0, 0, 0, 0, -t165, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t142, t192 + t198, t36, t191, -t124 * t192 - t151, t124 * t193 - t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, -t106, -t196, -t25 * qJD(2) - t38 * qJD(3) - t210, 0, 0, 0, 0, 0, t102 - t45 - t205, t157 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, -t187, 0, -t203, 0, 0, 0, 0, 0, t165, -t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t143, -t56 * qJD(3) + t114 * t167, t52 * qJD(3) + t114 * t164, -t103, -t27 * qJD(2) - t8 * qJD(3) + t52 * qJD(4) - t212, t20 * qJD(2) - t6 * qJD(3) + t136 * t185 - t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t142, -t198, t199, -t191, t151, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
