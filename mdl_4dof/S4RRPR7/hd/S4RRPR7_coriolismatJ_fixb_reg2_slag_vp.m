% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:47
% EndTime: 2019-12-31 17:06:50
% DurationCPUTime: 1.77s
% Computational Cost: add. (2194->175), mult. (4484->266), div. (0->0), fcn. (4701->6), ass. (0->152)
t125 = sin(qJ(4));
t126 = sin(qJ(2));
t223 = qJ(3) + pkin(5);
t109 = t223 * t126;
t124 = sin(pkin(7));
t128 = cos(qJ(2));
t174 = t223 * t128;
t199 = cos(pkin(7));
t75 = t199 * t109 + t124 * t174;
t225 = t75 * t125;
t127 = cos(qJ(4));
t224 = t75 * t127;
t221 = -t124 * t109 + t199 * t174;
t205 = t221 * t125;
t204 = t221 * t127;
t148 = t199 * t126;
t194 = t124 * t128;
t105 = t148 + t194;
t193 = t125 * t127;
t147 = 0.2e1 * t105 * t193;
t103 = t124 * t126 - t199 * t128;
t101 = t103 ^ 2;
t102 = t105 ^ 2;
t222 = -t102 - t101;
t160 = t102 - t101;
t123 = t127 ^ 2;
t220 = -t123 / 0.2e1;
t219 = t126 * pkin(2);
t217 = qJD(2) * pkin(2);
t65 = t105 * pkin(3) + t103 * pkin(6) + t219;
t215 = t125 * t65;
t36 = -t224 + t215;
t209 = t36 * t125;
t213 = t127 * t65;
t35 = t213 + t225;
t210 = t35 * t127;
t118 = -t128 * pkin(2) - pkin(1);
t132 = t103 * pkin(3) - t105 * pkin(6) + t118;
t33 = -t127 * t132 + t205;
t34 = t125 * t132 + t204;
t1 = (t209 + t210) * t105 + (-t34 * t125 + t127 * t33) * t103;
t216 = t1 * qJD(1);
t2 = t75 * t221 - t33 * t35 + t34 * t36;
t211 = t2 * qJD(1);
t4 = (-t33 + t205) * t105 + (t35 - t225) * t103;
t208 = t4 * qJD(1);
t5 = (-t34 + t204) * t105 + (-t36 - t224) * t103;
t207 = t5 * qJD(1);
t203 = t75 * t105;
t6 = t203 + (-t125 * t33 - t34 * t127) * t103;
t206 = t6 * qJD(1);
t122 = t125 ^ 2;
t116 = t124 * pkin(2) + pkin(6);
t196 = t116 * t103;
t117 = -t199 * pkin(2) - pkin(3);
t197 = t105 * t117;
t129 = (t220 - t122 / 0.2e1) * t196 + t197 / 0.2e1;
t140 = -t210 / 0.2e1 - t209 / 0.2e1;
t8 = t129 + t140;
t200 = t8 * qJD(1);
t62 = t125 * t105;
t21 = t33 * t103 - t75 * t62;
t192 = t21 * qJD(1);
t22 = -t34 * t103 + t105 * t224;
t191 = t22 * qJD(1);
t23 = t118 * t219;
t190 = t23 * qJD(1);
t24 = -t103 * t221 + t203;
t189 = t24 * qJD(1);
t39 = t160 * t125;
t188 = t39 * qJD(1);
t40 = t222 * t125;
t187 = t40 * qJD(1);
t41 = t160 * t127;
t186 = t41 * qJD(1);
t131 = -t124 * t103 / 0.2e1 - t199 * t105 / 0.2e1;
t46 = (-t126 / 0.2e1 + t131) * pkin(2);
t185 = t46 * qJD(1);
t184 = t160 * qJD(1);
t51 = t103 * t219 + t118 * t105;
t183 = t51 * qJD(1);
t52 = -t118 * t103 + t105 * t219;
t182 = t52 * qJD(1);
t58 = (t122 / 0.2e1 + t220) * t105;
t181 = t58 * qJD(4);
t59 = t125 * t103;
t180 = t59 * qJD(1);
t179 = t62 * qJD(1);
t64 = t127 * t103;
t178 = t64 * qJD(1);
t97 = t122 * t103;
t98 = t123 * t103;
t66 = t97 + t98;
t177 = t66 * qJD(1);
t68 = t222 * t127;
t176 = t68 * qJD(1);
t175 = t222 * qJD(1);
t112 = t123 - t122;
t173 = qJD(1) * t128;
t172 = qJD(2) * t127;
t171 = qJD(4) * t125;
t170 = qJD(4) * t127;
t100 = t148 / 0.2e1 + t194 / 0.2e1;
t169 = t100 * qJD(1);
t168 = t103 * qJD(1);
t167 = t103 * qJD(3);
t166 = t105 * qJD(1);
t165 = t105 * qJD(2);
t164 = t105 * qJD(3);
t113 = -t126 ^ 2 + t128 ^ 2;
t163 = t113 * qJD(1);
t162 = t126 * qJD(2);
t161 = t128 * qJD(2);
t159 = pkin(1) * t126 * qJD(1);
t158 = pkin(1) * t173;
t155 = t103 * t170;
t154 = t103 * t166;
t153 = t103 * t165;
t152 = t125 * t170;
t151 = t125 * t172;
t150 = t126 * t161;
t149 = t127 * t166;
t146 = -qJD(4) - t168;
t145 = t102 * t152;
t144 = qJD(2) * t147;
t143 = -t35 * t125 + t36 * t127;
t142 = -t103 * t117 - t105 * t116;
t141 = t146 * t127;
t139 = t196 / 0.2e1 - t197 / 0.2e1;
t130 = t65 / 0.2e1 + t139;
t19 = t130 * t127;
t138 = -t117 * t125 * qJD(2) + t19 * qJD(1);
t17 = t130 * t125;
t137 = -t17 * qJD(1) - t117 * t172;
t43 = -t58 * qJD(1) + t151;
t136 = t105 * t141;
t135 = t100 * qJD(4) + t154;
t37 = t102 * qJD(1) * t193 + t58 * qJD(2);
t67 = t112 * t102;
t134 = t67 * qJD(1) + t144;
t133 = qJD(1) * t147 - t112 * qJD(2);
t114 = t126 * t173;
t99 = t103 * qJD(2);
t96 = t100 * qJD(2);
t95 = t127 * t165;
t55 = t59 * qJD(4);
t45 = t219 / 0.2e1 + t131 * pkin(2);
t44 = -t171 - t180;
t20 = t225 + t213 / 0.2e1 - t139 * t127;
t18 = t224 - t215 / 0.2e1 + t139 * t125;
t7 = t129 - t140;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t113 * qJD(2), 0, -t150, 0, 0, -pkin(1) * t162, -pkin(1) * t161, 0, 0, -t153, -t160 * qJD(2), 0, t153, 0, 0, t51 * qJD(2), t52 * qJD(2), -qJD(3) * t222, t23 * qJD(2) + t24 * qJD(3), -t123 * t153 - t145, -t67 * qJD(4) + t103 * t144, -t103 * t105 * t171 + t41 * qJD(2), -t122 * t153 + t145, -t39 * qJD(2) - t105 * t155, t153, t4 * qJD(2) - t40 * qJD(3) + t22 * qJD(4), t5 * qJD(2) - t68 * qJD(3) + t21 * qJD(4), -t1 * qJD(2), t2 * qJD(2) + t6 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t163, t161, -t114, -t162, 0, -pkin(5) * t161 - t159, pkin(5) * t162 - t158, 0, 0, -t154, -t184, -t99, t154, -t165, 0, -qJD(2) * t221 + t183, qJD(2) * t75 + t182, (t199 * t103 - t105 * t124) * t217, t190 + (-t124 * t75 - t199 * t221) * t217 + t45 * qJD(3), -t181 + (-t123 * t166 - t151) * t103, (t97 - t98) * qJD(2) + (-qJD(4) + t168) * t147, t125 * t165 + t186, t181 + (-t122 * t166 + t151) * t103, t95 - t188, t135, t208 + (t125 * t142 - t204) * qJD(2) + t20 * qJD(4), t207 + (t127 * t142 + t205) * qJD(2) + t18 * qJD(4), qJD(2) * t143 - t216, t211 + (t116 * t143 + t117 * t221) * qJD(2) + t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t45 * qJD(2) + t189, 0, 0, 0, 0, 0, 0, -t187, -t176, 0, t7 * qJD(2) + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t134, t146 * t62, t37, t136, t96, t20 * qJD(2) - t34 * qJD(4) + t191, t18 * qJD(2) + t33 * qJD(4) + t192, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t163, 0, t114, 0, 0, t159, t158, 0, 0, t154, t184, 0, -t154, 0, 0, -t164 - t183, t167 - t182, 0, t46 * qJD(3) - t190, t123 * t154 - t181, 0.2e1 * t125 * t136, t64 * qJD(4) - t186, t122 * t154 + t181, -t55 + t188, -t135, -t19 * qJD(4) - t127 * t164 - t208, t62 * qJD(3) + t17 * qJD(4) - t207, -t66 * qJD(3) + t216, t8 * qJD(3) - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t112 * qJD(4), 0, -t152, 0, 0, t117 * t171, t117 * t170, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t168, 0, t185, 0, 0, 0, 0, 0, 0, -t149, t179, -t177, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t133, t170 + t178, -t43, t44, -t169, -t116 * t170 - t138, t116 * t171 - t137, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t99, t175, -t46 * qJD(2) - t189, 0, 0, 0, 0, 0, 0, -t55 + t95 + t187, -t62 * qJD(2) - t155 + t176, t66 * qJD(2), -t8 * qJD(2) - t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, -t168, 0, -t185, 0, 0, 0, 0, 0, 0, t149, -t179, t177, -t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t134, -t64 * qJD(2) + t125 * t154, -t37, t59 * qJD(2) + t103 * t149, t96, t19 * qJD(2) + t59 * qJD(3) - t191, -t17 * qJD(2) + t127 * t167 - t192, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t133, -t178, t43, t180, t169, t138, t137, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t127 * t168, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
