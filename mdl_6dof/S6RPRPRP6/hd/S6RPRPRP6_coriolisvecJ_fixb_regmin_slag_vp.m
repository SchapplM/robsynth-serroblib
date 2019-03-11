% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:45
% EndTime: 2019-03-09 03:19:51
% DurationCPUTime: 2.00s
% Computational Cost: add. (3044->284), mult. (7725->363), div. (0->0), fcn. (5616->6), ass. (0->151)
t121 = sin(pkin(9));
t122 = cos(pkin(9));
t124 = sin(qJ(3));
t194 = cos(qJ(3));
t99 = t194 * t121 + t124 * t122;
t205 = t99 * qJD(1);
t208 = qJD(5) + t205;
t123 = sin(qJ(5));
t125 = cos(qJ(5));
t163 = t194 * t122;
t149 = qJD(1) * t163;
t178 = t121 * t124;
t162 = qJD(1) * t178;
t88 = -t149 + t162;
t66 = qJD(3) * t123 - t125 * t88;
t155 = t208 * t66;
t93 = t99 * qJD(3);
t129 = qJD(1) * t93;
t171 = qJD(5) * t125;
t172 = qJD(5) * t123;
t46 = qJD(3) * t172 - t123 * t129 - t88 * t171;
t209 = t46 - t155;
t130 = t205 * qJD(3);
t111 = qJD(3) * t149;
t173 = qJD(3) * t124;
t161 = t121 * t173;
t75 = qJD(1) * t161 - t111;
t31 = pkin(3) * t130 + t75 * qJ(4) - qJD(4) * t205;
t17 = pkin(8) * t129 + t31;
t188 = pkin(7) + qJ(2);
t105 = t188 * t121;
t100 = qJD(1) * t105;
t106 = t188 * t122;
t101 = qJD(1) * t106;
t160 = qJD(2) * t194;
t148 = qJD(1) * t160;
t168 = qJD(1) * qJD(2);
t158 = t124 * t168;
t159 = qJD(3) * t194;
t36 = -t100 * t173 + t101 * t159 + t121 * t148 + t122 * t158;
t21 = -pkin(4) * t75 + t36;
t60 = t194 * t100 + t124 * t101;
t204 = qJD(4) + t60;
t142 = pkin(4) * t205 + t204;
t197 = pkin(3) + pkin(8);
t32 = -t197 * qJD(3) + t142;
t167 = t123 * t21 + t125 * t17 + t32 * t171;
t116 = -pkin(2) * t122 - pkin(1);
t102 = t116 * qJD(1) + qJD(2);
t131 = -qJ(4) * t205 + t102;
t26 = t197 * t88 + t131;
t136 = -t26 * t172 + t167;
t68 = qJD(3) * t125 + t123 * t88;
t72 = t125 * t129;
t47 = t68 * qJD(5) - t72;
t2 = -qJ(6) * t47 - qJD(6) * t66 + t136;
t11 = -t123 * t26 + t125 * t32;
t6 = -qJ(6) * t68 + t11;
t5 = pkin(5) * t208 + t6;
t202 = t208 * t5 - t2;
t12 = t123 * t32 + t125 * t26;
t156 = -t123 * t17 + t125 * t21;
t128 = -t12 * qJD(5) + t156;
t1 = -pkin(5) * t75 + qJ(6) * t46 - qJD(6) * t68 + t128;
t7 = -qJ(6) * t66 + t12;
t203 = t208 * t7 + t1;
t211 = -t202 * t123 + t203 * t125;
t206 = t125 * t208;
t210 = t68 * t206;
t154 = t123 * t208;
t138 = -t125 * t75 - t154 * t208;
t120 = qJD(3) * qJ(4);
t61 = -t124 * t100 + t194 * t101;
t45 = -pkin(4) * t88 + t61;
t37 = t120 + t45;
t201 = t197 * t75 + t208 * t37;
t50 = (qJD(2) * t121 + qJD(3) * t106) * t124 + t105 * t159 - t122 * t160;
t200 = t68 ^ 2;
t199 = t88 ^ 2;
t198 = t205 ^ 2;
t196 = t5 - t6;
t169 = qJD(6) * t125;
t176 = qJ(6) + t197;
t39 = t125 * t45;
t184 = qJ(4) * t88;
t40 = t197 * t205 + t184;
t195 = t176 * t172 - t169 + pkin(5) * t88 - t39 - (-qJ(6) * t205 - t40) * t123;
t49 = pkin(3) * t88 + t131;
t193 = t49 * t205;
t192 = t66 * t88;
t191 = t68 * t88;
t98 = -t163 + t178;
t190 = t75 * t98;
t189 = t205 * t88;
t187 = t123 * t45 + t125 * t40;
t139 = -qJ(4) * t99 + t116;
t43 = t197 * t98 + t139;
t63 = t194 * t105 + t124 * t106;
t54 = t99 * pkin(4) + t63;
t186 = t123 * t54 + t125 * t43;
t104 = t176 * t125;
t180 = qJ(6) * t125;
t185 = -qJD(5) * t104 - qJD(6) * t123 - t180 * t205 - t187;
t64 = t124 * t105 - t194 * t106;
t183 = t125 * t46;
t179 = qJD(3) * t50;
t51 = t99 * qJD(2) - t64 * qJD(3);
t177 = t51 * qJD(3);
t174 = t121 ^ 2 + t122 ^ 2;
t170 = qJD(5) * t197;
t92 = -t122 * t159 + t161;
t143 = qJ(4) * t92 - qJD(4) * t99;
t23 = t197 * t93 + t143;
t28 = -t92 * pkin(4) + t51;
t166 = t123 * t28 + t125 * t23 + t54 * t171;
t165 = t98 * t172;
t164 = t98 * t171;
t157 = -qJ(6) * t98 - t43;
t151 = t174 * qJD(1) ^ 2;
t150 = -t100 * t159 - t101 * t173 - t121 * t158 + t122 * t148;
t119 = qJD(3) * qJD(4);
t33 = -t119 - t150;
t18 = -pkin(4) * t129 - t33;
t145 = -t18 * t98 - t37 * t93;
t144 = t208 * t93 - t190;
t141 = 0.2e1 * t174 * t168;
t137 = -t125 * t93 + t165;
t135 = -qJD(3) * t60 - t150;
t134 = qJD(3) * t61 - t36;
t133 = t123 * t75 - t206 * t208;
t27 = -pkin(4) * t93 - t50;
t13 = t47 * pkin(5) + t18;
t103 = t176 * t123;
t78 = qJD(3) * t88;
t65 = t66 ^ 2;
t62 = t75 * t99;
t59 = pkin(3) * t98 + t139;
t58 = pkin(3) * t205 + t184;
t57 = -t120 - t61;
t56 = -qJD(3) * pkin(3) + t204;
t55 = -pkin(4) * t98 - t64;
t53 = t125 * t54;
t42 = pkin(3) * t93 + t143;
t41 = t125 * t47;
t25 = t125 * t28;
t16 = pkin(5) * t66 + qJD(6) + t37;
t14 = t98 * t180 + t186;
t9 = pkin(5) * t99 + t157 * t123 + t53;
t4 = -t137 * qJ(6) + t98 * t169 - t43 * t172 + t166;
t3 = -pkin(5) * t92 + t25 + t157 * t171 + (-qJ(6) * t93 - qJD(5) * t54 - qJD(6) * t98 - t23) * t123;
t8 = [0, 0, 0, 0, 0, t141, qJ(2) * t141, -t205 * t92 - t62, -t99 * t130 - t205 * t93 + t92 * t88 + t190, -t92 * qJD(3), -t93 * qJD(3), 0, t102 * t93 + t116 * t129 - t177, -t102 * t92 - t116 * t75 + t179, t129 * t64 + t205 * t51 + t33 * t98 + t36 * t99 + t50 * t88 - t56 * t92 + t57 * t93 - t63 * t75, -t130 * t59 - t31 * t98 - t42 * t88 - t49 * t93 + t177, -t205 * t42 - t31 * t99 + t49 * t92 + t59 * t75 - t179, t31 * t59 + t33 * t64 + t36 * t63 + t42 * t49 + t50 * t57 + t51 * t56, t68 * t164 + (-t46 * t98 + t68 * t93) * t123 (-t123 * t66 + t125 * t68) * t93 + (-t123 * t47 - t183 + (-t123 * t68 - t125 * t66) * qJD(5)) * t98, t123 * t144 + t164 * t208 - t46 * t99 - t68 * t92, t125 * t144 - t165 * t208 - t47 * t99 + t66 * t92, -t208 * t92 - t62 (-t123 * t23 + t25) * t208 - (-t123 * t43 + t53) * t75 + t156 * t99 - t11 * t92 + t27 * t66 + t55 * t47 + t145 * t125 + (t37 * t123 * t98 - t12 * t99 - t186 * t208) * qJD(5), -t166 * t208 + t186 * t75 - t167 * t99 + t12 * t92 + t27 * t68 - t55 * t46 + t37 * t164 + ((t208 * t43 + t26 * t99) * qJD(5) - t145) * t123, -t14 * t47 - t3 * t68 - t4 * t66 + t46 * t9 + (-t123 * t5 + t125 * t7) * t93 + (-t1 * t123 + t125 * t2 + (-t123 * t7 - t125 * t5) * qJD(5)) * t98, t2 * t14 + t7 * t4 + t1 * t9 + t5 * t3 + t13 * ((-pkin(5) * t125 - pkin(4)) * t98 - t64) + t16 * (pkin(5) * t137 + t27); 0, 0, 0, 0, 0, -t151, -qJ(2) * t151, 0, 0, 0, 0, 0, 0.2e1 * t130, t111 + (-t88 - t162) * qJD(3), -t198 - t199, -0.2e1 * t130, t75 + t78, -t205 * t56 - t57 * t88 + t31, 0, 0, 0, 0, 0, t133 + t192, t191 - t138, -t123 * t209 + t210 - t41, -t123 * t203 - t125 * t202 + t16 * t88; 0, 0, 0, 0, 0, 0, 0, t189, t198 - t199, t111 + (t88 - t162) * qJD(3), 0, 0, -t102 * t205 + t134, t102 * t88 + t135, pkin(3) * t75 - qJ(4) * t129 - (t57 + t61) * t205 + (t56 - t204) * t88, t58 * t88 - t134 + t193, t205 * t58 - t49 * t88 + 0.2e1 * t119 - t135, -pkin(3) * t36 - qJ(4) * t33 - t204 * t57 - t49 * t58 - t56 * t61, -t154 * t68 - t183, -t41 - t210 + (t46 + t155) * t123, t138 + t191, t133 - t192, t208 * t88, qJ(4) * t47 + t11 * t88 - t39 * t208 + t142 * t66 + (t18 + (t40 + t170) * t208) * t123 + t201 * t125, -qJ(4) * t46 + t187 * t208 - t12 * t88 + t142 * t68 + (t170 * t208 + t18) * t125 - t201 * t123, t103 * t47 - t104 * t46 - t185 * t66 - t195 * t68 - t211, -t2 * t103 - t1 * t104 + t13 * (pkin(5) * t123 + qJ(4)) + t185 * t7 + t195 * t5 + (pkin(5) * t206 + t142) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 + t78, -t189, -qJD(3) ^ 2 - t198, qJD(3) * t57 + t193 + t36, 0, 0, 0, 0, 0, -qJD(3) * t66 + t138, -qJD(3) * t68 + t133, t209 * t125 + (t208 * t68 - t47) * t123, -qJD(3) * t16 + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t66, -t65 + t200, -t209, t72 + (-qJD(5) + t208) * t68, -t75, t12 * t208 - t37 * t68 + t128, t11 * t208 + t37 * t66 - t136, pkin(5) * t46 - t196 * t66, t196 * t7 + (-t16 * t68 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 - t200, t5 * t68 + t7 * t66 + t13;];
tauc_reg  = t8;
