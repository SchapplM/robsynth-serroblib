% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:27
% EndTime: 2022-01-20 11:49:31
% DurationCPUTime: 1.41s
% Computational Cost: add. (2445->215), mult. (4150->272), div. (0->0), fcn. (2668->6), ass. (0->161)
t122 = cos(qJ(4));
t119 = sin(qJ(4));
t123 = cos(qJ(3));
t116 = qJD(1) + qJD(2);
t121 = sin(qJ(2));
t191 = pkin(1) * qJD(1);
t166 = t121 * t191;
t215 = pkin(7) + pkin(8);
t155 = t215 * t116 + t166;
t61 = t155 * t123;
t53 = t119 * t61;
t120 = sin(qJ(3));
t60 = t155 * t120;
t56 = qJD(3) * pkin(3) - t60;
t151 = t122 * t56 - t53;
t86 = t119 * t123 + t122 * t120;
t70 = t86 * t116;
t185 = t70 * qJ(5);
t14 = t151 - t185;
t172 = qJD(4) * t119;
t145 = qJD(3) * t155;
t124 = cos(qJ(2));
t190 = pkin(1) * qJD(2);
t163 = qJD(1) * t190;
t147 = t124 * t163;
t38 = -t120 * t145 + t123 * t147;
t39 = -t120 * t147 - t123 * t145;
t133 = -(qJD(4) * t56 + t38) * t122 - t119 * t39 + t61 * t172;
t115 = qJD(3) + qJD(4);
t44 = t115 * t86;
t31 = t44 * t116;
t214 = -t31 * qJ(5) - t133;
t171 = qJD(4) * t122;
t174 = qJD(3) * t123;
t213 = -t122 * t174 - t123 * t171;
t165 = t124 * t191;
t179 = t122 * t123;
t181 = t119 * t120;
t85 = -t179 + t181;
t161 = qJD(3) * t215;
t87 = t120 * t161;
t88 = t123 * t161;
t98 = t215 * t120;
t113 = t123 * pkin(8);
t99 = t123 * pkin(7) + t113;
t212 = t119 * t88 + t122 * t87 - t85 * t165 + t98 * t171 + t99 * t172;
t141 = t119 * t98 - t122 * t99;
t211 = t141 * qJD(4) + t119 * t87 - t122 * t88 + t86 * t165;
t210 = t70 ^ 2;
t208 = t85 * pkin(4);
t207 = pkin(3) * t115;
t206 = t124 * pkin(1);
t162 = t116 * t181;
t68 = -t116 * t179 + t162;
t205 = t70 * t68;
t109 = -t123 * pkin(3) - pkin(2);
t72 = t109 * t116 - t165;
t204 = t72 * t70;
t106 = t121 * pkin(1) + pkin(7);
t203 = -pkin(8) - t106;
t196 = -t44 * qJ(5) - t85 * qJD(5);
t202 = t196 - t212;
t144 = t115 * t181;
t43 = t144 + t213;
t140 = t43 * qJ(5) - t86 * qJD(5);
t201 = t140 + t211;
t13 = t115 * pkin(4) + t14;
t200 = t13 - t14;
t104 = t121 * t163;
t175 = qJD(3) * t120;
t160 = t116 * t175;
t74 = pkin(3) * t160 + t104;
t24 = t31 * pkin(4) + t74;
t154 = t68 * pkin(4) + qJD(5);
t33 = t154 + t72;
t199 = t24 * t85 + t33 * t44;
t198 = t24 * t86 - t33 * t43;
t197 = t72 * t44 + t74 * t85;
t195 = -t72 * t43 + t74 * t86;
t194 = -t122 * t60 - t53;
t92 = -t116 * pkin(2) - t165;
t193 = t120 * t104 + t92 * t174;
t192 = t213 * t116;
t55 = t122 * t61;
t30 = t116 * t144 + t192;
t189 = t30 * qJ(5);
t187 = t68 * qJ(5);
t186 = t68 * t115;
t183 = t86 * qJ(5);
t182 = t116 * t120;
t180 = t121 * t123;
t125 = qJD(3) ^ 2;
t178 = t125 * t120;
t112 = t125 * t123;
t177 = qJD(5) + t33;
t176 = t120 ^ 2 - t123 ^ 2;
t173 = qJD(3) * t124;
t170 = -qJD(1) - t116;
t169 = -qJD(2) + t116;
t168 = pkin(3) * t182;
t167 = t124 * t190;
t111 = t121 * t190;
t110 = pkin(3) * t175;
t143 = -t119 * t56 - t55;
t15 = -t143 - t187;
t4 = -t68 * qJD(5) + t214;
t153 = -t119 * t38 + t122 * t39;
t132 = t143 * qJD(4) + t153;
t128 = t132 + t189;
t5 = -t70 * qJD(5) + t128;
t164 = t13 * t43 - t15 * t44 - t4 * t85 - t5 * t86;
t159 = t116 * t174;
t158 = t120 * t173;
t37 = t44 * pkin(4) + t110;
t150 = t119 * t60 - t55;
t149 = qJD(3) * t203;
t146 = t37 - t166;
t82 = t203 * t120;
t83 = t123 * t106 + t113;
t142 = -t119 * t82 - t122 * t83;
t95 = t109 - t206;
t57 = t120 * t149 + t123 * t167;
t58 = -t120 * t167 + t123 * t149;
t139 = t119 * t58 + t122 * t57 + t82 * t171 - t83 * t172;
t137 = -t166 + t110;
t136 = -t116 * t92 - t147;
t135 = -t121 * t182 + t123 * t173;
t134 = -t115 * t162 - t192;
t131 = t142 * qJD(4) - t119 * t57 + t122 * t58;
t129 = t72 * t68 + t133;
t127 = (-t55 + (-t56 - t207) * t119) * qJD(4) + t153;
t126 = t177 * t68 - t214;
t114 = t116 ^ 2;
t108 = -pkin(2) - t206;
t107 = t122 * pkin(3) + pkin(4);
t100 = t171 * t207;
t90 = 0.2e1 * t120 * t159;
t89 = t111 + t110;
t81 = t85 * qJ(5);
t75 = t92 * t175;
t67 = -0.2e1 * t176 * t116 * qJD(3);
t66 = t109 + t208;
t65 = t68 ^ 2;
t59 = t95 + t208;
t47 = t70 * pkin(4) + t168;
t41 = t44 * t115;
t40 = t43 * t115;
t35 = -t141 - t81;
t34 = -t119 * t99 - t122 * t98 - t183;
t32 = t111 + t37;
t27 = -t142 - t81;
t26 = -t119 * t83 + t122 * t82 - t183;
t25 = -t65 + t210;
t22 = t134 + t186;
t19 = -t185 + t194;
t18 = t150 + t187;
t10 = -t30 * t86 - t70 * t43;
t7 = t131 + t140;
t6 = t139 + t196;
t3 = t30 * t85 - t86 * t31 + t43 * t68 - t70 * t44;
t1 = [0, 0, 0, 0, -t116 * t111 - t104, t170 * t167, t90, t67, t112, -t178, 0, t108 * t160 - t106 * t112 + t75 + (t170 * t180 - t158) * t190, t106 * t178 + t108 * t159 - t135 * t190 + t193, t10, t3, -t40, -t41, 0, t131 * t115 + t95 * t31 + t89 * t68 + t197, -t139 * t115 - t95 * t30 + t89 * t70 + t195, t7 * t115 + t59 * t31 + t32 * t68 + t199, -t6 * t115 - t59 * t30 + t32 * t70 + t198, t26 * t30 - t27 * t31 - t6 * t68 - t7 * t70 + t164, t13 * t7 + t15 * t6 + t24 * t59 + t5 * t26 + t4 * t27 + t33 * t32; 0, 0, 0, 0, t116 * t166 - t104, t169 * t165, t90, t67, t112, -t178, 0, -pkin(2) * t160 - pkin(7) * t112 + t75 + (t169 * t180 + t158) * t191, -pkin(2) * t159 + pkin(7) * t178 + t135 * t191 + t193, t10, t3, -t40, -t41, 0, t109 * t31 + t211 * t115 + t137 * t68 + t197, -t109 * t30 + t212 * t115 + t137 * t70 + t195, t201 * t115 + t146 * t68 + t66 * t31 + t199, -t202 * t115 + t146 * t70 - t66 * t30 + t198, -t201 * t70 - t202 * t68 + t34 * t30 - t35 * t31 + t164, t201 * t13 + t146 * t33 + t202 * t15 + t24 * t66 + t5 * t34 + t4 * t35; 0, 0, 0, 0, 0, 0, -t120 * t114 * t123, t176 * t114, 0, 0, 0, t136 * t120, t136 * t123, t205, t25, t22, 0, 0, -t150 * t115 - t68 * t168 + t127 - t204, t194 * t115 - t70 * t168 - t100 + t129, -t18 * t115 - t177 * t70 - t47 * t68 + t127 + t189, t19 * t115 - t47 * t70 - t100 + t126, t107 * t30 + (t15 + t18) * t70 + (-t13 + t19) * t68 + (-t119 * t31 + (t119 * t70 - t122 * t68) * qJD(4)) * pkin(3), t5 * t107 - t13 * t18 - t15 * t19 - t33 * t47 + (t119 * t4 + (-t119 * t13 + t122 * t15) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, t25, t22, 0, 0, -t143 * t115 + t132 - t204, t151 * t115 + t129, t15 * t115 + (-t154 - t33) * t70 + t128, -t210 * pkin(4) + t14 * t115 + t126, t30 * pkin(4) - t200 * t68, t200 * t15 + (-t33 * t70 + t5) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t115 + t31, t134 - t186, -t65 - t210, t13 * t70 + t15 * t68 + t24;];
tauc_reg = t1;
