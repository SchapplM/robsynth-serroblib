% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:43
% EndTime: 2019-12-31 20:53:48
% DurationCPUTime: 1.46s
% Computational Cost: add. (1395->225), mult. (2512->251), div. (0->0), fcn. (1821->4), ass. (0->171)
t123 = cos(qJ(3));
t121 = sin(qJ(3));
t115 = t121 * pkin(3);
t185 = t123 * qJ(4);
t88 = t115 - t185;
t73 = t121 * qJ(5) + t88;
t191 = t73 * t123;
t124 = cos(qJ(2));
t206 = t124 * pkin(1);
t120 = pkin(3) + qJ(5);
t186 = t121 * qJ(4);
t210 = -t120 * t123 - t186;
t70 = -pkin(2) + t210;
t58 = t70 - t206;
t49 = t58 * t121;
t61 = t70 * t121;
t202 = t49 / 0.2e1 + t61 / 0.2e1;
t213 = t202 - t191;
t118 = t121 ^ 2;
t119 = t123 ^ 2;
t102 = t119 - t118;
t164 = qJD(1) + qJD(2);
t212 = t164 * t102;
t113 = t123 * qJD(4);
t142 = -t123 * pkin(3) - t186;
t211 = t142 * qJD(3) + t113;
t209 = pkin(4) + pkin(7);
t107 = -pkin(2) - t206;
t208 = t107 / 0.2e1;
t207 = -t121 / 0.2e1;
t3 = t58 * t73;
t85 = -pkin(2) + t142;
t71 = t85 - t206;
t205 = t71 * t88;
t4 = t73 * t70;
t204 = t88 * t85;
t122 = sin(qJ(2));
t106 = t122 * pkin(1) + pkin(7);
t203 = pkin(4) + t106;
t144 = (t118 + t119) * t124;
t72 = pkin(1) * t144;
t179 = t72 * qJD(1);
t63 = t72 * qJD(2);
t200 = t179 + t63;
t199 = pkin(1) * qJD(1);
t198 = pkin(1) * qJD(2);
t197 = t3 * qJD(1);
t196 = t58 * t123;
t195 = t70 * t123;
t194 = t71 * t121;
t193 = t71 * t123;
t192 = t73 * t121;
t190 = t85 * t121;
t189 = t85 * t123;
t188 = t88 * t121;
t74 = t203 * t121;
t75 = t203 * t123;
t9 = (t122 * t58 + (t121 * t74 + t123 * t75) * t124) * pkin(1);
t187 = t9 * qJD(1);
t22 = (t106 * t144 + t122 * t71) * pkin(1);
t184 = t22 * qJD(1);
t25 = t192 + t196;
t183 = t25 * qJD(1);
t26 = -t49 + t191;
t182 = t26 * qJD(1);
t31 = t188 + t193;
t181 = t31 * qJD(1);
t79 = t88 * t123;
t32 = t79 - t194;
t180 = t32 * qJD(1);
t178 = t75 * qJD(3);
t92 = t209 * t123;
t177 = t92 * qJD(3);
t167 = qJD(5) * t123;
t103 = t121 * t167;
t111 = t118 * qJD(4);
t176 = t103 + t111;
t163 = t122 * t198;
t101 = t123 * t163;
t104 = t121 * t113;
t175 = t104 - t101;
t174 = t119 * qJD(5) + t104;
t100 = t121 * t163;
t173 = t111 - t100;
t172 = qJD(1) * t121;
t171 = qJD(1) * t123;
t170 = qJD(2) * t121;
t169 = qJD(2) * t123;
t168 = qJD(4) * t121;
t166 = t120 * qJD(3);
t165 = t121 * qJD(3);
t114 = t123 * qJD(3);
t158 = t206 / 0.2e1;
t96 = t123 * t158;
t21 = -t195 / 0.2e1 - t196 / 0.2e1 + t96;
t94 = t121 * t158;
t24 = -t190 / 0.2e1 - t194 / 0.2e1 + t94;
t162 = pkin(7) * t165;
t161 = t122 * t199;
t160 = qJD(1) * t205;
t159 = -t206 / 0.2e1;
t157 = pkin(2) / 0.2e1 - t107 / 0.2e1;
t156 = t70 / 0.2e1 + t58 / 0.2e1;
t155 = t85 / 0.2e1 + t71 / 0.2e1;
t154 = t58 * t172;
t153 = t58 * t171;
t152 = t71 * t172;
t151 = t103 + t173;
t150 = -t101 + t174;
t149 = t107 * t172;
t148 = t107 * t171;
t147 = t106 * t165;
t146 = t185 / 0.2e1;
t145 = pkin(1) * t164;
t143 = t210 * qJD(3) - qJD(5) * t121 + t113;
t126 = (t120 * t207 + t146) * t206;
t1 = -t156 * t73 + t126;
t141 = t1 * qJD(1) - t4 * qJD(2);
t27 = t192 + t195;
t19 = t156 * t123 + t96;
t6 = t19 + t192;
t140 = t6 * qJD(1) + t27 * qJD(2);
t28 = -t61 + t191;
t95 = t121 * t159;
t5 = t95 - t213;
t139 = t5 * qJD(1) + t28 * qJD(2);
t23 = t155 * t121 + t94;
t14 = -t79 + t23;
t35 = t79 - t190;
t138 = t14 * qJD(1) - t35 * qJD(2);
t15 = t155 * t123 + t188 + t96;
t34 = t188 + t189;
t137 = t15 * qJD(1) + t34 * qJD(2);
t136 = qJD(3) * t88 - t168;
t135 = -t167 - t168;
t45 = t157 * t121 + t95;
t134 = pkin(2) * t170 + t45 * qJD(1);
t97 = t123 * t159;
t46 = t157 * t123 + t97;
t133 = pkin(2) * t169 + t46 * qJD(1);
t127 = (t146 - t115 / 0.2e1) * t206;
t11 = -t155 * t88 + t127;
t132 = t11 * qJD(1) - qJD(2) * t204;
t18 = t94 + t202;
t131 = t18 * qJD(1) + t70 * t170;
t130 = t23 * qJD(1) + t85 * t170;
t129 = t19 * qJD(1) + t70 * t169;
t117 = qJ(4) * qJD(4);
t116 = qJD(3) * qJ(4);
t110 = pkin(7) * t114;
t105 = t121 * t114;
t99 = t123 * t161;
t98 = t121 * t161;
t91 = t209 * t121;
t90 = t102 * qJD(3);
t87 = t164 * t119;
t86 = t164 * t118;
t84 = t106 * t114;
t80 = t91 * qJD(3);
t67 = t164 * t123 * t121;
t66 = t74 * qJD(3);
t48 = t97 + (-pkin(2) / 0.2e1 + t208) * t123;
t47 = pkin(2) * t207 + t121 * t208 + t95;
t20 = t94 - t202;
t17 = t79 + t24;
t16 = -t188 - t193 / 0.2e1 - t189 / 0.2e1 + t96;
t12 = t204 / 0.2e1 + t205 / 0.2e1 + t127;
t8 = t95 + t213;
t7 = t21 - t192;
t2 = t4 / 0.2e1 + t3 / 0.2e1 + t126;
t10 = [0, 0, 0, 0, -t163, -t124 * t198, t105, t90, 0, 0, 0, t107 * t165 - t101, t107 * t114 + t100, t63, t32 * qJD(3) - t175, -t31 * qJD(3) + t173, t22 * qJD(2) + t136 * t71, t63, -t25 * qJD(3) + t151, -t26 * qJD(3) + t150, t9 * qJD(2) + t3 * qJD(3) + t135 * t58; 0, 0, 0, 0, -t122 * t145, -t124 * t145, t105, t90, 0, 0, 0, t47 * qJD(3) - t101 - t99, t48 * qJD(3) + t100 + t98, t200, t17 * qJD(3) - t175 + t99, t16 * qJD(3) + t173 - t98, t184 + t12 * qJD(3) + t24 * qJD(4) + (pkin(7) * t144 + t122 * t85) * t198, t200, t7 * qJD(3) + t151 - t98, t8 * qJD(3) + t150 - t99, t187 + t2 * qJD(3) + t20 * qJD(4) + t21 * qJD(5) + (t122 * t70 + (t121 * t91 + t123 * t92) * t124) * t198; 0, 0, 0, 0, 0, 0, t67, t212, t114, -t165, 0, t47 * qJD(2) + t149 - t84, t48 * qJD(2) + t147 + t148, t211, t17 * qJD(2) + t180 + t84, t16 * qJD(2) - t147 - t181, t12 * qJD(2) + t106 * t211 + t160, t143, t7 * qJD(2) - t183 - t66, t8 * qJD(2) - t178 - t182, t197 + t2 * qJD(2) + (-t74 * qJ(4) - t75 * t120) * qJD(3) + t75 * qJD(4) - t74 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t67, t86, t24 * qJD(2) - t152 + t84, t114, t86, t67, t20 * qJD(2) - t154 + t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t67, t87, t21 * qJD(2) - t153 - t66; 0, 0, 0, 0, t161, t124 * t199, t105, t90, 0, 0, 0, -t45 * qJD(3) + t99, -t46 * qJD(3) - t98, -t179, -t14 * qJD(3) - t104 - t99, -t15 * qJD(3) + t111 + t98, -t11 * qJD(3) - t23 * qJD(4) - t184, -t179, -t6 * qJD(3) + t176 + t98, -t5 * qJD(3) + t174 + t99, -t1 * qJD(3) - t18 * qJD(4) - t19 * qJD(5) - t187; 0, 0, 0, 0, 0, 0, t105, t90, 0, 0, 0, -pkin(2) * t165, -pkin(2) * t114, 0, t35 * qJD(3) - t104, -t34 * qJD(3) + t111, t136 * t85, 0, -t27 * qJD(3) + t176, -t28 * qJD(3) + t174, t4 * qJD(3) + t135 * t70; 0, 0, 0, 0, 0, 0, t67, t212, t114, -t165, 0, -t110 - t134, -t133 + t162, t211, t110 - t138, -t137 - t162, pkin(7) * t211 - t132, t143, -t140 - t80, -t139 - t177, (-t91 * qJ(4) - t92 * t120) * qJD(3) + t92 * qJD(4) - t91 * qJD(5) - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t67, t86, t110 - t130, t114, t86, t67, -t131 + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t67, t87, -t129 - t80; 0, 0, 0, 0, 0, 0, -t67, -t212, 0, 0, 0, t45 * qJD(2) - t149, t46 * qJD(2) - t148, 0, t14 * qJD(2) - t180, t15 * qJD(2) + t181, t11 * qJD(2) - t160, 0, t6 * qJD(2) + t183, t5 * qJD(2) + t182, t1 * qJD(2) - t197; 0, 0, 0, 0, 0, 0, -t67, -t212, 0, 0, 0, t134, t133, 0, t138, t137, t132, 0, t140, t139, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t117, 0, qJD(4), qJD(5), t120 * qJD(5) + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t116, 0, qJD(3), 0, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t86, t23 * qJD(2) + t152, 0, -t86, -t67, t18 * qJD(2) + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t86, t130, 0, -t86, -t67, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t116, 0, -qJD(3), 0, -qJD(5) - t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t87, t19 * qJD(2) + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t87, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), qJD(4) - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
