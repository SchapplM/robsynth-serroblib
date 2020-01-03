% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:46
% EndTime: 2019-12-31 20:53:50
% DurationCPUTime: 2.29s
% Computational Cost: add. (1484->262), mult. (2793->257), div. (0->0), fcn. (2042->4), ass. (0->177)
t129 = cos(qJ(3));
t127 = sin(qJ(3));
t121 = t127 * pkin(3);
t195 = t129 * qJ(4);
t90 = t121 - t195;
t75 = t127 * qJ(5) + t90;
t201 = t75 * t129;
t130 = cos(qJ(2));
t215 = t130 * pkin(1);
t126 = pkin(3) + qJ(5);
t196 = t127 * qJ(4);
t220 = -t126 * t129 - t196;
t72 = -pkin(2) + t220;
t60 = t72 - t215;
t51 = t60 * t127;
t63 = t72 * t127;
t211 = t51 / 0.2e1 + t63 / 0.2e1;
t224 = t211 - t201;
t124 = t127 ^ 2;
t125 = t129 ^ 2;
t223 = t124 + t125;
t107 = t125 - t124;
t172 = qJD(1) + qJD(2);
t222 = t172 * t107;
t119 = t129 * qJD(4);
t148 = -t129 * pkin(3) - t196;
t221 = t148 * qJD(3) + t119;
t219 = pkin(4) + pkin(7);
t112 = -pkin(2) - t215;
t218 = t112 / 0.2e1;
t217 = -t127 / 0.2e1;
t128 = sin(qJ(2));
t216 = t128 * pkin(1);
t3 = t60 * t75;
t87 = -pkin(2) + t148;
t73 = t87 - t215;
t214 = t73 * t90;
t4 = t75 * t72;
t213 = t90 * t87;
t111 = pkin(7) + t216;
t212 = pkin(4) + t111;
t151 = t223 * t130;
t74 = pkin(1) * t151;
t188 = t74 * qJD(1);
t65 = t74 * qJD(2);
t32 = t188 + t65;
t209 = pkin(1) * qJD(1);
t208 = pkin(1) * qJD(2);
t207 = t3 * qJD(1);
t206 = t60 * t129;
t205 = t72 * t129;
t204 = t73 * t127;
t203 = t73 * t129;
t202 = t75 * t127;
t200 = t87 * t127;
t199 = t87 * t129;
t76 = t212 * t127;
t77 = t212 * t129;
t9 = (t128 * t60 + (t127 * t76 + t129 * t77) * t130) * pkin(1);
t198 = t9 * qJD(1);
t197 = t90 * t127;
t150 = t111 * t151;
t22 = (t128 * t73 + t150) * pkin(1);
t194 = t22 * qJD(1);
t25 = t202 + t206;
t193 = t25 * qJD(1);
t26 = -t51 + t201;
t192 = t26 * qJD(1);
t29 = (t112 * t128 + t150) * pkin(1);
t191 = t29 * qJD(1);
t33 = t197 + t203;
t190 = t33 * qJD(1);
t81 = t90 * t129;
t34 = t81 - t204;
t189 = t34 * qJD(1);
t187 = t77 * qJD(3);
t95 = t219 * t129;
t186 = t95 * qJD(3);
t170 = t128 * t208;
t104 = t129 * t170;
t109 = t127 * t119;
t185 = t104 - t109;
t184 = t223 * pkin(7) * t215;
t175 = qJD(5) * t129;
t108 = t127 * t175;
t116 = t124 * qJD(4);
t183 = t108 + t116;
t182 = t125 * qJD(5) + t109;
t103 = t127 * t170;
t181 = t116 - t103;
t180 = qJD(1) * t127;
t179 = qJD(1) * t129;
t178 = qJD(2) * t127;
t177 = qJD(2) * t129;
t176 = qJD(4) * t127;
t174 = t107 * qJD(3);
t173 = t126 * qJD(3);
t117 = t127 * qJD(3);
t120 = t129 * qJD(3);
t165 = t215 / 0.2e1;
t99 = t129 * t165;
t21 = -t205 / 0.2e1 - t206 / 0.2e1 + t99;
t97 = t127 * t165;
t24 = -t200 / 0.2e1 - t204 / 0.2e1 + t97;
t169 = pkin(7) * t117;
t168 = t128 * t209;
t167 = qJD(1) * t214;
t166 = -t215 / 0.2e1;
t164 = pkin(2) / 0.2e1 - t112 / 0.2e1;
t163 = t72 / 0.2e1 + t60 / 0.2e1;
t162 = t87 / 0.2e1 + t73 / 0.2e1;
t161 = t60 * t180;
t160 = t60 * t179;
t159 = t73 * t180;
t158 = -t104 + t182;
t101 = t127 * t168;
t157 = -t101 + t181;
t156 = t112 * t180;
t155 = t112 * t179;
t154 = t111 * t117;
t153 = t195 / 0.2e1;
t152 = pkin(1) * t172;
t149 = t220 * qJD(3) - qJD(5) * t127 + t119;
t132 = (t126 * t217 + t153) * t215;
t1 = -t163 * t75 + t132;
t147 = t1 * qJD(1) - t4 * qJD(2);
t27 = t202 + t205;
t19 = t163 * t129 + t99;
t6 = t19 + t202;
t146 = t6 * qJD(1) + t27 * qJD(2);
t28 = -t63 + t201;
t98 = t127 * t166;
t5 = t98 - t224;
t145 = t5 * qJD(1) + t28 * qJD(2);
t23 = t162 * t127 + t97;
t14 = -t81 + t23;
t37 = t81 - t200;
t144 = t14 * qJD(1) - t37 * qJD(2);
t15 = t162 * t129 + t197 + t99;
t36 = t197 + t199;
t143 = t15 * qJD(1) + t36 * qJD(2);
t142 = qJD(3) * t90 - t176;
t141 = -t175 - t176;
t47 = t164 * t127 + t98;
t140 = pkin(2) * t178 + t47 * qJD(1);
t100 = t129 * t166;
t48 = t164 * t129 + t100;
t139 = pkin(2) * t177 + t48 * qJD(1);
t133 = (t153 - t121 / 0.2e1) * t215;
t11 = -t162 * t90 + t133;
t138 = t11 * qJD(1) - qJD(2) * t213;
t18 = t97 + t211;
t137 = t18 * qJD(1) + t72 * t178;
t136 = t23 * qJD(1) + t87 * t178;
t135 = t19 * qJD(1) + t72 * t177;
t123 = qJ(4) * qJD(4);
t122 = qJD(3) * qJ(4);
t115 = pkin(7) * t120;
t110 = t127 * t120;
t102 = t129 * t168;
t94 = t219 * t127;
t89 = t172 * t125;
t88 = t172 * t124;
t86 = t111 * t120;
t82 = t94 * qJD(3);
t69 = t172 * t129 * t127;
t68 = t76 * qJD(3);
t50 = t100 + (-pkin(2) / 0.2e1 + t218) * t129;
t49 = pkin(2) * t217 + t127 * t218 + t98;
t20 = t97 - t211;
t17 = t81 + t24;
t16 = -t197 - t203 / 0.2e1 - t199 / 0.2e1 + t99;
t12 = t213 / 0.2e1 + t214 / 0.2e1 + t133;
t8 = t98 + t224;
t7 = t21 - t202;
t2 = t4 / 0.2e1 + t3 / 0.2e1 + t132;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t130 * t208, 0, 0, t110, t174, 0, -t110, 0, 0, t112 * t117 - t104, t112 * t120 + t103, t65, t29 * qJD(2), 0, 0, 0, t110, t174, -t110, t65, t34 * qJD(3) + t185, -t33 * qJD(3) + t181, t22 * qJD(2) + t142 * t73, 0, 0, 0, -t110, -t174, t110, t65, -t25 * qJD(3) + t108 + t181, -t26 * qJD(3) + t158, t9 * qJD(2) + t3 * qJD(3) + t141 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 * t152, -t130 * t152, 0, 0, t110, t174, 0, -t110, 0, 0, t49 * qJD(3) - t102 - t104, t50 * qJD(3) + t101 + t103, t32, t191 + (-pkin(2) * t216 + t184) * qJD(2), 0, 0, 0, t110, t174, -t110, t32, t17 * qJD(3) + t102 + t185, t16 * qJD(3) + t157, t194 + (t87 * t216 + t184) * qJD(2) + t12 * qJD(3) + t24 * qJD(4), 0, 0, 0, -t110, -t174, t110, t32, t7 * qJD(3) + t108 + t157, t8 * qJD(3) - t102 + t158, t198 + t2 * qJD(3) + t20 * qJD(4) + t21 * qJD(5) + (t128 * t72 + (t127 * t94 + t129 * t95) * t130) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t222, t120, -t69, -t117, 0, t49 * qJD(2) + t156 - t86, t50 * qJD(2) + t154 + t155, 0, 0, 0, -t120, t117, t69, t222, -t69, t221, t17 * qJD(2) + t189 + t86, t16 * qJD(2) - t154 - t190, t12 * qJD(2) + t111 * t221 + t167, 0, t117, t120, -t69, -t222, t69, t149, t7 * qJD(2) - t193 - t68, t8 * qJD(2) - t187 - t192, t207 + t2 * qJD(2) + (-t76 * qJ(4) - t77 * t126) * qJD(3) + t77 * qJD(4) - t76 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t69, t88, t24 * qJD(2) - t159 + t86, 0, 0, 0, 0, 0, 0, t120, t88, t69, t20 * qJD(2) - t161 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t69, t89, t21 * qJD(2) - t160 - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t130 * t209, 0, 0, t110, t174, 0, -t110, 0, 0, -t47 * qJD(3) + t102, -t48 * qJD(3) - t101, -t188, -t191, 0, 0, 0, t110, t174, -t110, -t188, -t14 * qJD(3) - t102 - t109, -t15 * qJD(3) + t101 + t116, -t11 * qJD(3) - t23 * qJD(4) - t194, 0, 0, 0, -t110, -t174, t110, -t188, -t6 * qJD(3) + t101 + t183, -t5 * qJD(3) + t102 + t182, -t1 * qJD(3) - t18 * qJD(4) - t19 * qJD(5) - t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t174, 0, -t110, 0, 0, -pkin(2) * t117, -pkin(2) * t120, 0, 0, 0, 0, 0, t110, t174, -t110, 0, t37 * qJD(3) - t109, -t36 * qJD(3) + t116, t142 * t87, 0, 0, 0, -t110, -t174, t110, 0, -t27 * qJD(3) + t183, -t28 * qJD(3) + t182, t4 * qJD(3) + t141 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t222, t120, -t69, -t117, 0, -t115 - t140, -t139 + t169, 0, 0, 0, -t120, t117, t69, t222, -t69, t221, t115 - t144, -t143 - t169, pkin(7) * t221 - t138, 0, t117, t120, -t69, -t222, t69, t149, -t146 - t82, -t145 - t186, (-t94 * qJ(4) - t95 * t126) * qJD(3) + t95 * qJD(4) - t94 * qJD(5) - t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t69, t88, t115 - t136, 0, 0, 0, 0, 0, 0, t120, t88, t69, -t137 + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t69, t89, -t135 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t222, 0, t69, 0, 0, t47 * qJD(2) - t156, t48 * qJD(2) - t155, 0, 0, 0, 0, 0, -t69, -t222, t69, 0, t14 * qJD(2) - t189, t15 * qJD(2) + t190, t11 * qJD(2) - t167, 0, 0, 0, t69, t222, -t69, 0, t6 * qJD(2) + t193, t5 * qJD(2) + t192, t1 * qJD(2) - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t222, 0, t69, 0, 0, t140, t139, 0, 0, 0, 0, 0, -t69, -t222, t69, 0, t144, t143, t138, 0, 0, 0, t69, t222, -t69, 0, t146, t145, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t123, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJD(5), t126 * qJD(5) + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t122, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t88, qJD(2) * t23 + t159, 0, 0, 0, 0, 0, 0, 0, -t88, -t69, qJD(2) * t18 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t88, t136, 0, 0, 0, 0, 0, 0, 0, -t88, -t69, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t122, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, -qJD(5) - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t89, qJD(2) * t19 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t89, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), qJD(4) - t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
