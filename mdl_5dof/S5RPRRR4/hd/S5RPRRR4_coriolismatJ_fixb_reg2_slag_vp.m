% Calculate inertial parameters regressor of coriolis matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:49
% EndTime: 2022-01-23 09:34:51
% DurationCPUTime: 1.97s
% Computational Cost: add. (2365->192), mult. (4607->247), div. (0->0), fcn. (3840->8), ass. (0->153)
t115 = sin(qJ(4));
t137 = cos(pkin(9)) * pkin(1) + pkin(2);
t187 = cos(qJ(3));
t123 = t187 * t137;
t144 = sin(pkin(9)) * pkin(1);
t186 = sin(qJ(3));
t130 = t186 * t144;
t89 = t130 - t123;
t179 = t115 * t89;
t117 = cos(qJ(4));
t90 = t137 * t186 + t144 * t187;
t86 = t117 * t90;
t67 = t86 - t179;
t204 = -t67 / 0.2e1;
t114 = sin(qJ(5));
t112 = t114 ^ 2;
t116 = cos(qJ(5));
t113 = t116 ^ 2;
t143 = t113 / 0.2e1 + t112 / 0.2e1;
t177 = t117 * t89;
t178 = t115 * t90;
t68 = -t177 - t178;
t200 = t143 * t68;
t203 = pkin(4) * t204 + t200 * pkin(8);
t103 = t113 - t112;
t154 = -qJD(3) - qJD(4);
t142 = qJD(1) - t154;
t201 = t142 * t103;
t160 = t112 + t113;
t199 = pkin(4) / 0.2e1;
t120 = pkin(3) - t89;
t85 = t117 * t120;
t65 = -t85 + t178;
t63 = -pkin(4) + t65;
t198 = -t63 / 0.2e1;
t119 = t115 * t120;
t66 = t119 + t86;
t197 = t66 / 0.2e1;
t196 = -t68 / 0.2e1;
t195 = t89 / 0.2e1;
t183 = t117 * pkin(3);
t106 = -pkin(4) - t183;
t193 = -t106 / 0.2e1;
t190 = -t114 / 0.2e1;
t189 = t114 / 0.2e1;
t188 = t116 / 0.2e1;
t185 = pkin(4) * t116;
t184 = t115 * pkin(3);
t60 = t66 * qJD(4);
t61 = t67 * qJD(3);
t182 = -t61 - t60;
t181 = pkin(3) * qJD(3);
t180 = pkin(3) * qJD(4);
t25 = t160 * t65;
t64 = pkin(8) + t66;
t3 = -t25 * t64 + t63 * t66;
t176 = t3 * qJD(1);
t30 = t160 * t68;
t4 = t30 * t64 + t63 * t67;
t175 = t4 * qJD(1);
t174 = t67 * t116;
t8 = t65 * t67 + t66 * t68;
t173 = t8 * qJD(1);
t172 = qJD(1) * t25;
t171 = qJD(1) * t30;
t170 = t106 * t116;
t138 = t184 / 0.2e1 + t197;
t129 = t204 + t138;
t15 = t129 * t116;
t169 = t15 * qJD(1);
t26 = (t195 + pkin(3) - t130 / 0.2e1 + t123 / 0.2e1) * t115;
t168 = t26 * qJD(1);
t28 = t85 / 0.2e1 + (t195 + pkin(3) / 0.2e1) * t117;
t167 = t28 * qJD(1);
t166 = t65 * qJD(1);
t165 = t66 * qJD(1);
t164 = t67 * qJD(1);
t163 = t68 * qJD(1);
t162 = t89 * qJD(1);
t161 = t90 * qJD(1);
t159 = qJD(1) * t114;
t158 = qJD(1) * t116;
t157 = qJD(3) * t114;
t156 = qJD(4) * t114;
t155 = t114 * qJD(5);
t111 = t116 * qJD(5);
t153 = t115 * t180;
t152 = t115 * t181;
t151 = -t183 / 0.2e1;
t150 = t63 * t159;
t149 = t63 * t158;
t148 = t66 * t159;
t147 = t67 * t159;
t141 = pkin(3) * t154;
t140 = t160 * t117;
t139 = t65 / 0.2e1 + t199 + t198;
t136 = t196 + t193 + t198;
t135 = t143 * t65;
t134 = t115 * t141;
t133 = t143 * t117;
t105 = pkin(8) + t184;
t118 = (t63 * t115 / 0.2e1 + t64 * t133) * pkin(3) - t105 * t135 + t106 * t197;
t2 = t118 - t203;
t69 = (t105 * t140 + t106 * t115) * pkin(3);
t132 = -t2 * qJD(1) - t69 * qJD(3);
t6 = t160 * (t183 / 0.2e1 - t65 / 0.2e1 + t196);
t92 = pkin(3) * t140;
t131 = -qJD(1) * t6 - qJD(3) * t92;
t128 = t151 + t199 + t193;
t9 = t136 * t114;
t126 = qJD(1) * t9 - t106 * t157;
t10 = t136 * t116;
t125 = qJD(1) * t10 - qJD(3) * t170;
t14 = t129 * t114;
t124 = -qJD(1) * t14 - t114 * t152;
t17 = t139 * t114;
t72 = t128 * t114;
t122 = pkin(4) * t156 + qJD(1) * t17 + qJD(3) * t72;
t18 = t139 * t116;
t73 = t128 * t116;
t121 = qJD(1) * t18 + qJD(3) * t73 + qJD(4) * t185;
t110 = -t185 / 0.2e1;
t109 = pkin(4) * t190;
t104 = t114 * t111;
t102 = t114 * t153;
t95 = t103 * qJD(5);
t94 = t170 / 0.2e1;
t93 = t106 * t189;
t91 = t92 * qJD(4);
t88 = t90 * qJD(3);
t87 = t89 * qJD(3);
t78 = t142 * t116 * t114;
t75 = t116 * t151 + t110 + t94;
t74 = t114 * t151 + t109 + t93;
t62 = t68 * qJD(3);
t59 = t65 * qJD(4);
t57 = t67 * t157;
t50 = t66 * t156;
t41 = t63 * t188;
t40 = t63 * t189;
t29 = t151 + t178 - t85 / 0.2e1 + t177 / 0.2e1;
t27 = -t184 / 0.2e1 - t86 - t119 / 0.2e1 + t179 / 0.2e1;
t24 = t30 * qJD(3);
t23 = t25 * qJD(4);
t20 = t188 * t65 + t110 + t41;
t19 = t189 * t65 + t109 + t40;
t16 = -t174 / 0.2e1 - t138 * t116;
t13 = t114 * t138 + t189 * t67;
t12 = t116 * t196 + t41 + t94;
t11 = t190 * t68 + t40 + t93;
t5 = pkin(3) * t133 - t135 + t200;
t1 = t118 + t203;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t87, 0, 0, 0, 0, 0, 0, 0, 0, t182, -t62 + t59, 0, qJD(3) * t8, t104, t95, 0, -t104, 0, 0, t116 * t182 + t155 * t63, t111 * t63 + t50 + t57, t24 - t23, qJD(3) * t4 + qJD(4) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 - t161, t87 + t162, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t27 - t164 - t61, qJD(4) * t29 - t163 - t62, 0, t173 + (t115 * t68 - t117 * t67) * t181, t104, t95, 0, -t104, 0, 0, t16 * qJD(4) + t11 * qJD(5) + (-qJD(1) - qJD(3)) * t174, qJD(4) * t13 + qJD(5) * t12 + t147 + t57, qJD(4) * t5 + t171 + t24, t175 + (t105 * t30 + t67 * t106) * qJD(3) + t1 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t27 - t165 - t60, qJD(3) * t29 + t166 + t59, 0, 0, t104, t95, 0, -t104, 0, 0, t16 * qJD(3) + t19 * qJD(5) + (-qJD(1) - qJD(4)) * t66 * t116, qJD(3) * t13 + qJD(5) * t20 + t148 + t50, qJD(3) * t5 - t172 - t23, t176 + t1 * qJD(3) + (-t66 * pkin(4) - pkin(8) * t25) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t201, t111, -t78, -t155, 0, qJD(3) * t11 + qJD(4) * t19 - t111 * t64 + t150, qJD(3) * t12 + qJD(4) * t20 + t155 * t64 + t149, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t111, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t162, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t26 + t164, -qJD(4) * t28 + t163, 0, -t173, t104, t95, 0, -t104, 0, 0, -qJD(4) * t15 - qJD(5) * t9 + t158 * t67, qJD(4) * t14 - qJD(5) * t10 - t147, qJD(4) * t6 - t171, qJD(4) * t2 - t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t117 * t180, 0, 0, t104, t95, 0, -t104, 0, 0, t106 * t155 - t116 * t153, t106 * t111 + t102, t91, t69 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 - t168, t117 * t141 - t167, 0, 0, t104, t95, 0, -t104, 0, 0, t74 * qJD(5) + t116 * t134 - t169, qJD(5) * t75 + t102 - t124, -t131 + t91, (-pkin(4) * t115 + pkin(8) * t140) * t180 - t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t201, t111, -t78, -t155, 0, qJD(4) * t74 - t105 * t111 - t126, qJD(4) * t75 + t105 * t155 - t125, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t26 + t165, qJD(3) * t28 - t166, 0, 0, t104, t95, 0, -t104, 0, 0, qJD(3) * t15 - qJD(5) * t17 + t158 * t66, -qJD(3) * t14 - qJD(5) * t18 - t148, -qJD(3) * t6 + t172, -qJD(3) * t2 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152 + t168, t117 * t181 + t167, 0, 0, t104, t95, 0, -t104, 0, 0, -qJD(5) * t72 + t116 * t152 + t169, -qJD(5) * t73 + t124, t131, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t95, 0, -t104, 0, 0, -pkin(4) * t155, -pkin(4) * t111, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t201, t111, -t78, -t155, 0, -pkin(8) * t111 - t122, pkin(8) * t155 - t121, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t201, 0, t78, 0, 0, qJD(3) * t9 + qJD(4) * t17 - t150, qJD(3) * t10 + qJD(4) * t18 - t149, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t201, 0, t78, 0, 0, qJD(4) * t72 + t126, qJD(4) * t73 + t125, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t201, 0, t78, 0, 0, t122, t121, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
