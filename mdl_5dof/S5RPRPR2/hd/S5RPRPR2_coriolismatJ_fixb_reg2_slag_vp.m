% Calculate inertial parameters regressor of coriolis matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:12
% EndTime: 2020-01-03 11:34:16
% DurationCPUTime: 1.41s
% Computational Cost: add. (2910->136), mult. (5522->170), div. (0->0), fcn. (5428->8), ass. (0->117)
t134 = cos(pkin(9));
t196 = cos(qJ(5));
t161 = t196 * t134;
t132 = sin(pkin(9));
t135 = sin(qJ(5));
t182 = t135 * t132;
t115 = -t161 + t182;
t162 = t196 * t132;
t181 = t135 * t134;
t117 = t162 + t181;
t166 = qJD(1) + qJD(3);
t205 = t166 * t117;
t209 = t115 * t205;
t202 = t117 ^ 2;
t203 = t115 ^ 2;
t74 = -t202 + t203;
t208 = t166 * t74;
t81 = t202 + t203;
t207 = t166 * t81;
t206 = t166 * t115;
t130 = t132 ^ 2;
t131 = t134 ^ 2;
t174 = -t131 - t130;
t204 = t166 * t174;
t119 = t174 * qJ(4);
t200 = -t115 / 0.2e1;
t199 = t117 / 0.2e1;
t194 = t134 * pkin(4);
t128 = -pkin(3) - t194;
t198 = t128 / 0.2e1;
t197 = cos(qJ(3));
t195 = sin(pkin(8)) * pkin(1);
t136 = sin(qJ(3));
t152 = cos(pkin(8)) * pkin(1) + pkin(2);
t105 = t136 * t195 - t197 * t152;
t60 = t117 * t105;
t61 = t115 * t105;
t18 = -t61 * t115 - t60 * t117;
t80 = t81 * qJD(4);
t192 = t18 * qJD(3) + t80;
t129 = t134 * pkin(7);
t120 = t134 * qJ(4) + t129;
t156 = (-pkin(7) - qJ(4)) * t132;
t88 = t135 * t120 - t196 * t156;
t189 = t117 * t88;
t146 = t136 * t152;
t153 = t197 * t195;
t106 = t153 + t146;
t141 = qJ(4) + t106;
t137 = (-pkin(7) - t141) * t132;
t90 = t134 * t141 + t129;
t49 = t135 * t90 - t196 * t137;
t188 = t49 * t117;
t50 = t135 * t137 + t196 * t90;
t187 = t50 * t115;
t89 = t196 * t120 + t135 * t156;
t185 = t89 * t115;
t121 = t174 * qJD(4);
t73 = t174 * t105;
t184 = t73 * qJD(3) - t121;
t15 = -t187 + t188;
t180 = t15 * qJD(1);
t17 = t61 * t199 + t60 * t200;
t179 = t17 * qJD(1);
t178 = t18 * qJD(1);
t148 = -pkin(3) + t105;
t138 = t131 * t141;
t139 = t130 * t141;
t62 = t138 + t139;
t19 = -t62 * t105 + t148 * t106;
t177 = t19 * qJD(1);
t176 = t62 * qJD(1);
t175 = t73 * qJD(1);
t173 = t105 * qJD(1);
t172 = t106 * qJD(1);
t104 = t106 * qJD(3);
t171 = t115 * qJD(1);
t170 = t115 * qJD(3);
t110 = t115 * qJD(5);
t169 = t117 * qJD(1);
t168 = t117 * qJD(3);
t167 = t117 * qJD(5);
t97 = t148 - t194;
t165 = t97 * t171;
t164 = t97 * t169;
t163 = t198 + t97 / 0.2e1;
t160 = t106 * t171;
t159 = t106 * t169;
t158 = t132 * t172;
t157 = t115 * t167;
t8 = t97 * t106 - t49 * t60 + t50 * t61;
t151 = t8 * qJD(1) + t17 * qJD(2);
t35 = -t185 + t189;
t140 = t153 / 0.2e1 + t146 / 0.2e1;
t9 = (-t88 / 0.2e1 - t49 / 0.2e1) * t117 + (t89 / 0.2e1 + t50 / 0.2e1) * t115 + t140;
t149 = -t9 * qJD(1) + t35 * qJD(3);
t44 = t119 + t140 + (-t131 / 0.2e1 - t130 / 0.2e1) * t106;
t147 = t44 * qJD(1) + t119 * qJD(3);
t142 = (t161 / 0.2e1 - t182 / 0.2e1) * t105;
t22 = t163 * t115 + t142;
t145 = t22 * qJD(1) + t128 * t170;
t143 = (t181 / 0.2e1 + t162 / 0.2e1) * t105;
t21 = -t163 * t117 + t143;
t144 = t21 * qJD(1) - t128 * t168;
t113 = t117 * qJD(4);
t109 = t115 * qJD(4);
t103 = t105 * qJD(3);
t98 = t132 * t104;
t78 = t106 * t168;
t77 = t106 * t170;
t64 = t74 * qJD(5);
t45 = t138 / 0.2e1 + t139 / 0.2e1 + t140 - t119 / 0.2e1;
t24 = t117 * t198 + t97 * t199 + t143;
t23 = t142 + (t128 + t97) * t200;
t10 = t140 + t189 / 0.2e1 - t185 / 0.2e1 + t188 / 0.2e1 - t187 / 0.2e1;
t1 = t17 * qJD(3);
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t103, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t104, t98, t184, t19 * qJD(3) + t62 * qJD(4), -t157, t64, 0, t157, 0, 0, t97 * t167 + t77, -t97 * t110 + t78, t192, t8 * qJD(3) + t15 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 - t172, t103 + t173, 0, 0, 0, 0, 0, 0, 0, 0, -t166 * t134 * t106, t98 + t158, t175 + t184, t177 + (-t106 * pkin(3) + qJ(4) * t73) * qJD(3) + t45 * qJD(4), -t157, t64, 0, t157, 0, 0, t24 * qJD(5) + t160 + t77, t23 * qJD(5) + t159 + t78, t178 + t192, (t106 * t128 - t60 * t88 + t61 * t89) * qJD(3) + t10 * qJD(4) + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t45 * qJD(3) + t176, 0, 0, 0, 0, 0, 0, 0, 0, t207, t10 * qJD(3) + t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t208, -t110, t209, -t167, 0, t24 * qJD(3) - t50 * qJD(5) + t164, t23 * qJD(3) + t49 * qJD(5) - t165, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t110, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t173, 0, 0, 0, 0, 0, 0, 0, 0, t134 * t172, -t158, -t121 - t175, -t44 * qJD(4) - t177, -t157, t64, 0, t157, 0, 0, -t21 * qJD(5) - t160, -t22 * qJD(5) - t159, t80 - t178, -t9 * qJD(4) - t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t119 * qJD(4), -t157, t64, 0, t157, 0, 0, t128 * t167, -t128 * t110, t80, t35 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, -t147, 0, 0, 0, 0, 0, 0, 0, 0, t207, t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t208, -t110, t209, -t167, 0, -t89 * qJD(5) - t144, t88 * qJD(5) - t145, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, t44 * qJD(3) - t176, 0, 0, 0, 0, 0, 0, t167, -t110, -t207, t9 * qJD(3) - t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, t147, 0, 0, 0, 0, 0, 0, t167, -t110, -t207, -t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, -t206, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -t208, 0, -t209, 0, 0, t21 * qJD(3) - t113 - t164, t22 * qJD(3) + t109 + t165, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -t208, 0, -t209, 0, 0, -t113 + t144, t109 + t145, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t206, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;