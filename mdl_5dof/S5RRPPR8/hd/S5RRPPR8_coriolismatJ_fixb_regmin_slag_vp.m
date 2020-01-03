% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:18
% EndTime: 2019-12-31 19:39:23
% DurationCPUTime: 2.11s
% Computational Cost: add. (1446->164), mult. (2807->241), div. (0->0), fcn. (3018->6), ass. (0->142)
t163 = qJD(2) - qJD(5);
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t129 = sin(qJ(2));
t232 = pkin(6) - qJ(4);
t108 = t232 * t129;
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t131 = cos(qJ(2));
t140 = t232 * t131;
t68 = t126 * t108 + t127 * t140;
t95 = t129 * t126 + t131 * t127;
t237 = -t95 * pkin(7) + t68;
t67 = -t127 * t108 + t126 * t140;
t97 = -t131 * t126 + t129 * t127;
t36 = -t97 * pkin(7) - t67;
t248 = t163 * (t128 * t36 + t130 * t237);
t247 = t163 * (-t128 * t237 + t130 * t36);
t210 = t130 * t97;
t216 = t128 * t95;
t137 = -t210 + t216;
t215 = t128 * t97;
t85 = t130 * t95;
t230 = t85 + t215;
t7 = t137 ^ 2 - t230 ^ 2;
t242 = t7 * qJD(1);
t224 = t126 / 0.2e1;
t235 = t137 * qJD(1);
t236 = t230 * t235;
t186 = qJD(5) * t137;
t147 = -t216 / 0.2e1;
t73 = t216 / 0.2e1;
t231 = t147 + t73;
t234 = t231 * qJD(2);
t164 = t131 * qJD(3);
t174 = -t131 * pkin(2) - t129 * qJ(3);
t233 = t174 * qJD(2) + t164;
t145 = t215 / 0.2e1;
t146 = -t215 / 0.2e1;
t43 = t146 + t145;
t99 = t130 * t126 + t128 * t127;
t177 = t99 * qJD(2);
t229 = -t99 * qJD(5) + t177;
t94 = t128 * t126 - t130 * t127;
t182 = t94 * qJD(2);
t228 = -t94 * qJD(5) + t182;
t227 = t68 / 0.2e1;
t226 = t85 / 0.2e1;
t132 = -pkin(2) - pkin(3);
t223 = t127 / 0.2e1;
t222 = t129 / 0.2e1;
t107 = -pkin(1) + t174;
t91 = t131 * pkin(3) - t107;
t121 = t131 * qJ(3);
t92 = t132 * t129 + t121;
t6 = t91 * t92;
t207 = t6 * qJD(1);
t58 = t95 * pkin(4) + t91;
t60 = -t97 * pkin(4) + t92;
t8 = t137 * t58 + t230 * t60;
t205 = t8 * qJD(1);
t9 = -t137 * t60 + t230 * t58;
t204 = t9 * qJD(1);
t202 = qJD(1) * t58;
t13 = t67 * t97 - t68 * t95;
t199 = t13 * qJD(1);
t14 = 0.2e1 * t73 - t210;
t198 = t14 * qJD(1);
t139 = 0.2e1 * t226;
t16 = 0.2e1 * t145 + t139;
t197 = t16 * qJD(1);
t133 = (-pkin(2) / 0.2e1 - pkin(3) / 0.2e1) * t129 + t121 / 0.2e1;
t105 = -t126 * qJ(3) + t127 * t132;
t106 = t127 * qJ(3) + t126 * t132;
t135 = t105 * t97 / 0.2e1 + t106 * t95 / 0.2e1;
t18 = t133 + t135;
t196 = t18 * qJD(1);
t20 = -t91 * t97 + t92 * t95;
t195 = t20 * qJD(1);
t21 = t91 * t95 + t92 * t97;
t194 = t21 * qJD(1);
t24 = 0.2e1 * t146 - t85;
t192 = t24 * qJD(1);
t134 = t97 * t223 + t95 * t224;
t39 = t222 + t134;
t191 = t39 * qJD(1);
t190 = t231 * qJD(1);
t44 = t226 - t85 / 0.2e1;
t189 = t44 * qJD(1);
t188 = t230 * qJD(5);
t187 = t230 * qJD(2);
t59 = t95 ^ 2 + t97 ^ 2;
t185 = t59 * qJD(1);
t109 = t129 * pkin(2) - t121;
t61 = t107 * t131 + t109 * t129;
t184 = t61 * qJD(1);
t62 = -t107 * t129 + t109 * t131;
t183 = t62 * qJD(1);
t181 = t94 * qJD(3);
t179 = t95 * qJD(1);
t178 = t97 * qJD(1);
t176 = t99 * qJD(3);
t173 = qJD(1) * t129;
t172 = qJD(1) * t131;
t171 = qJD(2) * qJ(3);
t170 = qJD(3) * t129;
t125 = t129 ^ 2;
t110 = t131 ^ 2 - t125;
t169 = t110 * qJD(1);
t168 = t125 * qJD(1);
t167 = t126 * qJD(2);
t166 = t127 * qJD(2);
t165 = t129 * qJD(2);
t118 = t131 * qJD(2);
t162 = pkin(1) * t173;
t161 = pkin(1) * t172;
t160 = pkin(6) * t165;
t159 = pkin(6) * t118;
t156 = t230 * t202;
t155 = t137 * t202;
t152 = t230 * t173;
t151 = t137 * t173;
t150 = t95 * t173;
t149 = t97 * t173;
t148 = t91 * t173;
t143 = t107 * t109 * qJD(1);
t142 = t107 * t173;
t138 = -pkin(4) + t105;
t10 = (t227 - t68 / 0.2e1) * t127;
t57 = -t105 * t126 + t106 * t127;
t136 = t10 * qJD(1) - t57 * qJD(2);
t111 = t129 * t172;
t56 = t130 * t106 + t128 * t138;
t55 = t128 * t106 - t130 * t138;
t40 = t222 - t134;
t26 = 0.2e1 * t147 + t210;
t25 = t139 + t215;
t22 = t126 * t97 - t127 * t95;
t19 = t133 - t135;
t17 = t44 + t43;
t11 = t127 * t227 + t68 * t223 + 0.2e1 * t67 * t224;
t1 = [0, 0, 0, t129 * t118, t110 * qJD(2), 0, 0, 0, -pkin(1) * t165, -pkin(1) * t118, -t62 * qJD(2) + t129 * t164, 0, -t61 * qJD(2) + t125 * qJD(3), (qJD(2) * t109 - t170) * t107, t20 * qJD(2) + t95 * t170, t21 * qJD(2) + t97 * t170, t59 * qJD(4), t6 * qJD(2) + t13 * qJD(4) + t170 * t91, -(t187 - t188) * t137, t163 * t7, 0, 0, 0, t8 * qJD(2) + t170 * t230 - t186 * t58, t9 * qJD(2) - t137 * t170 - t188 * t58; 0, 0, 0, t111, t169, t118, -t165, 0, -t159 - t162, t160 - t161, -t159 - t183, t233, -t160 - t184, t233 * pkin(6) + t143, -qJD(2) * t68 + t195, t67 * qJD(2) + t194, (-t105 * t95 + t106 * t97) * qJD(2) + t22 * qJD(3), t207 + (t105 * t68 + t67 * t106) * qJD(2) + t11 * qJD(3) + t19 * qJD(4), -t236, t242, t25 * qJD(5) - t187, qJD(2) * t137 + t26 * qJD(5), 0, t231 * qJD(4) + t205 - t248, t17 * qJD(4) + t204 - t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t118, t168, -t142 + t159, t150, t149, t22 * qJD(2), t11 * qJD(2) + t40 * qJD(4) + t148, 0, 0, 0, 0, 0, t152, -t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t19 * qJD(2) + t40 * qJD(3) + t199, 0, 0, 0, 0, 0, t234, t17 * qJD(2) + t43 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t242, t25 * qJD(2) - t188, t26 * qJD(2) + t186, 0, -t155 + t248, t43 * qJD(4) - t156 + t247; 0, 0, 0, -t111, -t169, 0, 0, 0, t162, t161, t183, 0, t184, -t143, t97 * qJD(4) - t195, -t95 * qJD(4) - t194, 0, -t10 * qJD(3) - t18 * qJD(4) - t207, t236, -t242, -t44 * qJD(5), -t231 * qJD(5), 0, -t14 * qJD(4) - t205, -t16 * qJD(4) - t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), t126 * qJD(3), t127 * qJD(3), 0, t57 * qJD(3), 0, 0, 0, 0, 0, t56 * qJD(5) + t176, -t55 * qJD(5) - t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t171, t167, t166, 0, -t136, 0, 0, 0, 0, 0, t177, -t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, -t179, 0, -t196, 0, 0, 0, 0, 0, -t198, -t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t190, 0, t163 * t56, -t163 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, 0, -t168, t142, -t150, -t149, 0, t10 * qJD(2) - t39 * qJD(4) - t148, 0, 0, 0, 0, 0, -t152, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t171, -t167, -t166, 0, t136, 0, 0, 0, 0, 0, -t229, t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, -t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * qJD(2), t95 * qJD(2), -t185, t18 * qJD(2) + t39 * qJD(3) - t199, 0, 0, 0, 0, 0, t14 * qJD(2) - t186, t16 * qJD(2) + t24 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t179, 0, t196, 0, 0, 0, 0, 0, t198, t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t242, t44 * qJD(2), t234, 0, qJD(4) * t137 + t155, -t24 * qJD(4) + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t190, 0, -t56 * qJD(2) - t176, t55 * qJD(2) + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, -t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
