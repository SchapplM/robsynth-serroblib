% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:18
% EndTime: 2019-12-31 18:43:22
% DurationCPUTime: 1.40s
% Computational Cost: add. (1760->195), mult. (3445->304), div. (0->0), fcn. (2913->6), ass. (0->168)
t117 = sin(qJ(4));
t217 = 0.2e1 * t117;
t120 = cos(qJ(3));
t118 = sin(qJ(3));
t119 = cos(qJ(4));
t190 = t118 * t119;
t160 = qJ(5) * t190;
t105 = sin(pkin(8)) * pkin(1) + pkin(6);
t195 = t105 * t117;
t106 = -cos(pkin(8)) * pkin(1) - pkin(2);
t139 = -t120 * pkin(3) - t118 * pkin(7);
t72 = t106 + t139;
t60 = t119 * t72;
t38 = -t160 + t60 + (-pkin(4) - t195) * t120;
t188 = t120 * t105;
t162 = t117 * t188;
t50 = -t60 + t162;
t44 = -t50 - t160;
t216 = -t38 + t44;
t212 = t44 / 0.2e1;
t213 = -t38 / 0.2e1;
t163 = t212 + t213;
t201 = -qJ(5) - pkin(7);
t90 = t201 * t119;
t215 = t163 * t90;
t202 = t120 * pkin(7);
t204 = t118 * pkin(3);
t91 = -t202 + t204;
t84 = t119 * t91;
t192 = t117 * t118;
t86 = t105 * t192;
t46 = -t119 * t120 * qJ(5) + t118 * pkin(4) + t84 + t86;
t211 = t46 / 0.2e1;
t112 = t117 ^ 2;
t210 = t112 / 0.2e1;
t114 = t119 ^ 2;
t209 = -t114 / 0.2e1;
t208 = -t117 / 0.2e1;
t207 = -t118 / 0.2e1;
t206 = t119 / 0.2e1;
t205 = t117 * pkin(4);
t203 = t120 * pkin(4);
t200 = pkin(4) * qJD(4);
t164 = t203 / 0.2e1;
t5 = (t164 + t163) * t119;
t199 = t5 * qJD(1);
t144 = t105 + t205;
t63 = t144 * t118;
t198 = t63 * t117;
t165 = -t203 / 0.2e1;
t140 = t213 + t165;
t7 = (t212 + t140) * t190;
t197 = t7 * qJD(1);
t10 = t216 * t192;
t196 = t10 * qJD(1);
t158 = t119 * t188;
t51 = t117 * t72 + t158;
t45 = -qJ(5) * t192 + t51;
t161 = t105 * t190;
t191 = t117 * t120;
t83 = t117 * t91;
t49 = -qJ(5) * t191 - t161 + t83;
t11 = (t46 * t118 + t38 * t120) * t119 + (t49 * t118 + t45 * t120) * t117;
t194 = t11 * qJD(1);
t193 = t112 * t120;
t113 = t118 ^ 2;
t189 = t119 * t113;
t14 = (t45 * t117 + t38 * t119) * t118;
t187 = t14 * qJD(1);
t15 = t50 * t118 + (-t86 + t84) * t120;
t186 = t15 * qJD(1);
t16 = t83 * t120 + (-t51 + t158) * t118;
t185 = t16 * qJD(1);
t25 = -t113 * t195 - t50 * t120;
t184 = t25 * qJD(1);
t26 = -t105 * t189 - t51 * t120;
t183 = t26 * qJD(1);
t98 = t114 + t112;
t79 = t98 * t113;
t182 = t79 * qJD(1);
t115 = t120 ^ 2;
t100 = t115 - t113;
t81 = t100 * t117;
t181 = t81 * qJD(1);
t82 = t119 * t115 - t189;
t180 = t82 * qJD(1);
t179 = t98 * qJD(3);
t99 = t114 - t112;
t178 = qJD(3) * t117;
t177 = qJD(3) * t119;
t176 = qJD(4) * t117;
t175 = qJD(4) * t119;
t174 = qJD(4) * t120;
t173 = t100 * qJD(1);
t172 = t118 * qJD(1);
t171 = t118 * qJD(3);
t170 = t118 * qJD(4);
t169 = t120 * qJD(1);
t168 = t120 * qJD(3);
t167 = pkin(4) * t190;
t166 = t90 * t192;
t110 = -t119 * pkin(4) - pkin(3);
t159 = t110 * t190;
t157 = t38 * t208;
t156 = t119 * t172;
t155 = t117 * t170;
t154 = t117 * t174;
t153 = t119 * t170;
t152 = t119 * t174;
t151 = t106 * t172;
t150 = t106 * t169;
t149 = t117 * t175;
t148 = t117 * t177;
t147 = t118 * t168;
t146 = t118 * t169;
t145 = t119 * t171;
t143 = pkin(4) * t153;
t142 = -qJD(4) + t169;
t141 = t117 * t145;
t64 = t144 * t120;
t4 = (t45 * t206 + t157 - t64 / 0.2e1) * t120 + (t49 * t206 + t46 * t208 + t63 / 0.2e1) * t118;
t8 = t38 * t46 + t45 * t49 + t63 * t64;
t138 = t8 * qJD(1) + t4 * qJD(2);
t9 = t63 * t167 + t216 * t45;
t137 = t9 * qJD(1) + t7 * qJD(2);
t89 = t201 * t117;
t55 = -t89 * t117 - t90 * t119;
t56 = (-0.1e1 + t98) * t120 * t118;
t136 = t4 * qJD(1) + t56 * qJD(2);
t18 = t110 * t205;
t2 = t215 + (t211 - t159 / 0.2e1 - t198 / 0.2e1) * pkin(4);
t135 = -t2 * qJD(1) + t18 * qJD(3);
t134 = t142 * t118;
t133 = t202 / 0.2e1 - t204 / 0.2e1;
t126 = t133 * t119;
t53 = -t84 / 0.2e1 + t126;
t132 = pkin(3) * t178 - t53 * qJD(1);
t127 = t133 * t117;
t52 = t83 / 0.2e1 - t127;
t131 = pkin(3) * t177 - t52 * qJD(1);
t74 = (t210 + t209) * t118;
t130 = -t74 * qJD(1) + t148;
t129 = t119 * t134;
t128 = t156 + t178;
t125 = t117 * qJD(1) * t189 + t74 * qJD(3);
t80 = t99 * t113;
t124 = t80 * qJD(1) + 0.2e1 * t141;
t123 = -t99 * qJD(3) + t156 * t217;
t68 = t166 / 0.2e1;
t122 = t68 + (t45 / 0.2e1 + t89 * t207) * t119;
t13 = -t188 / 0.2e1 + t140 * t117 + t122;
t61 = t112 * t207 + (0.1e1 / 0.2e1 + t209) * t118;
t121 = t13 * qJD(1) - t61 * qJD(2) + t55 * qJD(3);
t108 = t114 * t120;
t107 = t171 / 0.2e1;
t103 = t117 * t171;
t78 = (t169 - qJD(4) / 0.2e1) * t118;
t73 = t128 * pkin(4);
t70 = t74 * qJD(4);
t62 = (t114 / 0.2e1 + t210 + 0.1e1 / 0.2e1) * t118;
t41 = -pkin(4) * t191 - t166 / 0.2e1 + t68;
t40 = t86 + t84 / 0.2e1 + t126;
t39 = t161 - t83 / 0.2e1 - t127;
t12 = t157 + t188 / 0.2e1 + t117 * t164 + t122;
t6 = (t163 + t165) * t119;
t3 = pkin(4) * t211 - t215 + (t159 + t198) * pkin(4) / 0.2e1;
t1 = t4 * qJD(3) + t7 * qJD(4);
t17 = [0, 0, 0, 0, t147, t100 * qJD(3), 0, 0, 0, t106 * t171, t106 * t168, -t113 * t149 + t114 * t147, -t80 * qJD(4) - 0.2e1 * t120 * t141, -t82 * qJD(3) + t118 * t154, t81 * qJD(3) + t118 * t152, -t147, -t15 * qJD(3) - t26 * qJD(4), t16 * qJD(3) + t25 * qJD(4), -t11 * qJD(3) - t10 * qJD(4) + t79 * qJD(5), t8 * qJD(3) + t9 * qJD(4) - t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, t146, t173, t168, -t171, 0, -t105 * t168 + t151, t105 * t171 + t150, -t70 + (t114 * t172 + t148) * t120, (t108 - t193) * qJD(3) + (-qJD(4) - t169) * t190 * t217, t103 - t180, t145 + t181, -t78, -t186 + (t139 * t117 - t158) * qJD(3) + t40 * qJD(4), t185 + (t139 * t119 + t162) * qJD(3) + t39 * qJD(4), -t194 + ((-t120 * t89 + t49) * t119 + (t120 * t90 - t46) * t117) * qJD(3) + t6 * qJD(4), (t64 * t110 + t46 * t89 - t49 * t90) * qJD(3) + t3 * qJD(4) + t12 * qJD(5) + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124, t117 * t134, t129, t107, t40 * qJD(3) - t51 * qJD(4) - t183, t39 * qJD(3) + t50 * qJD(4) + t184, pkin(4) * t155 + t6 * qJD(3) - t196, t3 * qJD(3) - t200 * t45 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t12 * qJD(3) - t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t168, 0, 0, 0, 0, 0, -t145 - t154, t103 - t152, (t108 + t193) * qJD(3), (t118 * t110 + t120 * t55) * qJD(3) + t41 * qJD(4) + t62 * qJD(5) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 * t168 - t153, -t119 * t168 + t155, 0, t41 * qJD(3) - t143 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * qJD(3); 0, 0, 0, 0, -t146, -t173, 0, 0, 0, -t151, -t150, -t114 * t146 - t70, t129 * t217, -t152 + t180, t154 - t181, t78, t53 * qJD(4) + t186, t52 * qJD(4) - t185, t5 * qJD(4) + t194, -t2 * qJD(4) + t13 * qJD(5) - t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * qJD(5) - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t99 * qJD(4), 0, 0, 0, -pkin(3) * t176, -pkin(3) * t175, t98 * qJD(5), t18 * qJD(4) + t55 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t123, -t142 * t119, t142 * t117, -t172 / 0.2e1, -pkin(7) * t175 - t132, pkin(7) * t176 - t131, -pkin(4) * t175 + t199, t200 * t90 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t124, (-t117 * t172 + t177) * t120, -t128 * t120, t107, -t53 * qJD(3) + t183, -t52 * qJD(3) - t184, -t5 * qJD(3) + t196, t2 * qJD(3) - qJD(5) * t167 - t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t123, t119 * t169, -t117 * t169, t172 / 0.2e1, t132, t131, -t199, -qJD(5) * t205 - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, -t13 * qJD(3) + t143 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, pkin(4) * t176 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t17;
