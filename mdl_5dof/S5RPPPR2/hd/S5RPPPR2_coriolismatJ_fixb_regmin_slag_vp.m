% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:48
% EndTime: 2022-01-23 08:59:50
% DurationCPUTime: 1.17s
% Computational Cost: add. (1040->179), mult. (2810->307), div. (0->0), fcn. (3066->8), ass. (0->170)
t117 = sin(pkin(9));
t120 = cos(pkin(9));
t122 = cos(pkin(7));
t188 = t122 * t120;
t119 = sin(pkin(7));
t121 = cos(pkin(8));
t191 = t121 * t119;
t80 = t117 * t191 + t188;
t78 = t80 ^ 2;
t217 = t80 / 0.2e1;
t216 = -t117 / 0.2e1;
t215 = t117 / 0.2e1;
t214 = -t120 / 0.2e1;
t123 = sin(qJ(5));
t118 = sin(pkin(8));
t124 = cos(qJ(5));
t195 = t118 * t124;
t189 = t122 * t117;
t82 = t120 * t191 - t189;
t55 = -t119 * t195 + t123 * t82;
t196 = t118 * t123;
t57 = t119 * t196 + t82 * t124;
t213 = t55 * t57;
t190 = t121 * t122;
t93 = -pkin(2) * t122 - t119 * qJ(3) - pkin(1);
t66 = qJ(2) * t190 + t118 * t93;
t58 = -t122 * qJ(4) + t66;
t75 = (pkin(3) * t118 - qJ(4) * t121 + qJ(2)) * t119;
t33 = t117 * t75 + t120 * t58;
t198 = t118 * t119;
t32 = -t117 * t58 + t120 * t75;
t29 = -pkin(4) * t198 - t32;
t197 = t118 * t122;
t65 = -qJ(2) * t197 + t121 * t93;
t133 = t122 * pkin(3) - t65;
t28 = t80 * pkin(4) - t82 * pkin(6) + t133;
t30 = pkin(6) * t198 + t33;
t6 = t123 * t30 - t124 * t28;
t1 = -t29 * t55 + t6 * t80;
t212 = t1 * qJD(1);
t7 = -t123 * t28 - t124 * t30;
t2 = t29 * t57 + t7 * t80;
t211 = t2 * qJD(1);
t192 = t120 * t119;
t79 = t121 * t189 - t192;
t194 = t119 * t117;
t84 = t121 * t188 + t194;
t3 = t133 * t197 - t32 * t79 + t33 * t84;
t210 = t3 * qJD(1);
t149 = t118 * t192;
t150 = t118 * t194;
t4 = -t133 * t191 + t33 * t149 - t32 * t150;
t209 = t4 * qJD(1);
t5 = -t32 * t82 - t33 * t80;
t208 = t5 * qJD(1);
t146 = t55 * t216;
t148 = t120 * t196;
t185 = t124 * t121;
t85 = t148 + t185;
t152 = t85 * t217;
t186 = t123 * t122;
t205 = t84 * t124;
t8 = t152 + t205 / 0.2e1 + (t146 + t186 / 0.2e1) * t118;
t207 = t8 * qJD(1);
t206 = t84 * t123;
t204 = qJD(1) * t80;
t203 = qJD(5) * t80;
t145 = t57 * t215;
t187 = t123 * t121;
t86 = -t120 * t195 + t187;
t151 = t86 * t217;
t184 = t124 * t122;
t10 = t151 + t206 / 0.2e1 + (t145 - t184 / 0.2e1) * t118;
t202 = t10 * qJD(1);
t113 = t118 ^ 2;
t114 = t119 ^ 2;
t201 = t113 * t114;
t200 = t113 * t119;
t199 = t117 * t118;
t12 = (t118 * t184 - t206) * t80 + t79 * t55;
t193 = t12 * qJD(1);
t147 = t118 * t186;
t13 = -(t147 + t205) * t80 + t79 * t57;
t183 = t13 * qJD(1);
t14 = t55 ^ 2 - t57 ^ 2;
t182 = t14 * qJD(1);
t15 = -t123 * t78 - t82 * t55;
t181 = t15 * qJD(1);
t126 = (t148 / 0.2e1 + t185 / 0.2e1) * t119;
t144 = t80 * t216;
t129 = t124 * t144 + t57 * t214;
t16 = t126 - t129;
t180 = t16 * qJD(1);
t139 = t195 / 0.2e1;
t125 = (t120 * t139 - t187 / 0.2e1) * t119;
t128 = t123 * t144 + t55 * t214;
t17 = t125 + t128;
t179 = t17 * qJD(1);
t20 = (-t55 * t199 + t85 * t80) * t119;
t178 = t20 * qJD(1);
t21 = (-t57 * t199 - t86 * t80) * t119;
t177 = t21 * qJD(1);
t115 = t121 ^ 2;
t105 = -t115 * t119 / 0.2e1;
t127 = t105 + (-t120 ^ 2 / 0.2e1 - t117 ^ 2 / 0.2e1) * t200;
t130 = t79 * t120 / 0.2e1 + t84 * t216;
t23 = t127 + t130;
t176 = t23 * qJD(1);
t24 = t79 * t82 - t84 * t80;
t175 = t24 * qJD(1);
t25 = t124 * t78 + t82 * t57;
t174 = t25 * qJD(1);
t27 = t114 * qJ(2) + (-t118 * t65 + t121 * t66) * t122;
t173 = t27 * qJD(1);
t31 = (t118 * t66 + t121 * t65) * t119;
t172 = t31 * qJD(1);
t131 = t80 * t214 + t82 * t215;
t35 = (-t122 / 0.2e1 + t131) * t118;
t171 = t35 * qJD(1);
t36 = -t80 * t149 + t82 * t150;
t170 = t36 * qJD(1);
t132 = t82 * t214 + t144;
t38 = -t191 / 0.2e1 + t132;
t169 = t38 * qJD(1);
t39 = (-t119 * t79 + t122 * t80) * t118;
t168 = t39 * qJD(1);
t40 = (-t119 * t84 + t122 * t82) * t118;
t167 = t40 * qJD(1);
t42 = t123 * t80;
t166 = t42 * qJD(1);
t46 = t82 ^ 2 + t78;
t165 = t46 * qJD(1);
t52 = t117 * t201 + t80 * t191;
t164 = t52 * qJD(1);
t53 = t120 * t201 + t82 * t191;
t163 = t53 * qJD(1);
t76 = (0.1e1 / 0.2e1 + t113 / 0.2e1 + t115 / 0.2e1) * t119;
t162 = t76 * qJD(1);
t88 = (t113 + t115) * t114;
t161 = t88 * qJD(1);
t102 = t122 ^ 2 + t114;
t89 = t102 * t118;
t160 = t89 * qJD(1);
t90 = t102 * t121;
t159 = t90 * qJD(1);
t158 = qJD(1) * t119;
t157 = qJD(5) * t123;
t156 = qJD(5) * t124;
t100 = t102 * qJ(2);
t155 = t100 * qJD(1);
t154 = t102 * qJD(1);
t153 = qJD(1) * t213;
t143 = t118 * t158;
t142 = t122 * t158;
t141 = qJD(3) * t119 * t122;
t140 = qJD(4) * t198;
t138 = -qJD(5) - t204;
t137 = t80 * t143;
t136 = t82 * t143;
t135 = t118 * t142;
t134 = t121 * t142;
t77 = -t200 / 0.2e1 + t105 + t119 / 0.2e1;
t37 = t191 / 0.2e1 + t132;
t34 = t197 / 0.2e1 + t131 * t118;
t22 = t127 - t130;
t19 = t126 + t129;
t18 = t125 - t128;
t11 = t151 + t118 * t145 - t206 / 0.2e1 + t122 * t139;
t9 = t152 + t118 * t146 - t205 / 0.2e1 - t147 / 0.2e1;
t26 = [0, 0, 0, 0, 0, t102 * qJD(2), t100 * qJD(2), t89 * qJD(2) + t121 * t141, t90 * qJD(2) - t118 * t141, t88 * qJD(3), t27 * qJD(2) - t31 * qJD(3), t39 * qJD(2) + t52 * qJD(3) - t82 * t140, t40 * qJD(2) + t53 * qJD(3) + t80 * t140, t24 * qJD(2) - t36 * qJD(3) + t46 * qJD(4), t3 * qJD(2) - t4 * qJD(3) + t5 * qJD(4), -qJD(5) * t213, t14 * qJD(5), -t55 * t203, -t57 * t203, 0, t12 * qJD(2) + t20 * qJD(3) - t15 * qJD(4) + t2 * qJD(5), t13 * qJD(2) + t21 * qJD(3) + t25 * qJD(4) + t1 * qJD(5); 0, 0, 0, 0, 0, t154, t155, t160, t159, 0, t77 * qJD(3) + t173, t168, t167, t175, t210 + t22 * qJD(3) + t34 * qJD(4) + (t117 * t79 + t120 * t84 - t190) * qJD(2) * t118, 0, 0, 0, 0, 0, t11 * qJD(5) + t193, t9 * qJD(5) + t183; 0, 0, 0, 0, 0, 0, 0, t134, -t135, t161, t77 * qJD(2) - t172, t164, t163, -t170, t22 * qJD(2) + t37 * qJD(4) - t209, 0, 0, 0, 0, 0, t19 * qJD(5) + t178, t18 * qJD(5) + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t137, t165, t34 * qJD(2) + t37 * qJD(3) + t208, 0, 0, 0, 0, 0, -t181, t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t182, t138 * t55, t138 * t57, 0, t11 * qJD(2) + t19 * qJD(3) + t7 * qJD(5) + t211, t9 * qJD(2) + t18 * qJD(3) + t6 * qJD(5) + t212; 0, 0, 0, 0, 0, -t154, -t155, -t160, -t159, 0, -t76 * qJD(3) - t173, -t168, -t167, -t175, t23 * qJD(3) + t35 * qJD(4) - t210, 0, 0, 0, 0, 0, t10 * qJD(5) - t193, t8 * qJD(5) - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, 0, 0, 0, t176, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * qJD(5) + t202, t85 * qJD(5) + t207; 0, 0, 0, 0, 0, 0, 0, -t134, t135, -t161, t76 * qJD(2) + t172, -t164, -t163, t170, -t23 * qJD(2) + t38 * qJD(4) + t209, 0, 0, 0, 0, 0, -t16 * qJD(5) - t178, -t17 * qJD(5) - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, 0, 0, 0, -t176, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 * t156 - t180, t117 * t157 - t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t137, -t165, -t35 * qJD(2) - t38 * qJD(3) - t208, 0, 0, 0, 0, 0, -t42 * qJD(5) + t181, -t156 * t80 - t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157 - t166, t138 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t182, t55 * t204, t57 * t204, 0, -t10 * qJD(2) + t16 * qJD(3) + t42 * qJD(4) - t211, t124 * t80 * qJD(4) - t8 * qJD(2) + t17 * qJD(3) - t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, -t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t124 * t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t26;
