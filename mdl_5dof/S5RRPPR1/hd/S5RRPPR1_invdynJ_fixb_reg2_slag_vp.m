% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:48
% DurationCPUTime: 1.52s
% Computational Cost: add. (2420->242), mult. (3882->295), div. (0->0), fcn. (2553->16), ass. (0->164)
t143 = qJ(1) + qJ(2);
t132 = pkin(8) + t143;
t117 = cos(t132);
t107 = g(2) * t117;
t138 = qJDD(1) + qJDD(2);
t145 = sin(pkin(8));
t147 = cos(pkin(8));
t153 = cos(qJ(2));
t205 = pkin(1) * qJD(2);
t187 = qJD(1) * t205;
t150 = sin(qJ(2));
t192 = qJDD(1) * t150;
t225 = pkin(1) * t192 + t153 * t187;
t212 = t153 * pkin(1);
t125 = qJDD(1) * t212;
t66 = pkin(2) * t138 - t150 * t187 + t125;
t190 = t145 * t225 - t147 * t66;
t176 = qJDD(4) + t190;
t33 = -t138 * pkin(3) + t176;
t226 = t33 + t107;
t206 = pkin(1) * qJD(1);
t189 = t150 * t206;
t100 = t145 * t189;
t188 = t153 * t206;
t71 = t147 * t188 - t100;
t196 = qJD(4) - t71;
t116 = sin(t132);
t108 = g(1) * t116;
t222 = t108 - t107;
t146 = cos(pkin(9));
t152 = cos(qJ(5));
t197 = t152 * t146;
t144 = sin(pkin(9));
t149 = sin(qJ(5));
t200 = t149 * t144;
t80 = -t197 + t200;
t81 = t152 * t144 + t149 * t146;
t142 = qJD(1) + qJD(2);
t63 = t81 * t142;
t139 = t144 ^ 2;
t140 = t146 ^ 2;
t193 = t139 + t140;
t224 = t142 * t193;
t223 = g(1) * t117 + g(2) * t116;
t133 = sin(t143);
t134 = cos(t143);
t221 = g(1) * t133 - g(2) * t134;
t220 = t63 ^ 2;
t219 = pkin(1) * t150;
t218 = pkin(2) * t133;
t216 = t145 * pkin(2);
t215 = t146 * pkin(4);
t135 = t146 * pkin(7);
t214 = t147 * pkin(2);
t151 = sin(qJ(1));
t213 = t151 * pkin(1);
t186 = t142 * t200;
t61 = -t142 * t197 + t186;
t211 = t63 * t61;
t114 = qJ(4) + t216;
t77 = (-pkin(7) - t114) * t144;
t78 = t114 * t146 + t135;
t45 = -t149 * t78 + t152 * t77;
t210 = qJD(5) * t45 - t196 * t80;
t46 = t149 * t77 + t152 * t78;
t209 = -qJD(5) * t46 - t196 * t81;
t208 = t226 * t144;
t173 = t80 * t138;
t76 = t81 * qJD(5);
t39 = t142 * t76 + t173;
t182 = qJD(5) * t197;
t183 = qJD(5) * t200;
t75 = -t182 + t183;
t207 = -t39 * t81 + t61 * t75;
t201 = t147 * t150;
t157 = pkin(1) * (t145 * t153 + t201);
t69 = qJD(1) * t157;
t204 = t69 * t142;
t70 = qJD(2) * t157;
t203 = t70 * t142;
t124 = pkin(2) + t212;
t74 = pkin(1) * t201 + t124 * t145;
t41 = t145 * t66 + t147 * t225;
t31 = qJ(4) * t138 + qJD(4) * t142 + t41;
t21 = qJDD(3) * t144 + t146 * t31;
t85 = pkin(2) * t142 + t188;
t56 = t145 * t85 + t147 * t189;
t52 = qJ(4) * t142 + t56;
t43 = qJD(3) * t144 + t146 * t52;
t202 = t146 * t138;
t194 = g(1) * t134 + g(2) * t133;
t191 = t138 * t81 + t142 * t182;
t112 = t145 * t219;
t123 = pkin(2) * t134;
t184 = pkin(3) * t117 + qJ(4) * t116 + t123;
t118 = pkin(3) + t215;
t127 = t146 * qJDD(3);
t15 = t127 + (-pkin(7) * t138 - t31) * t144;
t16 = pkin(7) * t202 + t21;
t180 = -t149 * t16 + t15 * t152;
t55 = t147 * t85 - t100;
t179 = t193 * t138;
t73 = t124 * t147 - t112;
t178 = qJD(1) * (-qJD(2) + t142);
t177 = qJD(2) * (-qJD(1) - t142);
t68 = -pkin(3) - t73;
t174 = t125 + t221;
t148 = -pkin(7) - qJ(4);
t172 = -t116 * t148 + t117 * t118 + t123;
t171 = -t41 + t223;
t170 = -t190 + t222;
t154 = cos(qJ(1));
t168 = g(1) * t151 - g(2) * t154;
t38 = t142 * t183 - t191;
t167 = -t38 * t80 + t63 * t76;
t166 = qJD(4) - t55;
t129 = t146 * qJD(3);
t36 = t129 + (-pkin(7) * t142 - t52) * t144;
t37 = t135 * t142 + t43;
t10 = -t149 * t37 + t152 * t36;
t11 = t149 * t36 + t152 * t37;
t163 = t149 * t15 + t152 * t16;
t3 = t10 * qJD(5) + t163;
t4 = -t11 * qJD(5) + t180;
t165 = t10 * t75 - t11 * t76 - t3 * t80 - t4 * t81 - t223;
t164 = t138 * t68 + t203;
t67 = qJ(4) + t74;
t53 = (-pkin(7) - t67) * t144;
t54 = t146 * t67 + t135;
t22 = -t149 * t54 + t152 * t53;
t23 = t149 * t53 + t152 * t54;
t119 = -pkin(3) - t214;
t162 = t119 * t138 - t204;
t20 = -t144 * t31 + t127;
t161 = -t20 * t144 + t146 * t21 - t223;
t72 = t147 * t153 * t205 - qJD(2) * t112;
t160 = -pkin(3) * t116 + qJ(4) * t117 - t218;
t141 = pkin(9) + qJ(5);
t130 = sin(t141);
t24 = -t118 * t138 + t176;
t44 = -t118 * t142 + t166;
t159 = -t130 * t222 + t24 * t81 - t44 * t75;
t131 = cos(t141);
t158 = t131 * t222 + t24 * t80 + t44 * t76;
t156 = -t116 * t118 - t117 * t148 - t218;
t136 = t154 * pkin(1);
t111 = t140 * t138;
t110 = t139 * t138;
t90 = t146 * t108;
t88 = -t118 - t214;
t86 = 0.2e1 * t144 * t202;
t65 = qJD(4) + t72;
t60 = t61 ^ 2;
t57 = t68 - t215;
t51 = -pkin(3) * t142 + t166;
t49 = -qJD(5) * t76 - qJDD(5) * t80;
t48 = -qJD(5) * t75 + qJDD(5) * t81;
t42 = -t144 * t52 + t129;
t13 = t39 * t80 + t61 * t76;
t12 = -t38 * t81 - t63 * t75;
t9 = -qJD(5) * t23 - t65 * t81;
t8 = qJD(5) * t22 - t65 * t80;
t5 = -t167 + t207;
t1 = [0, 0, 0, 0, 0, qJDD(1), t168, g(1) * t154 + g(2) * t151, 0, 0, 0, 0, 0, 0, 0, t138, (t138 * t153 + t150 * t177) * pkin(1) + t174, ((-qJDD(1) - t138) * t150 + t153 * t177) * pkin(1) + t194, 0, (t168 + (t150 ^ 2 + t153 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t138, t138 * t73 + t170 - t203, -t138 * t74 - t142 * t72 + t171, 0, t41 * t74 + t56 * t72 - t190 * t73 - t55 * t70 - g(1) * (-t213 - t218) - g(2) * (t123 + t136), t110, t86, 0, t111, 0, 0, t90 + (-t164 - t226) * t146, (t164 - t108) * t144 + t208, t179 * t67 + t224 * t65 + t161, t33 * t68 + t51 * t70 - g(1) * (t160 - t213) - g(2) * (t136 + t184) + (t21 * t67 + t43 * t65) * t146 + (-t20 * t67 - t42 * t65) * t144, t12, t5, t48, t13, t49, 0, qJD(5) * t9 + qJDD(5) * t22 + t39 * t57 + t61 * t70 + t158, -qJD(5) * t8 - qJDD(5) * t23 - t38 * t57 + t63 * t70 + t159, t22 * t38 - t23 * t39 - t61 * t8 - t63 * t9 + t165, t3 * t23 + t11 * t8 + t4 * t22 + t10 * t9 + t24 * t57 + t44 * t70 - g(1) * (t156 - t213) - g(2) * (t136 + t172); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, t178 * t219 + t174, (t153 * t178 - t192) * pkin(1) + t194, 0, 0, 0, 0, 0, 0, 0, t138, t138 * t214 + t170 + t204, -t138 * t216 + t142 * t71 + t171, 0, t55 * t69 - t56 * t71 + (t145 * t41 - t147 * t190 + t221) * pkin(2), t110, t86, 0, t111, 0, 0, t90 + (-t162 - t226) * t146, (t162 - t108) * t144 + t208, t114 * t179 + t196 * t224 + t161, t33 * t119 - t51 * t69 - g(1) * t160 - g(2) * t184 + (t21 * t114 + t196 * t43) * t146 + (-t114 * t20 - t196 * t42) * t144, t12, t5, t48, t13, t49, 0, qJD(5) * t209 + t45 * qJDD(5) + t88 * t39 - t69 * t61 + t158, -qJD(5) * t210 - t46 * qJDD(5) - t88 * t38 - t69 * t63 + t159, -t209 * t63 - t210 * t61 + t45 * t38 - t46 * t39 + t165, -g(1) * t156 - g(2) * t172 + t10 * t209 + t11 * t210 + t24 * t88 + t3 * t46 + t4 * t45 - t44 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t21 + t146 * t20 - g(3), 0, 0, 0, 0, 0, 0, t49, -t48, t167 + t207, -t10 * t76 - t11 * t75 + t3 * t81 - t4 * t80 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t144 * t138, -t193 * t142 ^ 2, (t144 * t42 - t146 * t43) * t142 + t33 - t222, 0, 0, 0, 0, 0, 0, 0.2e1 * t63 * qJD(5) + t173, (-t61 - t186) * qJD(5) + t191, -t60 - t220, t10 * t63 + t11 * t61 - t222 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, -t60 + t220, (t61 - t186) * qJD(5) + t191, -t211, -t173, qJDD(5), -g(3) * t131 + t130 * t223 - t44 * t63 + t180, g(3) * t130 + t131 * t223 + t44 * t61 - t163, 0, 0;];
tau_reg = t1;
