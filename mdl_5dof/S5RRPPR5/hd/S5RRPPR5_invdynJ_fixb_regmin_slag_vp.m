% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:30
% EndTime: 2021-01-15 19:36:37
% DurationCPUTime: 1.76s
% Computational Cost: add. (1829->293), mult. (4210->364), div. (0->0), fcn. (3007->10), ass. (0->156)
t122 = qJD(2) - qJD(5);
t132 = sin(qJ(5));
t135 = cos(qJ(5));
t129 = sin(pkin(8));
t130 = cos(pkin(8));
t136 = cos(qJ(2));
t192 = t130 * t136;
t177 = qJD(1) * t192;
t133 = sin(qJ(2));
t186 = qJD(1) * t133;
t78 = t129 * t186 - t177;
t91 = t129 * t136 + t130 * t133;
t81 = t91 * qJD(1);
t220 = t132 * t78 + t135 * t81;
t196 = t220 * t122;
t181 = t136 * qJDD(1);
t182 = t133 * qJDD(1);
t165 = t129 * t182 - t130 * t181;
t80 = t91 * qJD(2);
t51 = qJD(1) * t80 + t165;
t183 = qJD(1) * qJD(2);
t175 = t133 * t183;
t147 = t91 * qJDD(1) - t129 * t175;
t174 = t136 * t183;
t52 = t130 * t174 + t147;
t5 = qJD(5) * t220 + t132 * t52 - t135 * t51;
t230 = t5 + t196;
t205 = t136 * pkin(2);
t117 = pkin(1) + t205;
t71 = pkin(2) * t175 - t117 * qJDD(1) + qJDD(3);
t229 = -t52 * qJ(4) + t71;
t98 = -t117 * qJD(1) + qJD(3);
t228 = -t81 * qJ(4) + t98;
t134 = sin(qJ(1));
t137 = cos(qJ(1));
t221 = g(1) * t134 - g(2) * t137;
t160 = t132 * t81 - t135 * t78;
t170 = t160 * qJD(5) - t132 * t51 - t135 * t52;
t195 = t160 * t122;
t227 = t170 + t195;
t167 = g(1) * t137 + g(2) * t134;
t226 = -t160 ^ 2 + t220 ^ 2;
t214 = -pkin(3) - pkin(4);
t14 = t214 * t78 - t228;
t131 = -qJ(3) - pkin(6);
t171 = qJD(2) * t131;
t149 = -t133 * qJD(3) + t136 * t171;
t176 = t131 * t133;
t47 = qJDD(2) * pkin(2) + t149 * qJD(1) + qJDD(1) * t176;
t100 = t131 * t136;
t74 = t136 * qJD(3) + t133 * t171;
t55 = t74 * qJD(1) - qJDD(1) * t100;
t203 = t129 * t55 - t130 * t47;
t178 = -qJDD(4) - t203;
t2 = -t52 * pkin(7) + t214 * qJDD(2) - t178;
t10 = t129 * t47 + t130 * t55;
t124 = qJDD(2) * qJ(4);
t125 = qJD(2) * qJD(4);
t7 = t124 + t125 + t10;
t3 = t51 * pkin(7) + t7;
t123 = qJ(2) + pkin(8);
t118 = sin(t123);
t119 = cos(t123);
t72 = -t118 * t132 - t119 * t135;
t73 = t118 * t135 - t119 * t132;
t225 = g(3) * t72 + t132 * t3 - t135 * t2 + t14 * t220 + t167 * t73;
t76 = t81 ^ 2;
t223 = -t78 ^ 2 - t76;
t222 = t220 * t160;
t96 = qJD(1) * t100;
t202 = t129 * t96;
t95 = qJD(1) * t176;
t57 = t130 * t95 + t202;
t188 = qJD(4) - t57;
t219 = qJD(5) + t122;
t61 = -t130 * t100 + t129 * t176;
t217 = -t61 * qJDD(2) - t118 * t221;
t216 = g(3) * t118 + t57 * qJD(2) + t167 * t119 - t10;
t215 = g(3) * t73 - t132 * t2 - t135 * t3 + t14 * t160 - t167 * t72;
t213 = t78 * pkin(7);
t212 = t81 * pkin(7);
t211 = pkin(2) * t133;
t207 = g(3) * t119;
t206 = g(3) * t136;
t27 = t78 * pkin(3) + t228;
t204 = t27 * t81;
t30 = t129 * t149 + t130 * t74;
t201 = t130 * t96;
t89 = qJD(2) * pkin(2) + t95;
t50 = t129 * t89 - t201;
t191 = t134 * t131;
t190 = t81 * qJD(4);
t189 = -t212 + t188;
t127 = t133 ^ 2;
t187 = -t136 ^ 2 + t127;
t185 = qJD(2) * t133;
t41 = qJD(2) * qJ(4) + t50;
t179 = pkin(2) * t185;
t115 = -t130 * pkin(2) - pkin(3);
t29 = t129 * t74 - t130 * t149;
t56 = t129 * t95 - t201;
t49 = t130 * t89 + t202;
t60 = -t129 * t100 - t130 * t176;
t169 = t122 ^ 2;
t168 = qJD(4) - t49;
t16 = t214 * qJD(2) + t168 - t212;
t22 = t41 + t213;
t164 = t132 * t22 - t135 * t16;
t163 = -t132 * t16 - t135 * t22;
t31 = -t91 * pkin(7) + t60;
t90 = t129 * t133 - t192;
t32 = t90 * pkin(7) + t61;
t162 = -t132 * t32 + t135 * t31;
t161 = t132 * t31 + t135 * t32;
t53 = t132 * t91 - t135 * t90;
t54 = t132 * t90 + t135 * t91;
t111 = -pkin(4) + t115;
t113 = t129 * pkin(2) + qJ(4);
t158 = t135 * t111 - t132 * t113;
t157 = t132 * t111 + t135 * t113;
t156 = t91 * qJ(4) + t117;
t8 = -qJDD(2) * pkin(3) - t178;
t155 = -pkin(2) * t186 - t78 * qJ(4);
t154 = -0.2e1 * pkin(1) * t183 - pkin(6) * qJDD(2);
t150 = -t60 * qJDD(2) + t221 * t119;
t83 = qJD(2) * t192 - t129 * t185;
t148 = t83 * qJ(4) + t91 * qJD(4) - t179;
t146 = t56 * qJD(2) + t167 * t118 - t203 - t207;
t138 = qJD(2) ^ 2;
t145 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t138 + t221;
t139 = qJD(1) ^ 2;
t144 = pkin(1) * t139 - pkin(6) * qJDD(1) + t167;
t142 = t29 * t81 - t30 * t78 - t61 * t51 + t60 * t52 - t167;
t141 = t51 * pkin(3) + t229;
t140 = 0.2e1 * t81 * qJD(2) + t165;
t121 = qJDD(2) - qJDD(5);
t116 = t131 * t137;
t99 = -t129 * pkin(3) + qJ(4) * t130;
t97 = pkin(3) * t130 + qJ(4) * t129 + pkin(2);
t58 = t99 * t133 + t97 * t136 + pkin(1);
t48 = t90 * pkin(3) - t156;
t33 = -qJD(2) * pkin(3) + t168;
t28 = t81 * pkin(3) - t155;
t25 = t56 + t213;
t24 = t214 * t90 + t156;
t23 = (-t78 + t177) * qJD(2) + t147;
t21 = t80 * pkin(3) - t148;
t19 = t80 * pkin(7) + t30;
t18 = -t83 * pkin(7) + t29;
t17 = t214 * t81 + t155;
t13 = t214 * t80 + t148;
t12 = t54 * qJD(5) + t132 * t83 - t135 * t80;
t11 = -t53 * qJD(5) + t132 * t80 + t135 * t83;
t6 = t141 - t190;
t1 = t214 * t51 + t190 - t229;
t4 = [qJDD(1), t221, t167, t127 * qJDD(1) + 0.2e1 * t133 * t174, 0.2e1 * t133 * t181 - 0.2e1 * t187 * t183, qJDD(2) * t133 + t138 * t136, qJDD(2) * t136 - t138 * t133, 0, t154 * t133 + t145 * t136, -t145 * t133 + t154 * t136, -t117 * t51 + t71 * t90 + t98 * t80 + (t78 * t211 - t29) * qJD(2) + t150, -t117 * t52 + t71 * t91 + t98 * t83 + (t81 * t211 - t30) * qJD(2) + t217, -t10 * t90 + t203 * t91 - t49 * t83 - t50 * t80 + t142, t10 * t61 + t50 * t30 + t203 * t60 - t49 * t29 - t71 * t117 + t98 * t179 - g(1) * (-t134 * t117 - t116) - g(2) * (t137 * t117 - t191), -t29 * qJD(2) + t21 * t78 + t27 * t80 + t48 * t51 + t6 * t90 + t150, t33 * t83 - t41 * t80 - t7 * t90 + t8 * t91 + t142, t30 * qJD(2) - t21 * t81 - t27 * t83 - t48 * t52 - t6 * t91 - t217, t7 * t61 + t41 * t30 + t6 * t48 + t27 * t21 + t8 * t60 + t33 * t29 - g(1) * (-t58 * t134 - t116) - g(2) * (t58 * t137 - t191), t11 * t220 - t170 * t54, -t11 * t160 - t12 * t220 + t170 * t53 - t54 * t5, -t11 * t122 - t54 * t121, t12 * t122 + t53 * t121, 0, t13 * t160 + t24 * t5 + t1 * t53 + t14 * t12 - (-qJD(5) * t161 - t132 * t19 + t135 * t18) * t122 - t162 * t121 - t221 * t72, t13 * t220 - t24 * t170 + t1 * t54 + t14 * t11 + (qJD(5) * t162 + t132 * t18 + t135 * t19) * t122 + t161 * t121 + t221 * t73; 0, 0, 0, -t133 * t139 * t136, t187 * t139, t182, t181, qJDD(2), t144 * t133 - t206, g(3) * t133 + t144 * t136, -t98 * t81 + (qJDD(2) * t130 - t78 * t186) * pkin(2) + t146, t98 * t78 + (-qJDD(2) * t129 - t81 * t186) * pkin(2) + t216, (t50 - t56) * t81 + (-t49 + t57) * t78 + (-t129 * t51 - t130 * t52) * pkin(2), t49 * t56 - t50 * t57 + (-t206 + t10 * t129 - t130 * t203 + (-qJD(1) * t98 + t167) * t133) * pkin(2), -t204 - t28 * t78 - qJDD(4) + (pkin(3) - t115) * qJDD(2) + t146, -t113 * t51 + t115 * t52 + (t41 - t56) * t81 + (t33 - t188) * t78, t113 * qJDD(2) - t27 * t78 + t28 * t81 + t124 + 0.2e1 * t125 - t216, t7 * t113 + t8 * t115 - t27 * t28 - t33 * t56 - g(3) * (t119 * pkin(3) + t118 * qJ(4) + t205) - t167 * (-t97 * t133 + t99 * t136) + t188 * t41, -t222, -t226, t227, t230, t121, -t158 * t121 - t17 * t160 + (t132 * t189 + t135 * t25) * t122 + (t122 * t157 - t163) * qJD(5) + t225, t157 * t121 - t17 * t220 + (-t132 * t25 + t135 * t189) * t122 + (t122 * t158 - t164) * qJD(5) - t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t23, t223, t49 * t81 + t50 * t78 - t221 + t71, t140, t223, -t23, t41 * t78 + (-qJD(4) - t33) * t81 + t141 - t221, 0, 0, 0, 0, 0, -t5 + t196, t170 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t78 - qJDD(2), (t78 + t177) * qJD(2) + t147, -t76 - t138, -t41 * qJD(2) - t167 * t91 + t204 + t207 + t8, 0, 0, 0, 0, 0, -t135 * t121 - t132 * t169 - t160 * t81, t132 * t121 - t135 * t169 - t220 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t226, -t227, -t230, -t121, t163 * t219 - t225, t164 * t219 + t215;];
tau_reg = t4;
