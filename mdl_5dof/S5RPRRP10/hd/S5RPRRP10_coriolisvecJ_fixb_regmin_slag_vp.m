% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP10
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
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:48
% EndTime: 2021-01-15 19:15:01
% DurationCPUTime: 1.97s
% Computational Cost: add. (2508->249), mult. (6595->331), div. (0->0), fcn. (4830->6), ass. (0->123)
t175 = 2 * qJD(3);
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t172 = -t100 * t97 + t102 * t98;
t174 = t172 * qJD(1);
t149 = pkin(6) + qJ(2);
t83 = t149 * t97;
t80 = qJD(1) * t83;
t84 = t149 * t98;
t81 = qJD(1) * t84;
t52 = -t100 * t80 + t102 * t81;
t173 = qJD(3) * t52;
t101 = cos(qJ(4));
t130 = qJD(4) * t101;
t99 = sin(qJ(4));
t132 = qJD(3) * t99;
t79 = t100 * t98 + t102 * t97;
t74 = t79 * qJD(1);
t26 = t74 * t130 + (qJD(4) + t174) * t132;
t60 = t101 * t74 + t132;
t155 = t60 * t99;
t91 = -pkin(2) * t98 - pkin(1);
t82 = t91 * qJD(1) + qJD(2);
t29 = -pkin(3) * t174 - t74 * pkin(7) + t82;
t47 = qJD(3) * pkin(7) + t52;
t15 = t101 * t47 + t29 * t99;
t108 = t172 * qJD(2);
t168 = -t100 * t81 - t102 * t80;
t20 = qJD(1) * t108 + qJD(3) * t168;
t75 = t172 * qJD(3);
t106 = qJD(1) * t75;
t76 = t79 * qJD(3);
t64 = qJD(1) * t76;
t35 = t64 * pkin(3) - pkin(7) * t106;
t34 = t101 * t35;
t105 = -qJD(4) * t15 - t99 * t20 + t34;
t129 = t101 * qJD(3);
t131 = qJD(4) * t99;
t25 = -qJD(4) * t129 - t101 * t106 + t74 * t131;
t104 = t25 * qJ(5) + t105;
t163 = t64 * pkin(4);
t1 = -t60 * qJD(5) + t104 + t163;
t67 = qJD(4) - t174;
t58 = t74 * t99 - t129;
t9 = -qJ(5) * t58 + t15;
t162 = t67 * t9;
t170 = t1 + t162;
t169 = t67 * t155;
t55 = t100 * t84 + t102 * t83;
t166 = t60 ^ 2;
t14 = t101 * t29 - t47 * t99;
t8 = -qJ(5) * t60 + t14;
t5 = pkin(4) * t67 + t8;
t165 = t5 - t8;
t164 = pkin(4) * t58;
t161 = t25 * t99;
t160 = t58 * t67;
t159 = t58 * t174;
t158 = t58 * t74;
t157 = t60 * t67;
t156 = t60 * t74;
t153 = t174 * t99;
t152 = t79 * t64;
t151 = t79 * t99;
t150 = t99 * t64;
t148 = -qJ(5) - pkin(7);
t121 = qJD(4) * t148;
t133 = qJ(5) * t101;
t48 = pkin(3) * t74 - pkin(7) * t174;
t40 = t101 * t48;
t147 = pkin(4) * t74 - t101 * t121 - t174 * t133 + t40 + (qJD(5) - t168) * t99;
t141 = qJ(5) * t99;
t144 = t101 * t168 + t99 * t48;
t146 = -t101 * qJD(5) - t121 * t99 - t141 * t174 + t144;
t145 = -t58 * t130 - t99 * t26;
t50 = -pkin(3) * t172 - pkin(7) * t79 + t91;
t56 = -t100 * t83 + t102 * t84;
t53 = t101 * t56;
t143 = t99 * t50 + t53;
t142 = t97 ^ 2 + t98 ^ 2;
t122 = -qJD(5) - t164;
t46 = -qJD(3) * pkin(3) - t168;
t18 = -t122 + t46;
t138 = t101 * t18;
t128 = qJD(1) * qJD(2);
t30 = -t55 * qJD(3) + t108;
t49 = pkin(3) * t76 - pkin(7) * t75;
t127 = t101 * t30 + t50 * t130 + t99 * t49;
t126 = t79 * t131;
t124 = t79 * t130;
t123 = t142 * qJD(1) ^ 2;
t120 = -t101 * t20 - t29 * t130 + t47 * t131 - t99 * t35;
t21 = t79 * t128 + t173;
t119 = t101 * t67;
t111 = qJ(5) * t26 + t120;
t2 = -qJD(5) * t58 - t111;
t118 = -t67 * t5 + t2;
t117 = t21 * t79 + t46 * t75;
t115 = -qJ(5) * t75 - qJD(5) * t79;
t114 = 0.2e1 * t142 * t128;
t113 = t101 * t64 + (-t131 + t153) * t67;
t11 = pkin(4) * t26 + t21;
t112 = t75 * t99 + t124;
t110 = -pkin(7) * t64 + t67 * t46;
t31 = qJD(2) * t79 + qJD(3) * t56;
t93 = -pkin(4) * t101 - pkin(3);
t86 = t148 * t101;
t85 = t148 * t99;
t57 = t58 ^ 2;
t43 = t101 * t50;
t41 = t101 * t49;
t32 = pkin(4) * t151 + t55;
t22 = pkin(4) * t153 + t52;
t17 = pkin(4) * t112 + t31;
t16 = -t79 * t141 + t143;
t12 = -pkin(4) * t172 - t79 * t133 - t56 * t99 + t43;
t7 = -t67 ^ 2 * t101 - t150 - t156;
t6 = t113 - t158;
t4 = -qJ(5) * t124 + (-qJD(4) * t56 + t115) * t99 + t127;
t3 = t76 * pkin(4) - t99 * t30 + t41 + t115 * t101 + (-t53 + (qJ(5) * t79 - t50) * t99) * qJD(4);
t10 = [0, 0, 0, 0, t114, qJ(2) * t114, t106 * t79 + t74 * t75, t106 * t172 + t174 * t75 - t74 * t76 - t152, t75 * qJD(3), -t76 * qJD(3), 0, -qJD(3) * t31 + t64 * t91 + t76 * t82, t82 * t75 + (t174 * t91 - t30) * qJD(3), -t60 * t126 + (-t25 * t79 + t60 * t75) * t101, (-t101 * t58 - t155) * t75 + (-t101 * t26 + t161 + (-t101 * t60 + t58 * t99) * qJD(4)) * t79, -t67 * t126 + t25 * t172 + t60 * t76 + (t67 * t75 + t152) * t101, -t112 * t67 - t79 * t150 + t172 * t26 - t58 * t76, -t172 * t64 + t67 * t76, (-t56 * t130 + t41) * t67 + t43 * t64 - (-t47 * t130 + t34) * t172 + t14 * t76 + t31 * t58 + t55 * t26 + t46 * t124 + ((-qJD(4) * t50 - t30) * t67 - t56 * t64 - (-qJD(4) * t29 - t20) * t172 + t117) * t99, -(-t56 * t131 + t127) * t67 - t143 * t64 - t120 * t172 - t15 * t76 + t31 * t60 - t55 * t25 - t46 * t126 + t117 * t101, -t1 * t172 + t11 * t151 + t112 * t18 + t12 * t64 + t17 * t58 + t26 * t32 + t3 * t67 + t5 * t76, -t18 * t126 - t16 * t64 + t17 * t60 + t2 * t172 - t25 * t32 - t4 * t67 - t76 * t9 + (t11 * t79 + t18 * t75) * t101, t12 * t25 - t16 * t26 - t3 * t60 - t4 * t58 + (-t101 * t5 - t9 * t99) * t75 + (-t1 * t101 - t2 * t99 + (-t101 * t9 + t5 * t99) * qJD(4)) * t79, t1 * t12 + t11 * t32 + t16 * t2 + t17 * t18 + t3 * t5 + t4 * t9; 0, 0, 0, 0, -t123, -qJ(2) * t123, 0, 0, 0, 0, 0, t74 * t175, t174 * t175, 0, 0, 0, 0, 0, t6, t7, t6, t7, t169 + (t25 + t159) * t101 + t145, t170 * t101 + t118 * t99 - t18 * t74; 0, 0, 0, 0, 0, 0, -t74 * t174, -t174 ^ 2 + t74 ^ 2, 0, 0, 0, -t74 * t82 + t173 - t21, (-qJD(2) - t82) * t174, t119 * t60 - t161, -t169 + (-t25 + t159) * t101 + t145, t119 * t67 + t150 - t156, t113 + t158, -t67 * t74, -pkin(3) * t26 - t14 * t74 - t40 * t67 - t52 * t58 + (-pkin(7) * qJD(4) * t67 - t21) * t101 + (t168 * t67 + t110) * t99, pkin(3) * t25 + t15 * t74 + t21 * t99 - t52 * t60 + (pkin(7) * t131 + t144) * t67 + t110 * t101, -t101 * t11 - t22 * t58 + t26 * t93 - t5 * t74 + t64 * t85 - t147 * t67 + (-t18 * t174 + (t18 + t164) * qJD(4)) * t99, -t174 * t138 + t11 * t99 - t22 * t60 - t25 * t93 + t64 * t86 + t74 * t9 + t146 * t67 + (pkin(4) * t155 + t138) * qJD(4), t118 * t101 + t146 * t58 + t147 * t60 - t170 * t99 + t25 * t85 + t26 * t86, t1 * t85 + t11 * t93 - t2 * t86 - t146 * t9 - t147 * t5 + (pkin(4) * t131 - t22) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t57 + t166, -t25 + t160, t157 - t26, t64, t15 * t67 - t46 * t60 + t105, t14 * t67 + t46 * t58 + t120, 0.2e1 * t163 + t162 + (t122 - t18) * t60 + t104, -pkin(4) * t166 + t67 * t8 + (qJD(5) + t18) * t58 + t111, t25 * pkin(4) - t165 * t58, t165 * t9 + (-t18 * t60 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 + t157, -t25 - t160, -t57 - t166, t5 * t60 + t58 * t9 + t11;];
tauc_reg = t10;
