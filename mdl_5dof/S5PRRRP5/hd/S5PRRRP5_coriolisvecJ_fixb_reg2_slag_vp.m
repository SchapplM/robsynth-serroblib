% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:24
% DurationCPUTime: 1.29s
% Computational Cost: add. (1801->217), mult. (4653->283), div. (0->0), fcn. (3205->6), ass. (0->138)
t101 = sin(qJ(4));
t102 = sin(qJ(3));
t104 = cos(qJ(3));
t172 = cos(qJ(4));
t77 = t101 * t104 + t102 * t172;
t68 = t77 * qJD(2);
t157 = t68 * qJ(5);
t103 = sin(qJ(2));
t143 = t103 * qJD(1);
t160 = qJD(2) * pkin(6);
t88 = t143 + t160;
t129 = pkin(7) * qJD(2) + t88;
t61 = t129 * t104;
t52 = t101 * t61;
t60 = t129 * t102;
t55 = qJD(3) * pkin(3) - t60;
t26 = t172 * t55 - t52;
t14 = t26 - t157;
t98 = qJD(3) + qJD(4);
t177 = t98 * t77;
t35 = t177 * qJD(2);
t105 = cos(qJ(2));
t134 = t172 * t104;
t154 = t101 * t102;
t115 = t134 - t154;
t113 = t105 * t115;
t132 = qJD(4) * t172;
t144 = qJD(4) * t101;
t174 = -pkin(7) - pkin(6);
t135 = qJD(3) * t174;
t81 = t102 * t135;
t82 = t104 * t135;
t83 = t174 * t102;
t84 = t174 * t104;
t165 = -qJD(1) * t113 + t101 * t82 + t132 * t83 + t144 * t84 + t172 * t81;
t142 = t105 * qJD(1);
t48 = t101 * t83 - t172 * t84;
t164 = -qJD(4) * t48 - t101 * t81 + t77 * t142 + t172 * t82;
t106 = qJD(3) ^ 2;
t107 = qJD(2) ^ 2;
t176 = (t106 + t107) * t103;
t63 = t115 * t103;
t146 = qJD(3) * t102;
t116 = pkin(3) * t146 - t143;
t175 = t68 ^ 2;
t173 = qJ(5) * t177 - qJD(5) * t115 - t165;
t122 = qJD(2) * t134;
t148 = qJD(2) * t102;
t133 = t101 * t148;
t66 = -t122 + t133;
t171 = t66 * t98;
t170 = t68 * t66;
t169 = t68 * t98;
t97 = -pkin(3) * t104 - pkin(2);
t73 = qJD(2) * t97 - t142;
t168 = t73 * t68;
t120 = t98 * t154;
t42 = -qJD(3) * t134 - t104 * t132 + t120;
t167 = t42 * qJ(5) - t77 * qJD(5) + t164;
t13 = pkin(4) * t98 + t14;
t166 = t13 - t14;
t123 = pkin(3) * t132;
t163 = -pkin(3) * t101 * t35 - t123 * t66;
t29 = -t172 * t60 - t52;
t162 = t98 * t122;
t141 = qJD(2) * qJD(3);
t131 = t102 * t141;
t74 = pkin(3) * t131 + qJD(2) * t143;
t161 = qJD(2) * pkin(2);
t34 = qJD(2) * t120 - t162;
t159 = t34 * qJ(5);
t158 = t66 * qJ(5);
t100 = t104 ^ 2;
t99 = t102 ^ 2;
t156 = t100 - t99;
t155 = t100 + t99;
t153 = t106 * t102;
t152 = t106 * t104;
t151 = t107 * t105;
t130 = t66 * pkin(4) + qJD(5);
t41 = t130 + t73;
t150 = qJD(5) + t41;
t147 = qJD(2) * t103;
t145 = qJD(3) * t104;
t140 = t172 * pkin(3);
t139 = pkin(3) * t144;
t137 = pkin(3) * t148;
t54 = t172 * t61;
t136 = t102 * t107 * t104;
t44 = -t88 * t146 + (-pkin(7) * t146 + t104 * t142) * qJD(2);
t45 = -t88 * t145 + (-pkin(7) * t145 - t102 * t142) * qJD(2);
t128 = -t101 * t44 + t172 * t45;
t28 = t101 * t60 - t54;
t47 = t101 * t84 + t172 * t83;
t127 = -t101 * t45 - t132 * t55 + t144 * t61 - t172 * t44;
t126 = pkin(4) * t177 + t116;
t89 = -t142 - t161;
t125 = -t89 - t142;
t124 = -0.2e1 * t105 * t141;
t24 = pkin(4) * t35 + t74;
t121 = t104 * t131;
t119 = qJD(2) * t125;
t118 = t66 * t73 + t127;
t117 = t35 * qJ(5) + t127;
t27 = t101 * t55 + t54;
t114 = qJD(3) * (-t125 - t161);
t112 = -t133 * t98 + t162;
t111 = t150 * t66 + t117;
t8 = -qJD(4) * t27 + t128;
t110 = t8 + t159;
t109 = (-t54 + (-pkin(3) * t98 - t55) * t101) * qJD(4) + t128;
t96 = t140 + pkin(4);
t85 = t98 * t123;
t65 = t66 ^ 2;
t62 = t77 * t103;
t59 = -pkin(4) * t115 + t97;
t50 = pkin(4) * t68 + t137;
t37 = t177 * t98;
t36 = t42 * t98;
t31 = qJ(5) * t115 + t48;
t30 = -qJ(5) * t77 + t47;
t25 = -t65 + t175;
t23 = t169 - t35;
t22 = t112 + t171;
t19 = -t105 * t68 - t63 * t98;
t18 = qJD(2) * t113 - t103 * t177;
t17 = -t157 + t29;
t16 = t28 + t158;
t15 = t27 - t158;
t12 = -t115 * t35 + t177 * t66;
t11 = -t34 * t77 - t42 * t68;
t6 = -t105 * t35 + t147 * t66 + t19 * t98;
t5 = t105 * t34 + t147 * t68 - t18 * t98;
t4 = -t68 * qJD(5) + t110;
t3 = -qJD(5) * t66 - t117;
t2 = -t115 * t34 - t177 * t68 - t35 * t77 + t42 * t66;
t1 = -t18 * t66 - t19 * t68 - t34 * t62 - t35 * t63;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t103, -t151, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t124 - t104 * t176, t102 * t176 + t104 * t124, t155 * t151, (t103 * t89 + (-t143 + (t88 + t143) * t155) * t105) * qJD(2), 0, 0, 0, 0, 0, 0, t6, t5, t1, -t105 * t74 - t127 * t63 + t147 * t73 + t18 * t27 + t19 * t26 - t62 * t8, 0, 0, 0, 0, 0, 0, t6, t5, t1, -t105 * t24 + t13 * t19 + t147 * t41 + t15 * t18 + t3 * t63 - t4 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121, 0.2e1 * t156 * t141, t152, -0.2e1 * t121, -t153, 0, -pkin(6) * t152 + t102 * t114, pkin(6) * t153 + t104 * t114, 0, ((-t89 - t161) * t103 + (t160 - t88) * t105 * t155) * qJD(1), t11, t2, -t36, t12, -t37, 0, -t115 * t74 + t116 * t66 + t164 * t98 + t177 * t73 + t97 * t35, t116 * t68 - t165 * t98 - t97 * t34 - t73 * t42 + t74 * t77, -t115 * t127 - t164 * t68 - t165 * t66 - t177 * t27 + t26 * t42 + t47 * t34 - t48 * t35 - t8 * t77, t116 * t73 - t127 * t48 + t164 * t26 + t165 * t27 + t8 * t47 + t74 * t97, t11, t2, -t36, t12, -t37, 0, -t115 * t24 + t126 * t66 + t167 * t98 + t177 * t41 + t59 * t35, t126 * t68 + t173 * t98 + t24 * t77 - t59 * t34 - t41 * t42, t115 * t3 + t13 * t42 - t15 * t177 - t167 * t68 + t173 * t66 + t30 * t34 - t31 * t35 - t4 * t77, t126 * t41 + t13 * t167 - t15 * t173 + t24 * t59 + t3 * t31 + t4 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t156 * t107, 0, t136, 0, 0, t102 * t119, t104 * t119, 0, 0, t170, t25, t22, -t170, t23, 0, -t137 * t66 - t28 * t98 + t109 - t168, -t137 * t68 + t29 * t98 + t118 - t85, t34 * t140 + (-t26 + t29) * t66 + (t27 + t28 + t139) * t68 + t163, -t26 * t28 - t27 * t29 + (-t73 * t148 + t172 * t8 - t101 * t127 + (-t101 * t26 + t172 * t27) * qJD(4)) * pkin(3), t170, t25, t22, -t170, t23, 0, -t150 * t68 - t16 * t98 - t50 * t66 + t109 + t159, t17 * t98 - t50 * t68 + t111 - t85, t96 * t34 + (-t13 + t17) * t66 + (t15 + t16 + t139) * t68 + t163, -t13 * t16 - t15 * t17 + t4 * t96 - t41 * t50 + (t101 * t3 + (-t101 * t13 + t15 * t172) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t25, t22, -t170, t23, 0, t27 * t98 - t168 + t8, t26 * t98 + t118, 0, 0, t170, t25, t22, -t170, t23, 0, t15 * t98 + (-t130 - t41) * t68 + t110, -pkin(4) * t175 + t14 * t98 + t111, t34 * pkin(4) - t166 * t66, t166 * t15 + (-t41 * t68 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 + t169, t112 - t171, -t65 - t175, t13 * t68 + t15 * t66 + t24;];
tauc_reg = t7;
