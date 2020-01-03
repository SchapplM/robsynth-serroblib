% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:57
% DurationCPUTime: 0.96s
% Computational Cost: add. (1679->180), mult. (4567->261), div. (0->0), fcn. (3254->6), ass. (0->119)
t111 = sin(qJ(4));
t113 = cos(qJ(4));
t109 = sin(pkin(7));
t110 = cos(pkin(7));
t112 = sin(qJ(2));
t114 = cos(qJ(2));
t91 = t109 * t114 + t110 * t112;
t144 = qJD(1) * t91;
t143 = t110 * t114;
t130 = qJD(1) * t143;
t137 = qJD(1) * t112;
t79 = t109 * t137 - t130;
t120 = t111 * t79 - t113 * t144;
t81 = t91 * qJD(2);
t68 = qJD(1) * t81;
t134 = qJD(1) * qJD(2);
t129 = t112 * t134;
t100 = t109 * t129;
t128 = t114 * t134;
t69 = t110 * t128 - t100;
t117 = t120 * qJD(4) - t111 * t69 - t113 * t68;
t106 = qJD(2) + qJD(4);
t146 = t120 * t106;
t166 = t117 - t146;
t135 = qJD(4) * t113;
t136 = qJD(4) * t111;
t119 = -t111 * t68 + t113 * t69 - t79 * t135 - t136 * t144;
t34 = -t111 * t144 - t113 * t79;
t145 = t34 * t106;
t165 = t119 - t145;
t152 = t120 ^ 2;
t153 = t34 ^ 2;
t164 = t152 - t153;
t151 = t34 * t120;
t149 = -qJ(3) - pkin(5);
t125 = qJD(2) * t149;
t75 = t114 * qJD(3) + t112 * t125;
t54 = t75 * qJD(1);
t76 = -t112 * qJD(3) + t114 * t125;
t55 = t76 * qJD(1);
t23 = -t109 * t54 + t110 * t55;
t16 = -t69 * pkin(6) + t23;
t24 = t109 * t55 + t110 * t54;
t17 = -t68 * pkin(6) + t24;
t158 = t144 * pkin(6);
t99 = t149 * t114;
t96 = qJD(1) * t99;
t85 = t109 * t96;
t148 = qJD(2) * pkin(2);
t98 = t149 * t112;
t95 = qJD(1) * t98;
t89 = t95 + t148;
t36 = t110 * t89 + t85;
t21 = qJD(2) * pkin(3) - t158 + t36;
t159 = t79 * pkin(6);
t147 = t110 * t96;
t37 = t109 * t89 - t147;
t22 = t37 - t159;
t6 = t111 * t21 + t113 * t22;
t2 = -t6 * qJD(4) - t111 * t17 + t113 * t16;
t105 = -t114 * pkin(2) - pkin(1);
t138 = qJD(1) * t105;
t97 = qJD(3) + t138;
t42 = t79 * pkin(3) + t97;
t163 = t120 * t42 + t2;
t1 = (qJD(4) * t21 + t17) * t113 + t111 * t16 - t22 * t136;
t162 = -t42 * t34 - t1;
t161 = -0.2e1 * t134;
t160 = t144 ^ 2;
t40 = -t109 * t95 + t147;
t25 = t40 + t159;
t41 = t110 * t95 + t85;
t26 = t41 - t158;
t104 = t110 * pkin(2) + pkin(3);
t155 = pkin(2) * t109;
t73 = t113 * t104 - t111 * t155;
t157 = t73 * qJD(4) - t111 * t25 - t113 * t26;
t74 = t111 * t104 + t113 * t155;
t156 = -t74 * qJD(4) + t111 * t26 - t113 * t25;
t154 = pkin(2) * t112;
t150 = t144 * t79;
t28 = t109 * t76 + t110 * t75;
t45 = t109 * t98 - t110 * t99;
t116 = qJD(1) ^ 2;
t142 = t114 * t116;
t115 = qJD(2) ^ 2;
t141 = t115 * t112;
t140 = t115 * t114;
t139 = t112 ^ 2 - t114 ^ 2;
t133 = t112 * t148;
t132 = pkin(2) * t137;
t131 = t112 * t142;
t102 = pkin(2) * t129;
t43 = t68 * pkin(3) + t102;
t27 = -t109 * t75 + t110 * t76;
t44 = t109 * t99 + t110 * t98;
t124 = 0.2e1 * t144;
t122 = pkin(1) * t161;
t121 = t112 * t128;
t29 = -t91 * pkin(6) + t44;
t90 = t109 * t112 - t143;
t30 = -t90 * pkin(6) + t45;
t9 = -t111 * t30 + t113 * t29;
t10 = t111 * t29 + t113 * t30;
t39 = -t111 * t90 + t113 * t91;
t84 = t90 * qJD(2);
t77 = t79 ^ 2;
t57 = t90 * pkin(3) + t105;
t50 = t81 * pkin(3) + t133;
t49 = pkin(3) * t144 + t132;
t38 = t111 * t91 + t113 * t90;
t19 = -t81 * pkin(6) + t28;
t18 = t84 * pkin(6) + t27;
t14 = t39 * qJD(4) - t111 * t84 + t113 * t81;
t13 = t111 * t81 + t113 * t84 + t90 * t135 + t91 * t136;
t5 = -t111 * t22 + t113 * t21;
t4 = -t10 * qJD(4) - t111 * t19 + t113 * t18;
t3 = t9 * qJD(4) + t111 * t18 + t113 * t19;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121, t139 * t161, t140, -0.2e1 * t121, -t141, 0, -pkin(5) * t140 + t112 * t122, pkin(5) * t141 + t114 * t122, 0, 0, -t144 * t84 + t69 * t91, -t144 * t81 - t91 * t68 - t69 * t90 + t84 * t79, -t84 * qJD(2), t68 * t90 + t79 * t81, -t81 * qJD(2), 0, t105 * t68 + t97 * t81 + (t27 + (qJD(1) * t90 + t79) * t154) * qJD(2), t105 * t69 - t97 * t84 + (t124 * t154 - t28) * qJD(2), -t144 * t27 - t23 * t91 - t24 * t90 - t28 * t79 + t36 * t84 - t37 * t81 - t44 * t69 - t45 * t68, t23 * t44 + t24 * t45 + t36 * t27 + t37 * t28 + (t97 + t138) * t133, t119 * t39 + t120 * t13, t117 * t39 - t119 * t38 + t120 * t14 - t13 * t34, -t13 * t106, -t117 * t38 - t14 * t34, -t14 * t106, 0, t4 * t106 - t117 * t57 + t42 * t14 - t34 * t50 + t43 * t38, -t3 * t106 + t119 * t57 - t120 * t50 - t42 * t13 + t43 * t39, -t1 * t38 + t10 * t117 - t119 * t9 + t120 * t4 + t5 * t13 - t6 * t14 - t2 * t39 + t3 * t34, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4 + t42 * t50 + t43 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t139 * t116, 0, t131, 0, 0, t116 * pkin(1) * t112, pkin(1) * t142, 0, 0, t150, -t77 + t160, -t100 + (t79 + t130) * qJD(2), -t150, 0, 0, -t40 * qJD(2) - t79 * t132 - t144 * t97 + t23, t41 * qJD(2) - t132 * t144 + t97 * t79 - t24, (t37 + t40) * t144 + (-t36 + t41) * t79 + (-t109 * t68 - t110 * t69) * pkin(2), -t36 * t40 - t37 * t41 + (t109 * t24 + t110 * t23 - t97 * t137) * pkin(2), t151, t164, t165, -t151, t166, 0, t156 * t106 + t34 * t49 + t163, -t157 * t106 + t120 * t49 + t162, -t119 * t73 + t74 * t117 + (t157 + t5) * t34 + (t156 - t6) * t120, t1 * t74 + t156 * t5 + t157 * t6 + t2 * t73 - t42 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124 * qJD(2), -t100 + (-t79 + t130) * qJD(2), -t77 - t160, t144 * t36 + t37 * t79 + t102, 0, 0, 0, 0, 0, 0, -t117 - t146, t119 + t145, -t152 - t153, -t120 * t5 - t6 * t34 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t164, t165, -t151, t166, 0, t6 * t106 + t163, t5 * t106 + t162, 0, 0;];
tauc_reg = t7;
