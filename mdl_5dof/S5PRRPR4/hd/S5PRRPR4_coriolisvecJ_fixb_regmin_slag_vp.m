% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:22
% EndTime: 2019-12-05 16:23:27
% DurationCPUTime: 0.91s
% Computational Cost: add. (1010->162), mult. (2679->255), div. (0->0), fcn. (1998->8), ass. (0->109)
t101 = cos(qJ(5));
t102 = cos(qJ(3));
t96 = sin(pkin(9));
t97 = cos(pkin(9));
t99 = sin(qJ(3));
t80 = t97 * t102 - t96 * t99;
t73 = t80 * qJD(2);
t62 = t101 * t73;
t81 = t96 * t102 + t97 * t99;
t75 = t81 * qJD(2);
t98 = sin(qJ(5));
t33 = -t98 * t75 + t62;
t93 = qJD(3) + qJD(5);
t143 = t33 * t93;
t133 = qJD(5) * t98;
t74 = t81 * qJD(3);
t65 = qJD(2) * t74;
t126 = qJD(2) * qJD(3);
t120 = t102 * t126;
t122 = t99 * t126;
t66 = t97 * t120 - t96 * t122;
t4 = qJD(5) * t62 + t101 * t66 - t75 * t133 - t98 * t65;
t157 = t4 - t143;
t110 = -t101 * t75 - t98 * t73;
t106 = t110 * qJD(5) - t101 * t65 - t98 * t66;
t144 = t110 * t93;
t156 = t106 - t144;
t155 = t110 * t33;
t154 = t110 ^ 2 - t33 ^ 2;
t146 = t73 * pkin(7);
t100 = sin(qJ(2));
t128 = t100 * qJD(1);
t88 = qJD(2) * pkin(6) + t128;
t118 = qJ(4) * qJD(2) + t88;
t70 = t118 * t102;
t142 = t97 * t70;
t69 = t118 * t99;
t54 = qJD(3) * pkin(3) - t69;
t18 = t96 * t54 + t142;
t12 = t18 + t146;
t123 = -t102 * pkin(3) - pkin(2);
t103 = cos(qJ(2));
t127 = t103 * qJD(1);
t78 = t123 * qJD(2) + qJD(4) - t127;
t40 = -t73 * pkin(4) + t78;
t153 = t12 * t133 - t40 * t33;
t114 = qJD(4) + t127;
t134 = qJD(3) * t99;
t35 = -t88 * t134 + (-qJ(4) * t134 + t114 * t102) * qJD(2);
t129 = qJD(3) * t102;
t36 = -t88 * t129 + (-qJ(4) * t129 - t114 * t99) * qJD(2);
t6 = -t96 * t35 + t97 * t36;
t2 = -t66 * pkin(7) + t6;
t7 = t97 * t35 + t96 * t36;
t3 = -t65 * pkin(7) + t7;
t152 = t101 * t2 + t40 * t110 - t98 * t3;
t151 = -0.2e1 * t126;
t141 = -qJ(4) - pkin(6);
t119 = qJD(3) * t141;
t71 = t102 * qJD(4) + t99 * t119;
t72 = -t99 * qJD(4) + t102 * t119;
t140 = t81 * t127 - t96 * t71 + t97 * t72;
t139 = -t80 * t127 + t97 * t71 + t96 * t72;
t104 = qJD(3) ^ 2;
t105 = qJD(2) ^ 2;
t150 = t100 * (t104 + t105);
t77 = t80 * qJD(3);
t149 = qJD(5) - t93;
t148 = pkin(3) * t134 - t128;
t147 = pkin(3) * t96;
t145 = t75 * pkin(7);
t50 = t96 * t70;
t20 = -t97 * t69 - t50;
t85 = t141 * t99;
t86 = t141 * t102;
t42 = t96 * t85 - t97 * t86;
t79 = pkin(3) * t122 + qJD(2) * t128;
t138 = -t102 ^ 2 + t99 ^ 2;
t137 = qJD(2) * pkin(2);
t136 = t104 * t99;
t135 = qJD(2) * t99;
t132 = t104 * t102;
t130 = qJD(2) * t100;
t17 = t97 * t54 - t50;
t19 = t96 * t69 - t142;
t41 = t97 * t85 + t96 * t86;
t117 = t74 * pkin(4) + t148;
t115 = t103 * t151;
t113 = t77 * pkin(7) + qJD(5) * (t80 * pkin(7) + t42) - t140;
t112 = t74 * pkin(7) - qJD(5) * (-t81 * pkin(7) + t41) - t139;
t11 = qJD(3) * pkin(4) - t145 + t17;
t111 = -t101 * t12 - t98 * t11;
t39 = t101 * t81 + t98 * t80;
t109 = t101 * t80 - t98 * t81;
t108 = qJD(2) * t137;
t107 = -0.2e1 * qJD(3) * t137;
t91 = t97 * pkin(3) + pkin(4);
t68 = t80 * t100;
t67 = t81 * t100;
t58 = -t80 * pkin(4) + t123;
t45 = pkin(3) * t135 + t75 * pkin(4);
t37 = t65 * pkin(4) + t79;
t22 = -t100 * t74 + t103 * t73;
t21 = -t100 * t77 - t103 * t75;
t14 = t20 - t145;
t13 = t19 - t146;
t9 = t39 * qJD(5) + t101 * t74 + t98 * t77;
t8 = t109 * qJD(5) + t101 * t77 - t98 * t74;
t1 = [0, 0, -t105 * t100, -t105 * t103, 0, 0, 0, 0, 0, -t102 * t150 + t99 * t115, t102 * t115 + t99 * t150, -t21 * t75 + t22 * t73 - t68 * t65 + t67 * t66, -t79 * t103 + t78 * t130 + t17 * t21 + t18 * t22 - t6 * t67 + t7 * t68, 0, 0, 0, 0, 0, (t101 * t21 - t98 * t22 + (-t101 * t68 + t67 * t98) * qJD(5)) * t93 - t33 * t130 + t103 * t106, -(t101 * t22 + t98 * t21 + (-t101 * t67 - t68 * t98) * qJD(5)) * t93 - t110 * t130 - t103 * t4; 0, 0, 0, 0, 0.2e1 * t99 * t120, t138 * t151, t132, -t136, 0, -pkin(6) * t132 + t107 * t99, pkin(6) * t136 + t102 * t107, t139 * t73 - t140 * t75 - t17 * t77 - t18 * t74 - t41 * t66 - t42 * t65 - t6 * t81 + t7 * t80, t79 * t123 + t139 * t18 + t140 * t17 + t148 * t78 + t6 * t41 + t7 * t42, -t110 * t8 + t4 * t39, t106 * t39 + t109 * t4 + t110 * t9 + t33 * t8, t8 * t93, -t9 * t93, 0, -t37 * t109 + t40 * t9 - t58 * t106 + (-t113 * t101 + t112 * t98) * t93 - t117 * t33, t37 * t39 + t58 * t4 + t40 * t8 + (t112 * t101 + t113 * t98) * t93 - t117 * t110; 0, 0, 0, 0, -t99 * t105 * t102, t138 * t105, 0, 0, 0, t99 * t108, t102 * t108, (t18 + t19) * t75 + (t17 - t20) * t73 + (-t65 * t96 - t66 * t97) * pkin(3), -t17 * t19 - t18 * t20 + (-t78 * t135 + t6 * t97 + t7 * t96) * pkin(3), t155, t154, t157, t156, 0, -(t101 * t13 - t98 * t14) * t93 + t45 * t33 + ((-t101 * t147 - t91 * t98) * t93 + t111) * qJD(5) + t152, -t101 * t3 - t98 * t2 + (t101 * t14 + t98 * t13) * t93 + t45 * t110 + (-(t101 * t91 - t98 * t147) * t93 - t101 * t11) * qJD(5) + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 ^ 2 - t75 ^ 2, t17 * t75 - t18 * t73 + t79, 0, 0, 0, 0, 0, -t106 - t144, t4 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t154, t157, t156, 0, t149 * t111 + t152, (-t12 * t93 - t2) * t98 + (-t149 * t11 - t3) * t101 + t153;];
tauc_reg = t1;
