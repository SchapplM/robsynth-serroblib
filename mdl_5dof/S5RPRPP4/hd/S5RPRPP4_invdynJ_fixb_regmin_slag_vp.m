% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:52
% EndTime: 2019-12-31 18:14:54
% DurationCPUTime: 0.81s
% Computational Cost: add. (1315->211), mult. (2406->251), div. (0->0), fcn. (1433->8), ass. (0->114)
t100 = -pkin(1) - pkin(6);
t57 = t100 * qJD(1) + qJD(2);
t123 = -qJ(4) * qJD(1) + t57;
t134 = qJD(1) * qJD(3);
t98 = cos(qJ(3));
t125 = t98 * t134;
t96 = sin(qJ(3));
t136 = t96 * qJDD(1);
t167 = t125 + t136;
t102 = qJD(1) ^ 2;
t97 = sin(qJ(1));
t99 = cos(qJ(1));
t164 = g(1) * t97 - g(2) * t99;
t109 = -t102 * qJ(2) - t164;
t93 = sin(pkin(7));
t94 = cos(pkin(7));
t117 = t93 * t98 + t94 * t96;
t163 = t117 * qJD(1);
t166 = t163 * qJD(3);
t149 = qJD(1) * t96;
t130 = t93 * t149;
t148 = qJD(1) * t98;
t46 = t94 * t148 - t130;
t41 = t46 ^ 2;
t165 = -t163 ^ 2 - t41;
t121 = g(1) * t99 + g(2) * t97;
t89 = qJD(1) * qJD(2);
t129 = 0.2e1 * t89;
t88 = qJDD(1) * qJ(2);
t162 = -t121 + 0.2e1 * t88 + t129;
t142 = qJ(4) - t100;
t147 = qJD(3) * t96;
t107 = -t98 * qJD(4) + t142 * t147;
t124 = t142 * t98;
t37 = -qJD(3) * t124 - t96 * qJD(4);
t20 = -t107 * t94 + t93 * t37;
t21 = t107 * t93 + t94 * t37;
t126 = t96 * t134;
t135 = t98 * qJDD(1);
t131 = t93 * t135 + t167 * t94;
t27 = t93 * t126 - t131;
t119 = -t94 * t135 + t93 * t136;
t28 = t119 + t166;
t54 = t142 * t96;
t30 = t94 * t124 - t93 * t54;
t31 = -t93 * t124 - t94 * t54;
t161 = -t163 * t21 + t20 * t46 + t31 * t27 - t30 * t28;
t133 = qJD(1) * qJD(4);
t56 = t100 * qJDD(1) + qJDD(2);
t51 = t98 * t56;
t16 = -t98 * t133 - t57 * t147 + qJDD(3) * pkin(3) + t51 + (t126 - t135) * qJ(4);
t146 = qJD(3) * t98;
t22 = t123 * t146 + (-qJ(4) * qJDD(1) - t133 + t56) * t96;
t4 = t94 * t16 - t93 * t22;
t128 = -qJDD(5) + t4;
t53 = pkin(3) * t149 + qJD(1) * qJ(2) + qJD(4);
t19 = pkin(4) * t163 - t46 * qJ(5) + t53;
t86 = qJ(3) + pkin(7);
t77 = sin(t86);
t78 = cos(t86);
t160 = g(3) * t77 - t164 * t78 - t19 * t46 + t128;
t158 = g(3) * t96;
t82 = t96 * pkin(3);
t5 = t93 * t16 + t94 * t22;
t38 = t123 * t96;
t155 = t93 * t38;
t34 = t94 * t38;
t39 = t123 * t98;
t36 = qJD(3) * pkin(3) + t39;
t18 = t93 * t36 + t34;
t154 = t99 * pkin(1) + t97 * qJ(2);
t92 = t98 ^ 2;
t152 = t96 ^ 2 - t92;
t151 = qJ(2) + t82;
t150 = pkin(1) * qJDD(1);
t145 = qJDD(3) * pkin(4);
t143 = t53 * qJD(1);
t141 = pkin(3) * t146 + qJD(2);
t24 = t94 * t39 - t155;
t140 = qJD(5) - t24;
t101 = qJD(3) ^ 2;
t139 = -t101 - t102;
t137 = qJDD(3) * t96;
t132 = qJDD(3) * qJ(5) + t5;
t127 = -t97 * pkin(1) + t99 * qJ(2);
t122 = qJDD(2) - t150;
t118 = t77 * pkin(4) - t78 * qJ(5);
t17 = t94 * t36 - t155;
t116 = t167 * pkin(3) + qJDD(4) + t88 + t89;
t95 = -qJ(4) - pkin(6);
t115 = t99 * t82 + t97 * t95 + t127;
t114 = t97 * t82 - t99 * t95 + t154;
t44 = -t94 * t146 + t93 * t147;
t45 = -t93 * t146 - t94 * t147;
t49 = -t93 * t96 + t94 * t98;
t112 = t117 * t27 + t163 * t44 + t49 * t28 - t45 * t46;
t111 = 0.2e1 * qJ(2) * t134 + qJDD(3) * t100;
t106 = -t27 * pkin(4) + t28 * qJ(5) + t116;
t11 = -qJD(3) * pkin(4) + qJD(5) - t17;
t12 = qJD(3) * qJ(5) + t18;
t2 = qJD(3) * qJD(5) + t132;
t3 = -t128 - t145;
t105 = -t11 * t45 + t117 * t2 - t12 * t44 - t3 * t49 - t164;
t104 = t117 * t5 + t17 * t45 - t18 * t44 + t4 * t49 - t164;
t103 = -t100 * t101 + t162;
t79 = qJDD(3) * t98;
t72 = -t94 * pkin(3) - pkin(4);
t68 = t93 * pkin(3) + qJ(5);
t26 = pkin(4) * t117 - t49 * qJ(5) + t151;
t25 = pkin(3) * t148 + t46 * pkin(4) + qJ(5) * t163;
t23 = t93 * t39 + t34;
t8 = -t44 * pkin(4) - t45 * qJ(5) - t49 * qJD(5) + t141;
t1 = -t46 * qJD(5) + t106;
t6 = [qJDD(1), t164, t121, qJDD(2) - 0.2e1 * t150 - t164, t162, -t122 * pkin(1) - g(1) * t127 - g(2) * t154 + (t129 + t88) * qJ(2), t92 * qJDD(1) - 0.2e1 * t96 * t125, 0.2e1 * t152 * t134 - 0.2e1 * t96 * t135, -t101 * t96 + t79, -t101 * t98 - t137, 0, t103 * t96 + t111 * t98, t103 * t98 - t111 * t96, -t104 + t161, -g(1) * t115 - g(2) * t114 + t116 * t151 + t53 * t141 - t17 * t20 + t18 * t21 - t4 * t30 + t5 * t31, -t20 * qJD(3) - t30 * qJDD(3) + t1 * t117 - t121 * t77 + t163 * t8 - t19 * t44 - t26 * t27, -t105 + t161, t21 * qJD(3) + t31 * qJDD(3) - t1 * t49 + t121 * t78 - t19 * t45 + t26 * t28 - t8 * t46, t2 * t31 + t12 * t21 + t1 * t26 + t19 * t8 + t3 * t30 + t11 * t20 - g(1) * (t118 * t99 + t115) - g(2) * (t118 * t97 + t114); 0, 0, 0, qJDD(1), -t102, t122 + t109, 0, 0, 0, 0, 0, t139 * t96 + t79, t139 * t98 - t137, t112, t104 - t143, -qJD(1) * t163 + t45 * qJD(3) + t49 * qJDD(3), t112, qJD(1) * t46 - t44 * qJD(3) + qJDD(3) * t117, -t19 * qJD(1) + t105; 0, 0, 0, 0, 0, 0, t98 * t102 * t96, -t152 * t102, t135, -t136, qJDD(3), t109 * t98 + t158 + t51, g(3) * t98 + (-t109 - t56) * t96, (t18 - t23) * t46 + (-t17 + t24) * t163 + (t27 * t93 + t28 * t94) * pkin(3), t17 * t23 - t18 * t24 + (t158 + t4 * t94 + t5 * t93 + (-t164 - t143) * t98) * pkin(3), t23 * qJD(3) - t25 * t163 + (pkin(4) - t72) * qJDD(3) + t160, t68 * t27 - t72 * t28 + (t12 - t23) * t46 + (t11 - t140) * t163, -g(3) * t78 + t68 * qJDD(3) - t19 * t163 + t25 * t46 - t164 * t77 + (0.2e1 * qJD(5) - t24) * qJD(3) + t132, t2 * t68 + t3 * t72 - t19 * t25 - t11 * t23 - g(3) * (-t118 - t82) + t140 * t12 - t164 * (pkin(3) * t98 + pkin(4) * t78 + qJ(5) * t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t163 * t18 + t17 * t46 + t116 - t121, (t46 - t130) * qJD(3) + t131, t165, t119 + 0.2e1 * t166, t12 * t163 + (-qJD(5) - t11) * t46 + t106 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163 * t46 - qJDD(3), -t119, -t41 - t101, -t12 * qJD(3) - t145 - t160;];
tau_reg = t6;
