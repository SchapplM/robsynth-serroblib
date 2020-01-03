% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:44
% EndTime: 2019-12-31 16:52:46
% DurationCPUTime: 1.04s
% Computational Cost: add. (1981->95), mult. (3941->144), div. (0->0), fcn. (4479->6), ass. (0->93)
t153 = qJD(3) + qJD(4);
t141 = cos(qJ(3));
t82 = sin(pkin(7));
t83 = cos(pkin(7));
t85 = sin(qJ(3));
t72 = -t141 * t83 + t85 * t82;
t84 = sin(qJ(4));
t136 = t84 * t72;
t140 = cos(qJ(4));
t74 = t141 * t82 + t85 * t83;
t67 = t140 * t74;
t145 = t67 - t136;
t138 = t145 ^ 2;
t55 = t140 * t72 + t84 * t74;
t139 = t55 ^ 2;
t150 = -t138 + t139;
t152 = t150 * qJD(1);
t113 = t55 * qJD(4);
t27 = -t55 * qJD(3) - t113;
t135 = pkin(5) + qJ(2);
t93 = t135 * t83;
t94 = t135 * t82;
t60 = t141 * t94 + t85 * t93;
t44 = -t74 * pkin(6) - t60;
t61 = t141 * t93 - t85 * t94;
t45 = -t72 * pkin(6) + t61;
t25 = t140 * t44 - t84 * t45;
t149 = t153 * t25;
t115 = t145 * qJD(1);
t148 = t55 * t115;
t147 = qJD(1) * t55;
t146 = qJD(2) * t55;
t86 = t84 * t44;
t96 = -t86 / 0.2e1;
t95 = t140 * t45;
t26 = t95 + t86;
t144 = t74 ^ 2;
t89 = t67 / 0.2e1;
t143 = pkin(3) * t84;
t142 = t74 * pkin(3);
t77 = t82 ^ 2 + t83 ^ 2;
t134 = qJD(3) * pkin(3);
t79 = -t83 * pkin(2) - pkin(1);
t62 = t72 * pkin(3) + t79;
t3 = t62 * t142;
t132 = t3 * qJD(1);
t5 = -t145 * t25 - t26 * t55;
t131 = t5 * qJD(1);
t8 = t96 + t86 / 0.2e1;
t130 = t8 * qJD(1);
t128 = qJD(1) * t62;
t13 = t138 + t139;
t127 = t13 * qJD(1);
t16 = -t55 * t142 - t145 * t62;
t125 = t16 * qJD(1);
t17 = -t142 * t145 + t55 * t62;
t124 = t17 * qJD(1);
t87 = -t55 * t84 / 0.2e1 - t145 * t140 / 0.2e1;
t21 = (-t74 / 0.2e1 + t87) * pkin(3);
t122 = t21 * qJD(1);
t28 = t60 * t74 - t61 * t72;
t121 = t28 * qJD(1);
t34 = 0.2e1 * t89 - t136;
t119 = t34 * qJD(1);
t70 = t72 ^ 2;
t46 = t70 - t144;
t118 = t46 * qJD(1);
t53 = t89 - t67 / 0.2e1;
t117 = t53 * qJD(1);
t116 = t53 * qJD(4);
t114 = t145 * qJD(3);
t110 = t145 * qJD(4);
t59 = t70 + t144;
t109 = t59 * qJD(1);
t108 = t72 * qJD(1);
t69 = t72 * qJD(3);
t107 = t74 * qJD(1);
t106 = t74 * qJD(3);
t76 = t77 * qJ(2);
t105 = t76 * qJD(1);
t104 = t77 * qJD(1);
t101 = t145 * t147;
t100 = t55 * t128;
t99 = t145 * t128;
t98 = t72 * t107;
t97 = t72 * t106;
t92 = t140 * qJD(3);
t91 = t140 * qJD(4);
t90 = qJD(1) * t79 + qJD(2);
t88 = t34 * qJD(4) + t114;
t22 = t142 / 0.2e1 + t87 * pkin(3);
t11 = -t95 + 0.2e1 * t96;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * qJD(2), t76 * qJD(2), -t97, t46 * qJD(3), 0, t97, 0, 0, t79 * t106, -t79 * t69, t59 * qJD(2), t28 * qJD(2), t27 * t145, t153 * t150, 0, (t110 + t114) * t55, 0, 0, -t16 * qJD(3) + t62 * t110, -t17 * qJD(3) - t62 * t113, t13 * qJD(2), t5 * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t105, 0, 0, 0, 0, 0, 0, 0, 0, t109, t121, 0, 0, 0, 0, 0, 0, t116, 0, t127, t22 * qJD(3) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t118, -t69, t98, -t106, 0, -t61 * qJD(3) + t79 * t107, t60 * qJD(3) - t79 * t108, 0, 0, -t101, t152, t27, t148, -t88, 0, -qJD(3) * t26 + t11 * qJD(4) - t125, -t124 - t149, (t140 * t55 - t145 * t84) * t134, t132 + t22 * qJD(2) + (-t140 * t26 + t25 * t84) * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t152, t27, t148, -t34 * qJD(3) - t110, 0, t53 * qJD(2) + t11 * qJD(3) - t26 * qJD(4) + t99, -t100 - t149, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t105, 0, 0, 0, 0, 0, 0, t106, -t69, -t109, -t121, 0, 0, 0, 0, 0, 0, t88, t27, -t127, -t21 * qJD(3) - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t108, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t147, 0, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t147, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t118, 0, -t98, 0, 0, -t90 * t74, t90 * t72, 0, 0, t101, -t152, 0, -t148, -t116, 0, -qJD(2) * t145 - t8 * qJD(4) + t125, t124 + t146, 0, t21 * qJD(2) - t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, t108, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t147, 0, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t143, -pkin(3) * t91, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, 0, -t143 * t153 - t130, (-t92 - t91) * pkin(3), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t152, 0, -t148, t53 * qJD(3), 0, -t34 * qJD(2) + t8 * qJD(3) - t99, t100 + t146, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t147, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, t84 * t134 + t130, pkin(3) * t92, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
