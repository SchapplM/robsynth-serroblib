% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:01
% EndTime: 2021-01-15 11:24:05
% DurationCPUTime: 0.58s
% Computational Cost: add. (1009->150), mult. (2105->202), div. (0->0), fcn. (1235->4), ass. (0->94)
t76 = -pkin(1) - pkin(6);
t54 = qJD(1) * t76 + qJD(2);
t125 = -qJ(4) * qJD(1) + t54;
t72 = sin(pkin(7));
t73 = cos(pkin(7));
t74 = sin(qJ(3));
t75 = cos(qJ(3));
t48 = t72 * t75 + t73 * t74;
t123 = t48 * qJD(1);
t34 = qJD(3) * t123;
t107 = qJD(1) * t75;
t108 = qJD(1) * t74;
t92 = t72 * t108;
t44 = t107 * t73 - t92;
t39 = t44 ^ 2;
t124 = -t123 ^ 2 - t39;
t101 = t74 * qJD(4);
t109 = qJ(4) - t76;
t90 = t109 * t75;
t35 = -qJD(3) * t90 - t101;
t100 = t75 * qJD(4);
t106 = qJD(3) * t74;
t82 = t106 * t109 - t100;
t13 = t72 * t35 - t73 * t82;
t14 = t73 * t35 + t72 * t82;
t52 = t109 * t74;
t25 = -t72 * t52 + t73 * t90;
t26 = -t73 * t52 - t72 * t90;
t95 = qJD(1) * qJD(3);
t91 = t75 * t95;
t53 = t73 * t91;
t33 = qJD(3) * t92 - t53;
t122 = -t123 * t14 + t13 * t44 - t25 * t34 + t26 * t33;
t68 = qJD(1) * qJD(2);
t121 = 0.2e1 * t68;
t105 = qJD(3) * t75;
t96 = qJ(4) * qJD(3);
t27 = t54 * t105 + (-t75 * t96 - t101) * qJD(1);
t79 = -t54 * t106 + (t74 * t96 - t100) * qJD(1);
t3 = t72 * t27 - t73 * t79;
t120 = t3 * t25;
t47 = -t72 * t74 + t73 * t75;
t119 = t3 * t47;
t51 = pkin(3) * t108 + qJD(1) * qJ(2) + qJD(4);
t12 = pkin(4) * t123 - t44 * qJ(5) + t51;
t118 = t12 * t44;
t36 = t125 * t74;
t115 = t72 * t36;
t30 = t73 * t36;
t77 = qJD(3) ^ 2;
t114 = t77 * t74;
t113 = t77 * t75;
t4 = t73 * t27 + t72 * t79;
t37 = t125 * t75;
t32 = qJD(3) * pkin(3) + t37;
t11 = t72 * t32 + t30;
t50 = pkin(3) * t91 + t68;
t112 = t74 ^ 2 - t75 ^ 2;
t78 = qJD(1) ^ 2;
t111 = -t77 - t78;
t110 = t78 * qJ(2);
t63 = t74 * pkin(3) + qJ(2);
t104 = t13 * qJD(3);
t103 = t14 * qJD(3);
t102 = t51 * qJD(1);
t16 = t73 * t37 - t115;
t99 = qJD(5) - t16;
t55 = pkin(3) * t105 + qJD(2);
t98 = qJ(2) * qJD(3);
t94 = 0.2e1 * qJD(1);
t93 = pkin(3) * t107;
t10 = t73 * t32 - t115;
t43 = -t105 * t72 - t106 * t73;
t89 = -qJD(1) * t123 + t43 * qJD(3);
t42 = -t105 * t73 + t106 * t72;
t88 = qJD(1) * t44 - t42 * qJD(3);
t15 = t72 * t37 + t30;
t87 = t15 * qJD(3) - t3;
t85 = -t33 * pkin(4) + t34 * qJ(5) + t50;
t84 = t123 * t42 + t48 * t33 + t47 * t34 - t43 * t44;
t83 = t53 + (t44 - t92) * qJD(3);
t2 = qJD(3) * qJD(5) + t4;
t7 = -qJD(3) * pkin(4) + qJD(5) - t10;
t8 = qJD(3) * qJ(5) + t11;
t81 = t2 * t48 - t8 * t42 - t7 * t43 - t119;
t80 = t10 * t43 - t11 * t42 + t4 * t48 - t119;
t64 = -t73 * pkin(3) - pkin(4);
t61 = t72 * pkin(3) + qJ(5);
t19 = t48 * pkin(4) - t47 * qJ(5) + t63;
t18 = 0.2e1 * t34;
t17 = t44 * pkin(4) + qJ(5) * t123 + t93;
t6 = -t42 * pkin(4) - t43 * qJ(5) - t47 * qJD(5) + t55;
t1 = -t44 * qJD(5) + t85;
t5 = [0, 0, 0, 0, t121, qJ(2) * t121, -0.2e1 * t74 * t91, 0.2e1 * t112 * t95, -t114, -t113, 0, -t76 * t114 + (qJD(2) * t74 + t75 * t98) * t94, -t76 * t113 + (qJD(2) * t75 - t74 * t98) * t94, t123 * t55 - t63 * t33 - t51 * t42 + t50 * t48 - t104, -t63 * t34 + t51 * t43 + t55 * t44 + t50 * t47 - t103, -t80 + t122, -t10 * t13 + t11 * t14 + t4 * t26 + t50 * t63 + t51 * t55 + t120, t1 * t48 - t12 * t42 + t123 * t6 - t19 * t33 - t104, -t81 + t122, -t1 * t47 - t12 * t43 + t19 * t34 - t6 * t44 + t103, t1 * t19 + t12 * t6 + t7 * t13 + t8 * t14 + t2 * t26 + t120; 0, 0, 0, 0, -t78, -t110, 0, 0, 0, 0, 0, t111 * t74, t111 * t75, t89, -t88, t84, t80 - t102, t89, t84, t88, -t12 * qJD(1) + t81; 0, 0, 0, 0, 0, 0, t75 * t78 * t74, -t112 * t78, 0, 0, 0, -t75 * t110, t74 * t110, -t123 * t93 - t51 * t44 + t87, t16 * qJD(3) + t123 * t51 - t44 * t93 - t4, (t11 - t15) * t44 + (-t10 + t16) * t123 + (t33 * t72 + t34 * t73) * pkin(3), t10 * t15 - t11 * t16 + (-t102 * t75 - t3 * t73 + t4 * t72) * pkin(3), -t123 * t17 - t118 + t87, t61 * t33 - t64 * t34 + (-t15 + t8) * t44 + (t7 - t99) * t123, -t12 * t123 + t17 * t44 + (0.2e1 * qJD(5) - t16) * qJD(3) + t4, -t12 * t17 - t7 * t15 + t2 * t61 + t3 * t64 + t8 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t18, t124, t10 * t44 + t11 * t123 + t50, t83, t124, t18, t8 * t123 + (-qJD(5) - t7) * t44 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t123, 0, -t39 - t77, -t8 * qJD(3) + t118 + t3;];
tauc_reg = t5;
