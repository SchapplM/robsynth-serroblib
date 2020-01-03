% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:53
% DurationCPUTime: 0.51s
% Computational Cost: add. (690->135), mult. (1955->186), div. (0->0), fcn. (1271->4), ass. (0->91)
t88 = (qJD(1) * qJD(2));
t109 = -2 * t88;
t64 = sin(pkin(6));
t66 = cos(qJ(2));
t65 = sin(qJ(2));
t94 = cos(pkin(6));
t82 = t94 * t65;
t48 = t64 * t66 + t82;
t72 = qJD(1) * t48;
t108 = 0.2e1 * t72;
t105 = t72 ^ 2;
t81 = t94 * t66;
t77 = qJD(1) * t81;
t92 = qJD(1) * t65;
t38 = t64 * t92 - t77;
t35 = t38 ^ 2;
t107 = -t35 - t105;
t106 = -t35 + t105;
t104 = pkin(2) * t65;
t96 = -qJ(3) - pkin(5);
t53 = t96 * t66;
t22 = -t64 * t53 - t96 * t82;
t80 = qJD(2) * t96;
t34 = t66 * qJD(3) + t65 * t80;
t27 = t34 * qJD(1);
t73 = -t65 * qJD(3) + t66 * t80;
t28 = t73 * qJD(1);
t6 = t64 * t27 - t94 * t28;
t103 = t6 * t22;
t61 = -t66 * pkin(2) - pkin(1);
t93 = qJD(1) * t61;
t52 = qJD(3) + t93;
t9 = t38 * pkin(3) - qJ(4) * t72 + t52;
t102 = t9 * t72;
t101 = t72 * t38;
t51 = qJD(1) * t53;
t100 = t64 * t51;
t68 = qJD(1) ^ 2;
t99 = t66 * t68;
t67 = qJD(2) ^ 2;
t98 = t67 * t65;
t97 = t67 * t66;
t7 = t94 * t27 + t64 * t28;
t44 = t94 * t51;
t84 = t96 * t65;
t50 = qJD(1) * t84;
t46 = qJD(2) * pkin(2) + t50;
t19 = t64 * t46 - t44;
t95 = t65 ^ 2 - t66 ^ 2;
t40 = t48 * qJD(2);
t91 = qJD(2) * t40;
t90 = qJD(2) * t65;
t21 = t94 * t50 + t100;
t89 = qJD(4) - t21;
t87 = t65 * t99;
t86 = pkin(2) * t90;
t85 = pkin(2) * t92;
t83 = t65 * t88;
t79 = pkin(1) * t109;
t78 = t66 * t83;
t30 = qJD(1) * t40;
t47 = t64 * t65 - t81;
t76 = t30 * t47 + t38 * t40;
t54 = t64 * t83;
t31 = qJD(2) * t77 - t54;
t57 = pkin(2) * t83;
t75 = t30 * pkin(3) - t31 * qJ(4) + t57;
t20 = t64 * t50 - t44;
t74 = t20 * qJD(2) - t6;
t18 = t94 * t46 + t100;
t11 = t64 * t34 - t94 * t73;
t12 = t94 * t34 + t64 * t73;
t23 = -t94 * t53 + t64 * t84;
t71 = t11 * t72 - t12 * t38 + t22 * t31 - t23 * t30 + t6 * t48;
t43 = qJD(2) * t81 - t64 * t90;
t70 = t48 * t30 + t31 * t47 + t43 * t38 + t40 * t72;
t69 = t108 * qJD(2);
t60 = -t94 * pkin(2) - pkin(3);
t58 = t64 * pkin(2) + qJ(4);
t33 = t43 * qJD(2);
t17 = t47 * pkin(3) - t48 * qJ(4) + t61;
t16 = -t54 + (t77 + t38) * qJD(2);
t15 = -t54 + (t77 - t38) * qJD(2);
t14 = qJD(2) * qJ(4) + t19;
t13 = -qJD(2) * pkin(3) + qJD(4) - t18;
t10 = pkin(3) * t72 + t38 * qJ(4) + t85;
t5 = qJD(2) * qJD(4) + t7;
t4 = t40 * pkin(3) - t43 * qJ(4) - t48 * qJD(4) + t86;
t2 = t31 * t48 + t43 * t72;
t1 = -qJD(4) * t72 + t75;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, t95 * t109, t97, -0.2e1 * t78, -t98, 0, -pkin(5) * t97 + t65 * t79, pkin(5) * t98 + t66 * t79, 0, 0, t2, -t70, t33, t76, -t91, 0, t61 * t30 + t52 * t40 + (-t11 + (qJD(1) * t47 + t38) * t104) * qJD(2), t61 * t31 + t52 * t43 + (t108 * t104 - t12) * qJD(2), -t18 * t43 - t19 * t40 - t7 * t47 + t71, -t18 * t11 + t19 * t12 + t103 + t7 * t23 + (t52 + t93) * t86, t2, t33, t70, 0, t91, t76, -t11 * qJD(2) + t1 * t47 + t17 * t30 + t4 * t38 + t9 * t40, t13 * t43 - t14 * t40 - t5 * t47 + t71, t12 * qJD(2) - t1 * t48 - t17 * t31 - t4 * t72 - t9 * t43, t1 * t17 + t13 * t11 + t14 * t12 + t5 * t23 + t9 * t4 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t95 * t68, 0, t87, 0, 0, t68 * pkin(1) * t65, pkin(1) * t99, 0, 0, t101, t106, t16, -t101, 0, 0, -t38 * t85 - t52 * t72 + t74, t21 * qJD(2) + t52 * t38 - t72 * t85 - t7, (t19 - t20) * t72 + (-t18 + t21) * t38 + (-t30 * t64 - t94 * t31) * pkin(2), t18 * t20 - t19 * t21 + (-t52 * t92 - t94 * t6 + t64 * t7) * pkin(2), t101, t16, -t106, 0, 0, -t101, -t10 * t38 - t102 + t74, -t58 * t30 + t60 * t31 + (t14 - t20) * t72 + (t13 - t89) * t38, t10 * t72 - t9 * t38 + (0.2e1 * qJD(4) - t21) * qJD(2) + t7, -t9 * t10 - t13 * t20 + t89 * t14 + t5 * t58 + t6 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t15, t107, t18 * t72 + t19 * t38 + t57, 0, 0, 0, 0, 0, 0, t69, t107, -t15, t14 * t38 + (-qJD(4) - t13) * t72 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t16, -t105 - t67, -t14 * qJD(2) + t102 + t6;];
tauc_reg = t3;
