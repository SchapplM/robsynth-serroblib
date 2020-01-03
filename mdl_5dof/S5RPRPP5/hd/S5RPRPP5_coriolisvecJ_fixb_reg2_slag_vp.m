% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:44
% DurationCPUTime: 0.56s
% Computational Cost: add. (465->119), mult. (958->153), div. (0->0), fcn. (376->2), ass. (0->91)
t53 = -pkin(3) - pkin(4);
t54 = -pkin(1) - pkin(6);
t37 = t54 * qJD(1) + qJD(2);
t52 = cos(qJ(3));
t80 = qJ(5) * qJD(1);
t16 = (t37 + t80) * t52;
t83 = qJD(4) - t16;
t7 = qJD(3) * t53 + t83;
t98 = t52 * t37;
t26 = qJD(3) * t98;
t47 = qJD(3) * qJD(4);
t18 = t26 + t47;
t91 = qJD(3) * pkin(3);
t67 = -qJD(4) + t91;
t19 = -t67 - t98;
t51 = sin(qJ(3));
t28 = t51 * t37;
t48 = qJD(3) * qJ(4);
t22 = t28 + t48;
t99 = t22 * t52;
t104 = ((-t19 + t98) * t51 - t99) * qJD(3) - t18 * t51;
t76 = 2 * qJD(1);
t70 = t52 * qJ(4) - qJ(2);
t23 = t51 * t53 + t70;
t88 = qJD(1) * t23;
t10 = qJD(5) + t88;
t103 = (qJD(5) + t10) * t52;
t79 = qJD(1) * qJD(3);
t49 = t51 ^ 2;
t50 = t52 ^ 2;
t94 = t49 - t50;
t21 = 0.2e1 * t94 * t79;
t55 = qJD(3) ^ 2;
t97 = t55 * t51;
t96 = t55 * t52;
t87 = qJD(3) * t51;
t25 = t37 * t87;
t71 = qJ(5) * t79;
t95 = t51 * t71 + t25;
t56 = qJD(1) ^ 2;
t93 = t55 + t56;
t92 = qJ(4) * t51;
t90 = t56 * qJ(2);
t89 = qJ(5) + t54;
t31 = t51 * pkin(3) - t70;
t20 = qJD(1) * t31;
t44 = t51 * t80;
t11 = t44 + t22;
t86 = t11 * qJD(3);
t85 = t20 * qJD(1);
t84 = t52 * qJD(4);
t81 = qJ(2) * qJD(3);
t78 = qJD(1) * qJD(5);
t77 = t51 * t78 + t52 * t71 + t26;
t75 = qJD(2) * t76;
t74 = t52 * t79;
t43 = qJD(1) * t84;
t60 = t53 * t52 - t92;
t58 = t60 * qJD(3) - qJD(2);
t2 = t58 * qJD(1) + t43;
t9 = t58 + t84;
t73 = qJD(1) * t9 + t2;
t69 = -t10 - t88;
t68 = qJD(3) * t89;
t65 = t51 * t74;
t3 = t47 + t77;
t5 = -t52 * t78 + t95;
t64 = t3 * t51 + t7 * t87 + (-t5 + t86) * t52;
t63 = pkin(3) * t52 + t92;
t62 = 0.2e1 * qJD(3) * t20;
t59 = t63 * qJD(3) + qJD(2);
t14 = t59 - t84;
t8 = t59 * qJD(1) - t43;
t61 = qJD(1) * t14 - t54 * t55 + t8;
t46 = 0.2e1 * t47;
t45 = qJ(2) * t75;
t41 = t52 * t56 * t51;
t38 = -t50 * t56 - t55;
t36 = -0.2e1 * t65;
t35 = 0.2e1 * t65;
t34 = t93 * t52;
t33 = t93 * t51;
t32 = t94 * t56;
t30 = t89 * t52;
t29 = t89 * t51;
t27 = t63 * qJD(1);
t17 = t60 * qJD(1);
t15 = t28 + t44;
t13 = t51 * qJD(5) + t52 * t68;
t12 = -t52 * qJD(5) + t51 * t68;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t45, t36, t21, -t97, t35, -t96, 0, -t54 * t97 + (qJD(2) * t51 + t52 * t81) * t76, -t54 * t96 + (qJD(2) * t52 - t51 * t81) * t76, 0, t45, t36, -t97, -t21, 0, t96, t35, t61 * t51 + t52 * t62, t104, t51 * t62 - t61 * t52, -t104 * t54 + t20 * t14 + t8 * t31, t36, -t21, t97, t35, -t96, 0, -t73 * t51 + (t69 * t52 - t12) * qJD(3), t73 * t52 + (t69 * t51 + t13) * qJD(3), (-t12 * t52 + t13 * t51 + (t29 * t52 - t30 * t51) * qJD(3)) * qJD(1) + t64, t10 * t9 + t11 * t13 + t7 * t12 + t2 * t23 + t3 * t29 - t5 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t90, 0, 0, 0, 0, 0, 0, -t33, -t34, 0, -t90, 0, 0, 0, 0, 0, 0, -t33, 0, t34, -t104 - t85, 0, 0, 0, 0, 0, 0, -t33, t34, 0, t10 * qJD(1) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t32, 0, -t41, 0, 0, -t52 * t90, t51 * t90, 0, 0, t41, 0, t32, 0, 0, -t41, (-t20 * t52 - t27 * t51) * qJD(1), ((t22 - t48) * t52 + (t19 + t67) * t51) * qJD(1), t46 + (-t20 * t51 + t27 * t52) * qJD(1), t18 * qJ(4) + t22 * qJD(4) - t20 * t27 + (-t99 + (-t19 - t91) * t51) * t37, t41, t32, 0, -t41, 0, 0, t15 * qJD(3) + (t17 * t51 + t103) * qJD(1) - t95, -t16 * qJD(3) + t46 + (t10 * t51 - t17 * t52) * qJD(1) + t77, (-t11 + t15 + t48) * t52 * qJD(1), t3 * qJ(4) - t10 * t17 + t83 * t11 - t7 * t15 + t5 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t38, -t22 * qJD(3) + t52 * t85 + t25, 0, 0, 0, 0, 0, 0, t41, t38, 0, -qJD(1) * t103 - t86 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t74, -0.2e1 * t51 * t79, (-t49 - t50) * t56, t43 + (-t11 * t51 + t52 * t7 + t58) * qJD(1);];
tauc_reg = t1;
