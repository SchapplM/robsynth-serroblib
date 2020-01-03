% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:29
% EndTime: 2019-12-31 16:56:31
% DurationCPUTime: 0.72s
% Computational Cost: add. (820->148), mult. (1800->238), div. (0->0), fcn. (975->4), ass. (0->92)
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t47 = pkin(3) * t35 + pkin(6) * t33;
t18 = t47 * qJD(3) + qJD(2);
t14 = t18 * qJD(1);
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t23 = t33 * pkin(3) - t35 * pkin(6) + qJ(2);
t15 = t23 * qJD(1);
t36 = -pkin(1) - pkin(5);
t25 = t36 * qJD(1) + qJD(2);
t86 = t33 * t25;
t16 = qJD(3) * pkin(6) + t86;
t5 = t34 * t15 - t32 * t16;
t73 = qJD(3) * t35;
t61 = t25 * t73;
t1 = qJD(4) * t5 + t32 * t14 + t34 * t61;
t67 = t33 * qJD(1);
t26 = qJD(4) + t67;
t100 = -t5 * t26 + t1;
t62 = 0.2e1 * qJD(1);
t6 = t32 * t15 + t34 * t16;
t2 = -qJD(4) * t6 + t34 * t14 - t32 * t61;
t99 = -t6 * t26 - t2;
t68 = t32 * qJD(3);
t60 = t33 * t68;
t75 = qJD(1) * t35;
t56 = t34 * t75;
t21 = t56 + t68;
t72 = qJD(4) * t21;
t10 = -qJD(1) * t60 + t72;
t69 = qJD(4) * t35;
t58 = t32 * t69;
t66 = t34 * qJD(3);
t40 = t33 * t66 + t58;
t9 = t40 * qJD(1) - qJD(4) * t66;
t96 = t9 * t32;
t95 = t9 * t33;
t94 = t10 * t33;
t93 = t10 * t34;
t19 = t32 * t75 - t66;
t92 = t19 * t26;
t91 = t21 * t19;
t90 = t21 * t26;
t89 = t26 * t32;
t88 = t26 * t33;
t87 = t26 * t34;
t85 = t33 * t36;
t84 = t34 * t35;
t83 = t35 * t10;
t82 = t35 * t25;
t37 = qJD(3) ^ 2;
t81 = t37 * t33;
t80 = t37 * t35;
t31 = t35 ^ 2;
t79 = t33 ^ 2 - t31;
t38 = qJD(1) ^ 2;
t78 = -t37 - t38;
t77 = qJD(3) * pkin(3);
t76 = t38 * qJ(2);
t74 = qJD(3) * t33;
t71 = qJD(4) * t32;
t70 = qJD(4) * t34;
t65 = qJ(2) * qJD(3);
t64 = qJD(1) * qJD(3);
t63 = t35 * t38 * t33;
t59 = t36 * t73;
t57 = t34 * t69;
t55 = qJD(2) * t62;
t54 = t35 * t64;
t17 = -t77 - t82;
t53 = -t17 + t82;
t52 = t26 + t67;
t51 = t19 + t66;
t50 = -t21 + t68;
t49 = qJD(4) * t33 + qJD(1);
t48 = t33 * t54;
t46 = t32 * t6 + t34 * t5;
t45 = t32 * t5 - t34 * t6;
t44 = t53 * qJD(3);
t43 = qJD(1) * t31 - t88;
t12 = t32 * t23 + t34 * t85;
t11 = t34 * t23 - t32 * t85;
t42 = -pkin(6) * t73 + t17 * t33;
t39 = -t46 * qJD(4) + t1 * t34 - t2 * t32;
t28 = qJ(2) * t55;
t22 = t47 * qJD(1);
t8 = t32 * t22 + t34 * t82;
t7 = t34 * t22 - t32 * t82;
t4 = -t12 * qJD(4) + t34 * t18 - t32 * t59;
t3 = t11 * qJD(4) + t32 * t18 + t34 * t59;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t28, -0.2e1 * t48, 0.2e1 * t79 * t64, -t81, 0.2e1 * t48, -t80, 0, -t36 * t81 + (qJD(2) * t33 + t35 * t65) * t62, -t36 * t80 + (qJD(2) * t35 - t33 * t65) * t62, 0, t28, -t40 * t21 - t9 * t84, (t19 * t34 + t21 * t32) * t74 + (-t93 + t96 + (t19 * t32 - t21 * t34) * qJD(4)) * t35, -t26 * t58 - t95 + (t21 * t35 + t43 * t34) * qJD(3), t32 * t83 + (t57 - t60) * t19, -t26 * t57 - t94 + (-t19 * t35 - t43 * t32) * qJD(3), t52 * t73, t2 * t33 + t4 * t26 + (-t10 * t36 + t17 * t70) * t35 + ((qJD(1) * t11 + t5) * t35 + (t19 * t36 + t53 * t32) * t33) * qJD(3), -t1 * t33 - t3 * t26 + (-t17 * t71 + t36 * t9) * t35 + ((-qJD(1) * t12 - t6) * t35 + (t21 * t36 + t53 * t34) * t33) * qJD(3), -t12 * t10 + t11 * t9 - t3 * t19 - t4 * t21 + t46 * t74 + (t45 * qJD(4) - t1 * t32 - t2 * t34) * t35, t1 * t12 + t2 * t11 + t6 * t3 + t5 * t4 - t44 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t76, 0, 0, 0, 0, 0, 0, t78 * t33, t78 * t35, 0, -t76, 0, 0, 0, 0, 0, 0, -t83 - t49 * t87 + (-t52 * t35 * t32 + t19 * t33) * qJD(3), t35 * t9 + t49 * t89 + (-t26 * t84 + (t21 - t56) * t33) * qJD(3), (-t19 * t73 + t49 * t21 - t94) * t34 + (t49 * t19 + t21 * t73 - t95) * t32, -t45 * t73 - t46 * qJD(1) + (-t44 + t39) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t79 * t38, 0, -t63, 0, 0, -t35 * t76, t33 * t76, 0, 0, t21 * t87 - t96, (-t9 - t92) * t34 + (-t10 - t90) * t32, t26 * t70 + (t33 * t87 + t50 * t35) * qJD(1), t19 * t89 - t93, -t26 * t71 + (-t32 * t88 + t51 * t35) * qJD(1), -t26 * t75, -pkin(3) * t10 - t7 * t26 - t51 * t86 + (-pkin(6) * t87 + t17 * t32) * qJD(4) + (t42 * t32 - t35 * t5) * qJD(1), pkin(3) * t9 + t8 * t26 + t50 * t86 + (pkin(6) * t89 + t17 * t34) * qJD(4) + (t42 * t34 + t35 * t6) * qJD(1), t8 * t19 + t7 * t21 + ((-t10 + t72) * pkin(6) + t100) * t34 + ((qJD(4) * t19 - t9) * pkin(6) + t99) * t32, -t5 * t7 - t6 * t8 + (-t17 - t77) * t86 + t39 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t19 ^ 2 + t21 ^ 2, -t9 + t92, -t91, t90 - t10, t54, -t17 * t21 - t99, t17 * t19 - t100, 0, 0;];
tauc_reg = t13;
