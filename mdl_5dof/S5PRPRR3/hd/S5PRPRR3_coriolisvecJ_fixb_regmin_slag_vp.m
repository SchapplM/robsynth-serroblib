% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:32
% EndTime: 2019-12-05 15:47:36
% DurationCPUTime: 0.59s
% Computational Cost: add. (502->108), mult. (1263->169), div. (0->0), fcn. (934->8), ass. (0->84)
t62 = cos(qJ(2));
t86 = t62 * qJD(1);
t43 = qJD(2) * pkin(2) + t86;
t55 = sin(pkin(9));
t59 = sin(qJ(2));
t92 = qJD(1) * t59;
t44 = t55 * t92;
t56 = cos(pkin(9));
t20 = t56 * t43 - t44;
t61 = cos(qJ(4));
t80 = -t61 * pkin(4) - pkin(3);
t16 = t80 * qJD(2) - t20;
t29 = t56 * t86 - t44;
t104 = t16 + t29;
t57 = sin(qJ(5));
t58 = sin(qJ(4));
t60 = cos(qJ(5));
t38 = t57 * t61 + t60 * t58;
t52 = qJD(4) + qJD(5);
t103 = t38 * t52;
t12 = t103 * qJD(2);
t102 = qJD(5) - t52;
t37 = t57 * t58 - t60 * t61;
t67 = t37 * t52;
t45 = t56 * t92;
t21 = t55 * t43 + t45;
t76 = t21 + (pkin(6) + pkin(7)) * qJD(2);
t9 = t61 * qJD(3) - t76 * t58;
t10 = t58 * qJD(3) + t76 * t61;
t101 = t56 * pkin(2);
t49 = t55 * pkin(2) + pkin(6);
t100 = pkin(7) + t49;
t99 = t67 * t52;
t90 = qJD(2) * t61;
t81 = t60 * t90;
t91 = qJD(2) * t58;
t82 = t57 * t91;
t30 = -t81 + t82;
t32 = -t57 * t90 - t60 * t91;
t98 = t32 * t30;
t36 = t55 * t62 + t56 * t59;
t97 = t36 * t52;
t96 = t60 * t10;
t63 = qJD(4) ^ 2;
t95 = t63 * t58;
t94 = t63 * t61;
t93 = t58 ^ 2 - t61 ^ 2;
t89 = qJD(4) * t58;
t88 = qJD(4) * t61;
t85 = qJD(2) * qJD(4);
t84 = pkin(4) * t91;
t83 = pkin(4) * t89;
t8 = qJD(4) * pkin(4) + t9;
t79 = -pkin(4) * t52 - t8;
t78 = t61 * t85;
t77 = qJD(4) * t100;
t18 = -qJD(2) * pkin(3) - t20;
t35 = t55 * t59 - t56 * t62;
t28 = t35 * qJD(2);
t23 = qJD(1) * t28;
t75 = -t18 * qJD(2) + t23;
t27 = t55 * t86 + t45;
t74 = -t27 + t83;
t2 = t9 * qJD(4) - t61 * t23;
t3 = -t10 * qJD(4) + t58 * t23;
t71 = t16 * t32 - t57 * t2 + t60 * t3;
t26 = t36 * qJD(2);
t22 = qJD(1) * t26;
t70 = qJD(2) * t27 - t49 * t63 - t22;
t69 = qJD(4) * (qJD(2) * (-pkin(3) - t101) + t18 + t29);
t11 = qJD(5) * t81 - t52 * t82 + t60 * t78;
t66 = t16 * t30 + (t102 * t10 - t3) * t57;
t64 = qJD(2) ^ 2;
t40 = t80 - t101;
t34 = t100 * t61;
t33 = t100 * t58;
t25 = t61 * t77;
t24 = t58 * t77;
t17 = (t36 * qJD(1) + t83) * qJD(2);
t13 = t103 * t52;
t6 = -t30 ^ 2 + t32 ^ 2;
t5 = -t32 * t52 - t12;
t4 = t30 * t52 + t11;
t1 = [0, 0, -t64 * t59, -t64 * t62, -t20 * t26 - t21 * t28 + t22 * t35 - t23 * t36, 0, 0, 0, 0, 0, t28 * t89 - t36 * t94 + (-t26 * t61 + t35 * t89) * qJD(2), t28 * t88 + t36 * t95 + (t26 * t58 + t35 * t88) * qJD(2), 0, 0, 0, 0, 0, t103 * t28 + t35 * t12 + t26 * t30 + t97 * t67, t103 * t97 + t35 * t11 - t26 * t32 - t28 * t67; 0, 0, 0, 0, t20 * t27 - t21 * t29 + (-t22 * t56 - t23 * t55) * pkin(2), 0.2e1 * t58 * t78, -0.2e1 * t93 * t85, t94, -t95, 0, t58 * t69 + t70 * t61, -t70 * t58 + t61 * t69, t11 * t38 + t32 * t67, t103 * t32 - t11 * t37 - t38 * t12 + t30 * t67, -t99, -t13, 0, (t57 * t24 - t60 * t25 + (t33 * t57 - t34 * t60) * qJD(5)) * t52 + t40 * t12 + t17 * t37 + t74 * t30 + t104 * t103, -(-t60 * t24 - t57 * t25 + (-t33 * t60 - t34 * t57) * qJD(5)) * t52 + t40 * t11 + t17 * t38 - t74 * t32 - t104 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94, 0, 0, 0, 0, 0, -t13, t99; 0, 0, 0, 0, 0, -t58 * t64 * t61, t93 * t64, 0, 0, 0, t75 * t58, t75 * t61, -t98, t6, t4, t5, 0, -(-t57 * t9 - t96) * t52 - t30 * t84 + (t79 * t57 - t96) * qJD(5) + t71, t32 * t84 + (t79 * qJD(5) + t9 * t52 - t2) * t60 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t6, t4, t5, 0, t71 + t102 * (-t57 * t8 - t96), (-t102 * t8 - t2) * t60 + t66;];
tauc_reg = t1;
