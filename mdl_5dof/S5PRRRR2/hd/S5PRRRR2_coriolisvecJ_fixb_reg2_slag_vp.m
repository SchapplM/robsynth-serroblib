% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:51
% EndTime: 2019-12-05 17:04:53
% DurationCPUTime: 0.58s
% Computational Cost: add. (918->99), mult. (1944->145), div. (0->0), fcn. (1046->6), ass. (0->80)
t40 = cos(qJ(4));
t33 = qJD(2) + qJD(3);
t41 = cos(qJ(3));
t78 = pkin(2) * qJD(2);
t70 = t41 * t78;
t52 = t33 * pkin(3) + t70;
t37 = sin(qJ(4));
t38 = sin(qJ(3));
t72 = t38 * t78;
t65 = t37 * t72;
t14 = -t40 * t52 + t65;
t31 = qJD(4) + t33;
t62 = qJD(3) * t70;
t82 = (-qJD(3) - qJD(4)) * t65;
t5 = (qJD(4) * t52 + t62) * t40 + t82;
t68 = -t14 * t31 - t5;
t84 = t38 * t40;
t55 = t37 * t41 + t84;
t77 = qJD(4) * t37;
t71 = pkin(3) * t77;
t97 = t55 * qJD(3);
t6 = t33 * t71 + (t55 * qJD(4) + t97) * t78;
t39 = cos(qJ(5));
t64 = t40 * t72;
t15 = t37 * t52 + t64;
t11 = t31 * pkin(6) + t15;
t36 = sin(qJ(5));
t7 = t39 * qJD(1) - t36 * t11;
t2 = t7 * qJD(5) + t39 * t5;
t8 = t36 * qJD(1) + t39 * t11;
t3 = -t8 * qJD(5) - t36 * t5;
t98 = t2 * t39 - t3 * t36 + (-t36 * t8 - t39 * t7) * qJD(5);
t34 = t36 ^ 2;
t35 = t39 ^ 2;
t67 = (t34 + t35) * t31;
t95 = t40 * t6;
t4 = t6 * t36;
t74 = qJD(5) * t39;
t94 = t14 * t74 + t4;
t29 = t41 * pkin(2) + pkin(3);
t76 = qJD(4) * t40;
t10 = t29 * t77 + (t38 * t76 + t97) * pkin(2);
t93 = t10 * t31;
t18 = t55 * t78;
t92 = t14 * t18;
t90 = t14 * t37;
t89 = t15 * t31;
t88 = t18 * t31;
t28 = t37 * pkin(3) + pkin(6);
t42 = qJD(5) ^ 2;
t87 = t28 * t42;
t86 = t31 * t36;
t85 = t37 * t38;
t83 = t42 * t36;
t81 = pkin(2) * t84 + t37 * t29;
t80 = t34 - t35;
t75 = qJD(5) * t36;
t30 = t31 ^ 2;
t73 = t36 * t30 * t39;
t69 = t40 * t74;
t12 = t14 * t75;
t63 = t74 * t86;
t61 = t36 * t7 - t39 * t8;
t60 = (-qJD(3) + t33) * t78;
t59 = pkin(2) * qJD(3) * (-qJD(2) - t33);
t58 = pkin(6) * t42 - t89;
t20 = pkin(2) * t85 - t40 * t29;
t57 = t14 * t10 + t6 * t20;
t17 = pkin(6) + t81;
t56 = t17 * t42 + t93;
t54 = t40 * t41 - t85;
t9 = t29 * t76 + (t54 * qJD(3) - t38 * t77) * pkin(2);
t53 = qJD(5) * (t20 * t31 - t9);
t49 = -t70 + (-t31 - t33) * pkin(3);
t32 = t42 * t39;
t23 = -0.2e1 * t63;
t22 = 0.2e1 * t63;
t19 = t54 * t78;
t16 = -0.2e1 * t80 * t31 * qJD(5);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t32, 0, -t61 * qJD(5) + t2 * t36 + t3 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t59, t41 * t59, 0, 0, 0, 0, 0, 0, 0, 0, -t6 - t93, -t9 * t31 - t5, 0, t15 * t9 + t5 * t81 + t57, t22, t16, t32, t23, -t83, 0, t12 + t36 * t53 + (-t56 - t6) * t39, t56 * t36 + t39 * t53 + t94, t9 * t67 + t98, t17 * t98 - t61 * t9 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t60, t41 * t60, 0, 0, 0, 0, 0, 0, 0, 0, t88 - t97 * t78 + (t49 * t37 - t64) * qJD(4), t19 * t31 + (t49 * qJD(4) - t62) * t40 - t82, 0, -t92 - t15 * t19 + (t37 * t5 - t95 + (t15 * t40 + t90) * qJD(4)) * pkin(3), t22, t16, t32, t23, -t83, 0, t12 + (t19 + (-qJD(4) - t31) * t40 * pkin(3)) * t75 + (-t87 - t6 + (t18 - t71) * t31) * t39, t19 * t74 + (t87 - t88) * t36 + (-t31 * t69 + (t37 * t86 - t69) * qJD(4)) * pkin(3) + t94, t98 + (pkin(3) * t76 - t19) * t67, -t92 + t61 * t19 + t98 * t28 + (-t95 + (-t61 * t40 + t90) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 + t89, t68, 0, 0, t22, t16, t32, t23, -t83, 0, (-t58 - t6) * t39, t58 * t36 + t4, t14 * t67 + t98, (-t15 - t61) * t14 + t98 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t80 * t30, 0, t73, 0, 0, t68 * t36, t68 * t39, 0, 0;];
tauc_reg = t1;
