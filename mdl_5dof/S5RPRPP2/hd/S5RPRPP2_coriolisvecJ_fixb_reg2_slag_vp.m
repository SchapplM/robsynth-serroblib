% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:11
% EndTime: 2019-12-31 18:11:14
% DurationCPUTime: 0.57s
% Computational Cost: add. (546->126), mult. (1300->167), div. (0->0), fcn. (621->4), ass. (0->87)
t94 = pkin(3) + pkin(4);
t74 = t94 * qJD(3);
t40 = sin(pkin(7)) * pkin(1) + pkin(6);
t32 = t40 * qJD(1);
t53 = sin(qJ(3));
t29 = t53 * t32;
t54 = cos(qJ(3));
t44 = t54 * qJD(2);
t20 = -t29 + t44;
t77 = qJ(5) * qJD(1);
t8 = t53 * t77 + t20;
t84 = qJD(4) - t8;
t3 = -t74 + t84;
t96 = qJD(4) - t20;
t48 = qJD(3) * qJ(4);
t79 = t53 * qJD(2);
t9 = t79 + (t32 - t77) * t54;
t6 = t48 + t9;
t21 = t54 * t32 + t79;
t14 = t21 + t48;
t10 = -qJD(3) * pkin(3) + t96;
t95 = (t20 + t29) * qJD(3);
t76 = qJD(1) * qJD(3);
t49 = t53 ^ 2;
t50 = t54 ^ 2;
t87 = -t49 + t50;
t24 = 0.2e1 * t87 * t76;
t93 = t54 * t6;
t18 = t21 * qJD(3);
t92 = t18 * t53;
t91 = t18 * t54;
t56 = qJD(3) ^ 2;
t90 = t56 * t53;
t45 = t56 * t54;
t72 = t54 * t76;
t78 = t53 * qJD(4);
t89 = qJ(4) * t72 + qJD(1) * t78;
t43 = qJD(3) * t44;
t47 = qJD(3) * qJD(4);
t88 = t43 + 0.2e1 * t47;
t86 = qJ(4) * t54;
t85 = qJ(5) - t40;
t41 = -cos(pkin(7)) * pkin(1) - pkin(2);
t63 = t53 * qJ(4) - t41;
t19 = t94 * t54 + t63;
t82 = qJD(1) * t19;
t4 = qJD(5) + t82;
t83 = qJD(5) + t4;
t25 = -t54 * pkin(3) - t63;
t15 = qJD(1) * t25;
t33 = qJD(1) * t41;
t81 = qJD(1) * t53;
t80 = qJD(3) * t53;
t75 = qJD(1) * qJD(5);
t73 = t53 * t76;
t27 = t85 * t54;
t61 = -t94 * t53 + t86;
t11 = t61 * qJD(3) + t78;
t5 = -t94 * t73 + t89;
t70 = qJD(1) * t11 + t5;
t69 = t4 + t82;
t68 = 0.2e1 * t15;
t66 = t53 * t72;
t17 = -t32 * t80 + t43;
t65 = pkin(3) * t53 - t86;
t64 = 0.2e1 * qJD(3) * t33;
t7 = t17 + t47;
t16 = pkin(3) * t73 - t89;
t23 = t65 * qJD(3) - t78;
t62 = -qJD(1) * t23 - t40 * t56 - t16;
t59 = t92 + t7 * t54 + (t10 * t54 - t14 * t53) * qJD(3);
t58 = t17 * t54 + t92 + (-t20 * t54 - t21 * t53) * qJD(3);
t57 = qJD(1) ^ 2;
t39 = t53 * t57 * t54;
t37 = qJ(5) * t73;
t36 = -t49 * t57 - t56;
t35 = -0.2e1 * t66;
t34 = 0.2e1 * t66;
t31 = t87 * t57;
t30 = t65 * qJD(1);
t26 = t85 * t53;
t22 = t61 * qJD(1);
t13 = -qJD(3) * t27 - t53 * qJD(5);
t12 = -t54 * qJD(5) + t85 * t80;
t2 = t9 * qJD(3) - t53 * t75;
t1 = -t54 * t75 + t37 + t7;
t28 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t24, t45, t35, -t90, 0, -t40 * t45 + t53 * t64, t40 * t90 + t54 * t64, t58, t58 * t40, t34, t45, -t24, 0, t90, t35, t62 * t54 + t68 * t80, t59, -t68 * t54 * qJD(3) + t62 * t53, t15 * t23 + t16 * t25 + t59 * t40, t34, -t24, -t45, t35, -t90, 0, t70 * t54 + (-t69 * t53 - t13) * qJD(3), t70 * t53 + (t69 * t54 + t12) * qJD(3), -t1 * t54 - t2 * t53 + (-t3 * t54 + t53 * t6) * qJD(3) + (-t12 * t54 - t13 * t53 + (t26 * t54 - t27 * t53) * qJD(3)) * qJD(1), -t1 * t27 + t4 * t11 + t6 * t12 + t3 * t13 + t5 * t19 - t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t45, 0, t17 * t53 - t91 + (-t20 * t53 + t21 * t54) * qJD(3), 0, 0, 0, 0, 0, 0, -t90, 0, t45, -t91 + t7 * t53 + (t10 * t53 + t14 * t54) * qJD(3), 0, 0, 0, 0, 0, 0, -t90, t45, 0, t1 * t53 - t2 * t54 + (t3 * t53 + t93) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t31, 0, t39, 0, 0, -t33 * t81, -t33 * t54 * qJD(1) - t43 + t95, 0, 0, -t39, 0, t31, 0, 0, t39, (-t15 * t53 + t30 * t54) * qJD(1), 0, -t95 + (t15 * t54 + t30 * t53) * qJD(1) + t88, -t18 * pkin(3) + t7 * qJ(4) - t10 * t21 + t96 * t14 - t15 * t30, -t39, t31, 0, t39, 0, 0, (-t21 + t9) * qJD(3) + ((qJ(5) * qJD(3) - t22) * t54 + t83 * t53) * qJD(1), t37 + (-t8 - t29) * qJD(3) + (-t22 * t53 - t83 * t54) * qJD(1) + t88, 0, t1 * qJ(4) - t2 * t94 - t4 * t22 - t3 * t9 + t84 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, t36, t15 * t81 + (-t14 + t21) * qJD(3), 0, 0, 0, 0, 0, 0, -t39, t36, 0, -t83 * t81 + (-t6 + t9) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t73, 0.2e1 * t72, (-t49 - t50) * t57, (t93 + (t3 - t74) * t53) * qJD(1) + t89;];
tauc_reg = t28;
