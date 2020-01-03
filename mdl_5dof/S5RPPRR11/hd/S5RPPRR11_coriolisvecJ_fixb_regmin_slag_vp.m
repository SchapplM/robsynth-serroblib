% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:56
% DurationCPUTime: 0.66s
% Computational Cost: add. (427->136), mult. (930->211), div. (0->0), fcn. (518->4), ass. (0->81)
t33 = sin(qJ(4));
t65 = t33 * qJD(1);
t22 = qJD(5) + t65;
t88 = qJD(5) - t22;
t31 = pkin(1) + qJ(3);
t87 = qJD(1) * t31;
t32 = sin(qJ(5));
t34 = cos(qJ(5));
t35 = cos(qJ(4));
t72 = qJD(1) * t35;
t55 = t34 * t72;
t15 = t32 * qJD(4) + t55;
t61 = qJD(1) * qJD(4);
t53 = t33 * t61;
t4 = t15 * qJD(5) - t32 * t53;
t68 = qJD(5) * t32;
t57 = t35 * t68;
t64 = t34 * qJD(4);
t39 = -t33 * t64 - t57;
t3 = t39 * qJD(1) + qJD(5) * t64;
t86 = t3 * t32;
t24 = qJD(1) * qJ(2) + qJD(3);
t20 = -pkin(6) * qJD(1) + t24;
t62 = qJD(1) * qJD(2);
t70 = qJD(4) * t33;
t6 = t20 * t70 - t35 * t62;
t85 = t6 * t32;
t84 = t6 * t34;
t13 = t32 * t72 - t64;
t83 = t13 * t22;
t82 = t15 * t22;
t81 = t32 * t22;
t80 = t32 * t33;
t79 = t33 * t20;
t78 = t34 * t22;
t77 = t34 * t35;
t76 = t35 * t15;
t75 = t35 * t20;
t29 = t35 ^ 2;
t74 = t33 ^ 2 - t29;
t36 = qJD(4) ^ 2;
t37 = qJD(1) ^ 2;
t73 = -t36 - t37;
t30 = -pkin(6) + qJ(2);
t71 = qJD(4) * t30;
t69 = qJD(4) * t35;
t67 = qJD(5) * t34;
t11 = -qJD(4) * pkin(4) - t75;
t66 = t11 * qJD(5);
t21 = -qJD(2) + t87;
t63 = qJD(2) - t21;
t60 = t22 * t80;
t59 = t33 * t78;
t58 = t22 * t68;
t56 = t22 * t67;
t26 = 0.2e1 * t62;
t54 = 0.2e1 * qJD(3) * qJD(1);
t52 = t35 * t61;
t10 = qJD(4) * pkin(7) + t79;
t51 = t22 * t30 + t10;
t50 = t22 + t65;
t49 = t63 * qJD(1);
t48 = -0.2e1 * t52;
t47 = qJD(5) * t33 + qJD(1);
t46 = pkin(4) * t35 + pkin(7) * t33;
t45 = t32 * t52 + t56;
t17 = t33 * pkin(4) - t35 * pkin(7) + t31;
t8 = t17 * qJD(1) - qJD(2);
t2 = t34 * t10 + t32 * t8;
t44 = t32 * t10 - t34 * t8;
t43 = qJD(2) + t21 + t87;
t42 = qJD(1) * t29 - t22 * t33;
t41 = -t30 * t36 + t54;
t40 = -pkin(7) * t69 + t11 * t33;
t12 = t46 * qJD(4) + qJD(3);
t7 = t20 * t69 + t33 * t62;
t38 = -qJD(2) * t22 - t11 * qJD(4) - qJD(5) * t8 - t7;
t16 = t46 * qJD(1);
t9 = t12 * qJD(1);
t5 = t34 * t9;
t1 = [0, 0, 0, 0, t26, qJ(2) * t26, t26, t54, t24 * qJD(2) + t21 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t31) * qJD(1), t33 * t48, 0.2e1 * t74 * t61, -t36 * t33, -t36 * t35, 0, t41 * t33 + t43 * t69, t41 * t35 - t43 * t70, t39 * t15 + t3 * t77, (t13 * t34 + t15 * t32) * t70 + (-t86 - t34 * t4 + (t13 * t32 - t15 * t34) * qJD(5)) * t35, -t22 * t57 + t3 * t33 + (t42 * t34 + t76) * qJD(4), -t35 * t56 - t4 * t33 + (-t35 * t13 - t42 * t32) * qJD(4), t50 * t69, (t34 * t12 - t17 * t68) * t22 + (t13 * t71 + t38 * t32 - t51 * t67 + t5) * t33 + (t34 * t66 - qJD(2) * t13 - t30 * t4 + t85 + (-t30 * t81 + (t34 * t17 - t30 * t80) * qJD(1) - t44) * qJD(4)) * t35, -(t32 * t12 + t17 * t67) * t22 + (t15 * t71 + (t51 * qJD(5) - t9) * t32 + t38 * t34) * t33 + (-t32 * t66 - qJD(2) * t15 - t30 * t3 + t84 + (-t30 * t78 - (t34 * t33 * t30 + t32 * t17) * qJD(1) - t2) * qJD(4)) * t35; 0, 0, 0, 0, -t37, -t37 * qJ(2), -t37, 0, (-qJD(3) - t24) * qJD(1), 0, 0, 0, 0, 0, t48, 0.2e1 * t53, 0, 0, 0, 0, 0, t58 + (t60 + (t13 - t64) * t35) * qJD(1), (t59 + t76) * qJD(1) + t45; 0, 0, 0, 0, 0, 0, 0, -t37, t49, 0, 0, 0, 0, 0, t73 * t33, t73 * t35, 0, 0, 0, 0, 0, -t35 * t4 - t47 * t78 + (-t50 * t35 * t32 + t33 * t13) * qJD(4), -t35 * t3 + t47 * t81 + (-t22 * t77 + (t15 - t55) * t33) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t37 * t33, -t74 * t37, 0, 0, 0, t35 * t49, -t63 * t65, t15 * t78 + t86, (t3 - t83) * t34 + (-t4 - t82) * t32, (t59 - t76) * qJD(1) + t45, -t58 + (-t60 + (t13 + t64) * t35) * qJD(1), -t22 * t72, -pkin(4) * t4 - t84 - (t34 * t16 - t32 * t75) * t22 - t13 * t79 + (-pkin(7) * t78 + t11 * t32) * qJD(5) + (t40 * t32 + t35 * t44) * qJD(1), -pkin(4) * t3 + t85 + (t32 * t16 + t34 * t75) * t22 - t15 * t79 + (pkin(7) * t81 + t11 * t34) * qJD(5) + (t2 * t35 + t40 * t34) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t13, -t13 ^ 2 + t15 ^ 2, t3 + t83, -t4 + t82, t52, -t11 * t15 - t88 * t2 - t32 * t7 + t5, t11 * t13 - t32 * t9 - t34 * t7 + t88 * t44;];
tauc_reg = t1;
