% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:16
% EndTime: 2021-01-15 14:48:19
% DurationCPUTime: 0.44s
% Computational Cost: add. (407->99), mult. (1083->136), div. (0->0), fcn. (742->6), ass. (0->71)
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t83 = -t39 * t36 + t41 * t37;
t20 = t83 * qJD(1);
t25 = t41 * t36 + t39 * t37;
t23 = t25 * qJD(3);
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t21 = t25 * qJD(1);
t15 = qJD(3) * pkin(6) + t21;
t60 = qJ(5) * qJD(3);
t52 = t15 + t60;
t46 = t52 * t40;
t7 = t38 * qJD(2) + t46;
t33 = t40 * qJD(2);
t6 = -t52 * t38 + t33;
t67 = qJD(4) * pkin(4);
t5 = t6 + t67;
t82 = t5 - t6;
t34 = t38 ^ 2;
t81 = pkin(4) * t34;
t80 = t40 * pkin(4);
t43 = qJD(3) ^ 2;
t77 = t40 * t43;
t42 = qJD(4) ^ 2;
t74 = t42 * t38;
t73 = t42 * t40;
t72 = -qJ(5) - pkin(6);
t62 = t38 * qJD(4);
t65 = t21 * qJD(3);
t71 = t20 * t62 + t40 * t65;
t17 = qJD(1) * t23;
t35 = t40 ^ 2;
t70 = t34 - t35;
t69 = t34 + t35;
t68 = qJD(3) * pkin(3);
t14 = -t20 - t68;
t66 = t14 * qJD(3);
t22 = t83 * qJD(3);
t64 = t22 * qJD(3);
t61 = t40 * qJD(4);
t58 = qJD(3) * qJD(4);
t57 = qJD(4) * qJD(2);
t54 = t38 * t58;
t9 = pkin(4) * t54 + t17;
t56 = 0.2e1 * t58;
t32 = -pkin(3) - t80;
t55 = -pkin(6) * t42 - t17;
t53 = qJD(4) * t72;
t51 = t40 * t56;
t16 = qJD(1) * t22;
t50 = t16 + t57;
t49 = -qJD(3) * qJD(5) - t16;
t48 = t38 * t5 - t40 * t7;
t47 = qJD(4) * (t14 - t68);
t8 = qJD(3) * t32 + qJD(5) - t20;
t45 = (-qJD(5) - t8) * qJD(3) - t50;
t10 = t15 * t62;
t3 = -qJ(5) * t54 - t10 + (-t49 + t57) * t40;
t4 = -t7 * qJD(4) + t49 * t38;
t44 = t3 * t40 - t4 * t38 + (-t38 * t7 - t40 * t5) * qJD(4);
t27 = t72 * t40;
t26 = t72 * t38;
t19 = -t38 * qJD(5) + t40 * t53;
t18 = t40 * qJD(5) + t38 * t53;
t12 = t20 * t61;
t2 = -t22 * t62 - t25 * t73 + (-t23 * t40 - t62 * t83) * qJD(3);
t1 = -t22 * t61 + t25 * t74 + (t23 * t38 - t61 * t83) * qJD(3);
t11 = [0, 0, 0, -t23 * qJD(3), -t64, 0, 0, 0, 0, 0, t2, t1, t2, t1, t69 * t64, -t48 * t22 + t8 * t23 + t44 * t25 - t83 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73, -t74, -t73, 0, -t48 * qJD(4) + t3 * t38 + t4 * t40; 0, 0, 0, t65 - t17, 0, t38 * t51, -t70 * t56, t73, -t74, 0, t38 * t47 + t55 * t40 + t71, t12 + t40 * t47 + (-t55 - t65) * t38, -t9 * t40 + (t19 + (t8 + (t32 - t80) * qJD(3)) * t38) * qJD(4) + t71, t12 + (t9 - t65) * t38 + (t40 * t8 - t18 + (t32 * t40 + t81) * qJD(3)) * qJD(4), (t18 * t40 - t19 * t38 - t69 * t20 + (-t26 * t40 + t27 * t38) * qJD(4)) * qJD(3) + t44, t7 * t18 + t5 * t19 + t4 * t26 - t3 * t27 + t9 * t32 + (pkin(4) * t62 - t21) * t8 + t48 * t20; 0, 0, 0, 0, 0, -t38 * t77, t70 * t43, 0, 0, 0, (-t16 - t66) * t38, t10 + (-t38 * t15 + t33) * qJD(4) + (-t50 - t66) * t40, (t7 - t46) * qJD(4) + (pkin(4) * t77 + t45) * t38, -t43 * t81 + t10 + (t38 * t60 + t6) * qJD(4) + t45 * t40, (-t67 + t82) * t40 * qJD(3), t82 * t7 + (-t8 * t38 * qJD(3) + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t54, t51, -t69 * t43, t48 * qJD(3) + t9;];
tauc_reg = t11;
