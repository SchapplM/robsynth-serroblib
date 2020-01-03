% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:24
% EndTime: 2019-12-31 17:32:26
% DurationCPUTime: 0.38s
% Computational Cost: add. (556->95), mult. (1403->140), div. (0->0), fcn. (1033->6), ass. (0->67)
t50 = cos(qJ(3));
t61 = t50 * qJD(2);
t57 = qJD(4) - t61;
t46 = cos(pkin(8));
t49 = cos(qJ(5));
t45 = sin(pkin(8));
t47 = sin(qJ(5));
t76 = t47 * t45;
t29 = -t49 * t46 + t76;
t32 = (qJD(4) + t61) * qJD(3);
t81 = t29 * t32;
t48 = sin(qJ(3));
t18 = t29 * t48;
t30 = t49 * t45 + t47 * t46;
t69 = qJD(3) * t30;
t80 = t69 ^ 2;
t72 = pkin(6) + qJ(4);
t33 = t72 * t45;
t34 = t72 * t46;
t10 = -t47 * t33 + t49 * t34;
t79 = t10 * qJD(5) + t57 * t30;
t52 = t29 * t50;
t9 = -t49 * t33 - t47 * t34;
t78 = -qJD(2) * t52 + t29 * qJD(4) - t9 * qJD(5);
t67 = qJD(3) * t46;
t59 = t49 * t67;
t60 = qJD(3) * t76;
t22 = -t59 + t60;
t77 = t69 * t22;
t51 = qJD(3) ^ 2;
t74 = t51 * t48;
t73 = t51 * t50;
t71 = t45 ^ 2 + t46 ^ 2;
t70 = qJD(3) * pkin(3);
t41 = -t46 * pkin(4) - pkin(3);
t68 = qJD(3) * t41;
t66 = qJD(3) * t48;
t26 = t29 * qJD(5);
t65 = t26 * qJD(5);
t27 = t30 * qJD(5);
t64 = t27 * qJD(5);
t63 = t46 * qJD(1);
t62 = t48 * qJD(2);
t58 = t71 * t32;
t37 = qJD(3) * qJ(4) + t62;
t20 = -t45 * qJD(1) + t46 * t37;
t11 = -t63 + (-pkin(6) * qJD(3) - t37) * t45;
t12 = pkin(6) * t67 + t20;
t3 = t49 * t11 - t47 * t12;
t4 = t47 * t11 + t49 * t12;
t36 = qJD(5) * t59;
t15 = qJD(5) * t60 - t36;
t56 = t29 * t15 - t27 * t69;
t16 = qJD(3) * t27;
t55 = -t30 * t16 + t26 * t22;
t54 = (-t45 * t37 - t63) * t45 - t20 * t46;
t53 = t30 * t32;
t17 = t30 * t48;
t42 = qJD(3) * t62;
t35 = t57 - t70;
t28 = t57 + t68;
t21 = t22 ^ 2;
t8 = qJD(5) * t18 - t50 * t69;
t7 = -qJD(3) * t52 - qJD(5) * t17;
t2 = -t4 * qJD(5) - t53;
t1 = t3 * qJD(5) - t81;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t65, -t55 + t56, -t1 * t30 + t2 * t29 + t4 * t26 + t3 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t74, t45 * t74, t71 * t73, t48 * t58 + (t35 * t48 + (-t54 - t62) * t50) * qJD(3), 0, 0, 0, 0, 0, 0, t8 * qJD(5) - t50 * t16 + t22 * t66, -t7 * qJD(5) + t50 * t15 + t66 * t69, -t17 * t15 + t18 * t16 - t7 * t22 - t69 * t8, -t1 * t18 - t2 * t17 + t3 * t8 + t4 * t7 + (t28 - t61) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t57 * t71 + t58, -t54 * qJD(4) + qJ(4) * t58 + (t54 * t50 + (-t35 - t70) * t48) * qJD(2), -t15 * t30 - t26 * t69, t55 + t56, -t65, t16 * t29 + t22 * t27, -t64, 0, t41 * t16 + t28 * t27 - t79 * qJD(5) + (qJD(3) * t29 - t22) * t62, t78 * qJD(5) - t41 * t15 - t28 * t26, -t1 * t29 - t10 * t16 + t9 * t15 - t2 * t30 + t78 * t22 + t3 * t26 - t4 * t27 + t69 * t79, t1 * t10 + t2 * t9 - t78 * t4 - t79 * t3 + (-t28 + t68) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 * t51, t54 * qJD(3) + t42, 0, 0, 0, 0, 0, 0, 0.2e1 * t69 * qJD(5), t36 + (-t22 - t60) * qJD(5), -t21 - t80, t4 * t22 + t3 * t69 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t21 + t80, t36 + (t22 - t60) * qJD(5), -t77, 0, 0, -t28 * t69 - t53, t28 * t22 + t81, 0, 0;];
tauc_reg = t5;
