% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:14
% EndTime: 2019-12-31 18:01:15
% DurationCPUTime: 0.31s
% Computational Cost: add. (451->59), mult. (789->100), div. (0->0), fcn. (432->6), ass. (0->56)
t80 = qJD(1) - qJD(4);
t32 = t80 ^ 2;
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t22 = t36 * t39 - t37 * t41;
t67 = t80 * t22;
t23 = t36 * t41 + t37 * t39;
t43 = qJD(5) ^ 2;
t58 = t32 * t23;
t81 = -t23 * t43 - t58;
t42 = -pkin(1) - pkin(2);
t30 = qJD(1) * t42 + qJD(2);
t62 = qJD(1) * qJ(2);
t16 = t30 * t36 + t37 * t62;
t61 = qJD(1) * qJD(2);
t56 = t37 * t61;
t57 = t36 * t61;
t64 = qJD(4) * t41;
t26 = t37 * t30;
t65 = t36 * qJ(2);
t14 = t26 + (-pkin(3) - t65) * qJD(1);
t72 = t39 * t14;
t2 = qJD(4) * t72 + t16 * t64 + t39 * t56 + t41 * t57;
t79 = -t2 - (t16 * t41 + t72) * t80;
t55 = t37 * t42 - t65;
t24 = -pkin(3) + t55;
t25 = qJ(2) * t37 + t36 * t42;
t51 = t24 * t39 + t25 * t41;
t78 = t2 + (qJD(2) * t23 + qJD(4) * t51) * t80;
t1 = -(qJD(4) * t16 + t57) * t39 + t14 * t64 + t41 * t56;
t77 = t80 * pkin(4);
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t73 = t38 * t40;
t71 = t43 * t38;
t70 = t43 * t40;
t66 = t38 ^ 2 - t40 ^ 2;
t63 = qJD(5) * t80;
t60 = 0.2e1 * t61;
t7 = t14 * t41 - t16 * t39;
t3 = -t7 + t77;
t59 = t3 * t80 - t1;
t54 = t63 * t73;
t53 = (-t36 * t62 + t26) * t36 - t16 * t37;
t52 = t24 * t41 - t25 * t39;
t50 = pkin(7) * t43 - t79;
t49 = (-pkin(7) + t51) * t43 - t78;
t48 = qJD(5) * (t3 + t7 + t77);
t5 = -qJD(2) * t22 + qJD(4) * t52;
t47 = qJD(5) * (-t80 * (pkin(4) - t52) - t3 - t5);
t46 = -0.2e1 * qJD(5) * t67;
t44 = qJD(1) ^ 2;
t17 = t66 * t63;
t4 = [0, 0, 0, 0, t60, qJ(2) * t60, 0.2e1 * t57, 0.2e1 * t56, ((t37 * t25 - t36 * t55) * qJD(1) - t53) * qJD(2), 0, t78, t5 * t80 + t1, 0.2e1 * t54, -0.2e1 * t17, -t70, t71, 0, t38 * t47 - t40 * t49, t38 * t49 + t40 * t47; 0, 0, 0, 0, -t44, -t44 * qJ(2), -t36 * t44, -t37 * t44, t53 * qJD(1), 0, -t58, t67 * t80, 0, 0, 0, 0, 0, t38 * t46 + t40 * t81, -t38 * t81 + t40 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t7 * t80 - t1, -0.2e1 * t54, 0.2e1 * t17, t70, -t71, 0, t38 * t48 - t40 * t50, t38 * t50 + t40 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t73, t66 * t32, 0, 0, 0, t59 * t38, t59 * t40;];
tauc_reg = t4;
