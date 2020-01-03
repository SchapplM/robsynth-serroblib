% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:08
% DurationCPUTime: 0.39s
% Computational Cost: add. (288->85), mult. (787->151), div. (0->0), fcn. (538->6), ass. (0->66)
t37 = sin(qJ(4));
t38 = sin(qJ(3));
t40 = cos(qJ(4));
t41 = cos(qJ(3));
t21 = t37 * t41 + t40 * t38;
t34 = qJD(3) + qJD(4);
t83 = t21 * t34;
t5 = t83 * qJD(2);
t79 = qJD(4) - t34;
t62 = qJD(2) * qJD(3);
t84 = -0.2e1 * t62;
t20 = t37 * t38 - t40 * t41;
t47 = t20 * t34;
t39 = sin(qJ(2));
t43 = qJD(3) ^ 2;
t44 = qJD(2) ^ 2;
t82 = (t43 + t44) * t39;
t67 = qJD(3) * t41;
t81 = -qJD(4) * t41 - t67;
t68 = qJD(3) * t38;
t80 = qJD(4) * t38 + t68;
t78 = pkin(5) + pkin(6);
t69 = qJD(2) * t41;
t58 = t40 * t69;
t70 = qJD(2) * t38;
t59 = t37 * t70;
t15 = -t58 + t59;
t17 = -t37 * t69 - t40 * t70;
t77 = t17 * t15;
t64 = t39 * qJD(1);
t28 = qJD(2) * pkin(5) + t64;
t54 = pkin(6) * qJD(2) + t28;
t14 = t54 * t41;
t76 = t40 * t14;
t75 = t43 * t38;
t74 = t43 * t41;
t73 = t38 ^ 2 - t41 ^ 2;
t71 = qJD(2) * pkin(2);
t42 = cos(qJ(2));
t63 = t42 * qJD(1);
t61 = pkin(3) * t70;
t60 = pkin(3) * t68;
t33 = -t41 * pkin(3) - pkin(2);
t13 = t54 * t38;
t12 = qJD(3) * pkin(3) - t13;
t57 = -pkin(3) * t34 - t12;
t56 = qJD(3) * t78;
t55 = t41 * t62;
t52 = t42 * t84;
t50 = qJD(2) * t71;
t10 = -t28 * t67 + (-pkin(6) * t67 - t38 * t63) * qJD(2);
t18 = t33 * qJD(2) - t63;
t9 = -t28 * t68 + (-pkin(6) * t68 + t41 * t63) * qJD(2);
t49 = t40 * t10 + t18 * t17 - t37 * t9;
t46 = -0.2e1 * qJD(3) * t71;
t4 = qJD(4) * t58 - t34 * t59 + t40 * t55;
t45 = t18 * t15 + (t79 * t14 - t10) * t37;
t25 = t78 * t41;
t24 = t78 * t38;
t23 = t41 * t56;
t22 = t38 * t56;
t19 = (t60 + t64) * qJD(2);
t3 = -t15 ^ 2 + t17 ^ 2;
t2 = -t17 * t34 - t5;
t1 = t15 * t34 + t4;
t6 = [0, 0, -t44 * t39, -t44 * t42, 0, 0, 0, 0, 0, t38 * t52 - t41 * t82, t38 * t82 + t41 * t52, 0, 0, 0, 0, 0, -0.2e1 * t42 * t5 + ((t37 * t80 + t40 * t81) * t34 + qJD(2) * t15) * t39, (qJD(2) * t47 - t4) * t42 + (-(t37 * t81 - t40 * t80) * t34 - qJD(2) * t17) * t39; 0, 0, 0, 0, 0.2e1 * t38 * t55, t73 * t84, t74, -t75, 0, -pkin(5) * t74 + t38 * t46, pkin(5) * t75 + t41 * t46, t17 * t47 + t4 * t21, t15 * t47 + t17 * t83 - t4 * t20 - t21 * t5, -t47 * t34, -t83 * t34, 0, (t37 * t22 - t40 * t23 + (t24 * t37 - t25 * t40) * qJD(4)) * t34 + t15 * t60 + t33 * t5 + t19 * t20 + t18 * t83 + (-t39 * t15 + t42 * t83) * qJD(1), -(-t40 * t22 - t37 * t23 + (-t24 * t40 - t25 * t37) * qJD(4)) * t34 - t17 * t60 + t33 * t4 + t19 * t21 - t18 * t47 + (t39 * t17 - t42 * t47) * qJD(1); 0, 0, 0, 0, -t38 * t44 * t41, t73 * t44, 0, 0, 0, t38 * t50, t41 * t50, -t77, t3, t1, t2, 0, -(t37 * t13 - t76) * t34 - t15 * t61 + (t57 * t37 - t76) * qJD(4) + t49, t17 * t61 + (t57 * qJD(4) - t13 * t34 - t9) * t40 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t3, t1, t2, 0, t49 + t79 * (-t37 * t12 - t76), (-t12 * t79 - t9) * t40 + t45;];
tauc_reg = t6;
