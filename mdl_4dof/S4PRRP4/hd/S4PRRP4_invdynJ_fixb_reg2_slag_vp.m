% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:00
% EndTime: 2019-12-31 16:28:01
% DurationCPUTime: 0.40s
% Computational Cost: add. (362->111), mult. (745->142), div. (0->0), fcn. (393->4), ass. (0->68)
t39 = pkin(6) + qJ(2);
t34 = sin(t39);
t35 = cos(t39);
t58 = g(1) * t35 + g(2) * t34;
t44 = cos(qJ(3));
t43 = sin(qJ(3));
t76 = qJD(2) * t43;
t65 = pkin(5) * t76;
t14 = t44 * qJD(1) - t65;
t88 = qJD(4) - t14;
t8 = -qJD(3) * pkin(3) + t88;
t15 = t44 * qJD(2) * pkin(5) + t43 * qJD(1);
t10 = qJD(3) * qJ(4) + t15;
t57 = t44 * pkin(3) + t43 * qJ(4);
t16 = -pkin(2) - t57;
t56 = pkin(3) * t43 - qJ(4) * t44;
t7 = t56 * qJD(3) - t43 * qJD(4);
t1 = qJD(2) * t7 + qJDD(2) * t16;
t75 = qJDD(3) * pkin(3);
t87 = qJDD(4) - t75;
t77 = pkin(5) * qJDD(3);
t9 = qJD(2) * t16;
t86 = 0.2e1 * qJD(3) * t9 - t77;
t85 = g(1) * t34;
t82 = g(2) * t35;
t46 = qJD(2) ^ 2;
t81 = t43 * t46;
t80 = t35 * pkin(2) + t34 * pkin(5);
t40 = t43 ^ 2;
t41 = t44 ^ 2;
t79 = -t40 + t41;
t72 = t40 * qJDD(2);
t71 = t41 * qJDD(2);
t36 = t43 * qJDD(2);
t70 = t44 * qJDD(2);
t69 = qJD(1) * qJD(3);
t68 = qJD(2) * qJD(3);
t67 = qJDD(3) * qJ(4);
t66 = pkin(5) * t70 + t43 * qJDD(1) + t44 * t69;
t64 = t44 * t68;
t62 = -t44 * qJDD(1) + t43 * t69 + (t36 + t64) * pkin(5);
t60 = t43 * t64;
t45 = qJD(3) ^ 2;
t59 = pkin(5) * t45 + t82;
t55 = g(3) * t43 - t66;
t54 = -t58 + (t71 + t72) * pkin(5);
t53 = -0.2e1 * pkin(2) * t68 - t77;
t52 = 0.2e1 * qJDD(2) * pkin(2) - t59;
t51 = -g(3) * t44 + t58 * t43 - t62;
t50 = t15 * qJD(3) + t51;
t49 = -0.2e1 * t1 - t59;
t2 = t67 + (qJD(4) - t65) * qJD(3) + t66;
t3 = t62 + t87;
t48 = t2 * t44 + t3 * t43 + (-t10 * t43 + t44 * t8) * qJD(3);
t4 = -qJD(3) * t65 + t66;
t47 = t4 * t44 + t62 * t43 + (-t14 * t44 - t15 * t43) * qJD(3);
t42 = qJDD(1) - g(3);
t26 = t35 * pkin(5);
t23 = t44 * t81;
t21 = t44 * t85;
t19 = t79 * t46;
t18 = qJDD(3) * t44 - t45 * t43;
t17 = qJDD(3) * t43 + t45 * t44;
t13 = t56 * qJD(2);
t12 = -0.2e1 * t60 + t71;
t11 = 0.2e1 * t60 + t72;
t6 = t43 * t70 + t79 * t68;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t4 * t43 - t62 * t44 - g(3) + (-t14 * t43 + t15 * t44) * qJD(3), 0, 0, 0, 0, 0, 0, t18, 0, t17, t2 * t43 - t3 * t44 - g(3) + (t10 * t44 + t43 * t8) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t82 + t85, t58, 0, 0, t11, 0.2e1 * t6, t17, t12, t18, 0, t53 * t43 + t52 * t44 + t21, t53 * t44 + (-t52 - t85) * t43, t47 + t54, qJDD(2) * pkin(2) ^ 2 - g(1) * (-t34 * pkin(2) + t26) - g(2) * t80 + t47 * pkin(5), t11, t17, -0.2e1 * t6, 0, -t18, t12, t86 * t43 + t49 * t44 + t21, t48 + t54, -t86 * t44 + (t49 + t85) * t43, t9 * t7 - g(1) * t26 - g(2) * (t57 * t35 + t80) + t48 * pkin(5) + (t1 - t85) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t19, t36, t23, t70, qJDD(3), pkin(2) * t81 + t50, (t14 + t65) * qJD(3) + (pkin(2) * t46 + t58) * t44 + t55, 0, 0, -t23, t36, t19, qJDD(3), -t70, t23, 0.2e1 * t75 - qJDD(4) + (t13 * t44 - t43 * t9) * qJD(2) + t50, -t56 * qJDD(2), 0.2e1 * t67 - t58 * t44 + (0.2e1 * qJD(4) - t14) * qJD(3) + (t44 * t9 + (-pkin(5) * qJD(3) + t13) * t43) * qJD(2) - t55, -t3 * pkin(3) - g(3) * t57 + t2 * qJ(4) + t88 * t10 - t9 * t13 - t8 * t15 + t58 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t23, t36, -t40 * t46 - t45, -t10 * qJD(3) + t9 * t76 - t51 + t87;];
tau_reg = t5;
