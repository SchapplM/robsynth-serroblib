% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:59
% EndTime: 2019-12-31 16:28:00
% DurationCPUTime: 0.31s
% Computational Cost: add. (294->89), mult. (565->118), div. (0->0), fcn. (296->4), ass. (0->56)
t27 = pkin(6) + qJ(2);
t22 = sin(t27);
t23 = cos(t27);
t44 = g(1) * t23 + g(2) * t22;
t31 = cos(qJ(3));
t30 = sin(qJ(3));
t59 = qJD(2) * t30;
t50 = pkin(5) * t59;
t9 = t31 * qJD(1) - t50;
t71 = qJD(4) - t9;
t43 = t31 * pkin(3) + t30 * qJ(4);
t40 = pkin(2) + t43;
t6 = t40 * qJD(2);
t42 = pkin(3) * t30 - qJ(4) * t31;
t4 = t42 * qJD(3) - t30 * qJD(4);
t1 = qJD(2) * t4 - t40 * qJDD(2);
t5 = -qJD(3) * pkin(3) + t71;
t10 = t31 * qJD(2) * pkin(5) + t30 * qJD(1);
t7 = qJD(3) * qJ(4) + t10;
t64 = g(2) * t23;
t67 = g(1) * t22;
t70 = -t64 + t67;
t58 = qJDD(3) * pkin(3);
t69 = qJDD(4) - t58;
t60 = pkin(5) * qJDD(3);
t68 = -0.2e1 * t6 * qJD(3) - t60;
t33 = qJD(2) ^ 2;
t63 = t30 * t33;
t28 = t30 ^ 2;
t29 = t31 ^ 2;
t62 = t28 - t29;
t24 = t30 * qJDD(2);
t56 = t31 * qJDD(2);
t55 = qJD(1) * qJD(3);
t54 = qJD(2) * qJD(3);
t53 = qJDD(3) * qJ(4);
t52 = t31 * t63;
t51 = pkin(5) * t56 + t30 * qJDD(1) + t31 * t55;
t49 = t31 * t54;
t47 = -t31 * qJDD(1) + t30 * t55 + (t24 + t49) * pkin(5);
t32 = qJD(3) ^ 2;
t45 = pkin(5) * t32 + t64;
t41 = g(3) * t30 - t51;
t39 = -0.2e1 * pkin(2) * t54 - t60;
t38 = 0.2e1 * qJDD(2) * pkin(2) - t45;
t37 = -g(3) * t31 + t44 * t30 - t47;
t36 = t10 * qJD(3) + t37;
t35 = -0.2e1 * t1 - t45;
t2 = t53 + (qJD(4) - t50) * qJD(3) + t51;
t3 = t47 + t69;
t34 = t2 * t31 + t3 * t30 + (-t30 * t7 + t31 * t5) * qJD(3) - t44;
t15 = t31 * t67;
t13 = qJDD(3) * t31 - t32 * t30;
t12 = qJDD(3) * t30 + t32 * t31;
t8 = t42 * qJD(2);
t11 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t13, 0, t12, t2 * t30 - t3 * t31 - g(3) + (t30 * t5 + t31 * t7) * qJD(3); 0, qJDD(2), t70, t44, t28 * qJDD(2) + 0.2e1 * t30 * t49, 0.2e1 * t30 * t56 - 0.2e1 * t62 * t54, t12, t13, 0, t39 * t30 + t38 * t31 + t15, t39 * t31 + (-t38 - t67) * t30, t68 * t30 + t35 * t31 + t15, (t28 + t29) * qJDD(2) * pkin(5) + t34, -t68 * t31 + (t35 + t67) * t30, t34 * pkin(5) - t6 * t4 + (-t1 + t70) * t40; 0, 0, 0, 0, -t52, t62 * t33, t24, t56, qJDD(3), pkin(2) * t63 + t36, (t9 + t50) * qJD(3) + (pkin(2) * t33 + t44) * t31 + t41, 0.2e1 * t58 - qJDD(4) + (t30 * t6 + t31 * t8) * qJD(2) + t36, -t42 * qJDD(2), 0.2e1 * t53 - t44 * t31 + (0.2e1 * qJD(4) - t9) * qJD(3) + (-t31 * t6 + (-pkin(5) * qJD(3) + t8) * t30) * qJD(2) - t41, -t3 * pkin(3) - g(3) * t43 + t2 * qJ(4) - t5 * t10 + t44 * t42 + t6 * t8 + t71 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t52, t24, -t28 * t33 - t32, -t7 * qJD(3) - t6 * t59 - t37 + t69;];
tau_reg = t11;
