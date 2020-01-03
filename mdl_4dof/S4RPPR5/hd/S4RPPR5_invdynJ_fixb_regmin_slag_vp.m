% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tau_reg [4x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:50
% DurationCPUTime: 0.23s
% Computational Cost: add. (231->86), mult. (405->110), div. (0->0), fcn. (243->6), ass. (0->56)
t41 = -pkin(1) - pkin(2);
t20 = t41 * qJD(1) + qJD(2);
t36 = cos(pkin(6));
t73 = t36 * t20;
t19 = t41 * qJDD(1) + qJDD(2);
t35 = sin(pkin(6));
t59 = qJ(2) * qJDD(1);
t72 = t35 * t19 + t36 * t59;
t15 = t36 * qJ(2) + t35 * t41;
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t71 = t40 * pkin(1) + t38 * qJ(2);
t70 = g(1) * t38 - g(2) * t40;
t37 = sin(qJ(4));
t33 = t37 ^ 2;
t39 = cos(qJ(4));
t69 = -t39 ^ 2 + t33;
t42 = qJD(4) ^ 2;
t43 = qJD(1) ^ 2;
t68 = t42 + t43;
t67 = t35 * qJ(2);
t66 = pkin(1) * qJDD(1);
t65 = qJDD(3) + g(3);
t64 = qJ(2) * qJD(1);
t63 = qJDD(4) * t37;
t62 = qJDD(4) * t39;
t61 = t37 * qJDD(1);
t60 = t39 * qJDD(1);
t58 = qJD(1) * qJD(2);
t57 = qJD(1) * qJD(4);
t23 = t36 * t58;
t4 = t23 + t72;
t56 = 0.2e1 * t58;
t55 = 0.2e1 * t57;
t54 = t35 * t58;
t53 = t36 * t19 - t35 * t59;
t52 = t39 * t55;
t51 = qJDD(2) - t66;
t10 = t38 * t35 + t40 * t36;
t9 = t40 * t35 - t38 * t36;
t50 = g(1) * t9 + g(2) * t10;
t49 = g(1) * t10 - g(2) * t9;
t48 = g(1) * t40 + g(2) * t38;
t47 = (-t35 * t64 + t73) * t35 - (t35 * t20 + t36 * t64) * t36;
t14 = t36 * t41 - t67;
t3 = t53 - t54;
t5 = -t73 + (pkin(3) + t67) * qJD(1);
t46 = qJDD(1) * pkin(5) + t5 * qJD(1) - t4 + t49;
t11 = pkin(3) - t14;
t12 = -pkin(5) + t15;
t45 = -qJDD(4) * t12 + (-qJD(1) * t11 - qJD(2) * t36 - t5) * qJD(4);
t44 = -t12 * t42 - t3 - t50 + t54 + (pkin(3) + t11) * qJDD(1);
t29 = t40 * qJ(2);
t18 = -t42 * t37 + t62;
t17 = -t42 * t39 - t63;
t1 = [qJDD(1), t70, t48, -qJDD(2) + 0.2e1 * t66 + t70, -t48 + t56 + 0.2e1 * t59, -t51 * pkin(1) - g(1) * (-t38 * pkin(1) + t29) - g(2) * t71 + (t56 + t59) * qJ(2), -t14 * qJDD(1) - t50 - t53 + 0.2e1 * t54, t15 * qJDD(1) + 0.2e1 * t23 - t49 + t72, t4 * t15 + t3 * t14 - g(1) * (t41 * t38 + t29) - g(2) * (t40 * pkin(2) + t71) - t47 * qJD(2), t33 * qJDD(1) + t37 * t52, 0.2e1 * t37 * t60 - 0.2e1 * t69 * t57, t17, -t18, 0, t45 * t37 + t44 * t39, -t44 * t37 + t45 * t39; 0, 0, 0, -qJDD(1), -t43, -t43 * qJ(2) + t51 - t70, -t36 * qJDD(1) - t35 * t43, t35 * qJDD(1) - t36 * t43, t47 * qJD(1) + t3 * t36 + t4 * t35 - t70, 0, 0, 0, 0, 0, (t37 * t55 - t60) * t36 + (-t68 * t39 - t63) * t35, (t52 + t61) * t36 + (t68 * t37 - t62) * t35; 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t43 * t39, t69 * t43, -t61, -t60, qJDD(4), t46 * t37 + t65 * t39, -t65 * t37 + t46 * t39;];
tau_reg = t1;
