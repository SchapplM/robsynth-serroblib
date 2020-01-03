% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:57
% EndTime: 2019-12-31 16:37:58
% DurationCPUTime: 0.28s
% Computational Cost: add. (297->85), mult. (589->120), div. (0->0), fcn. (423->12), ass. (0->58)
t48 = sin(pkin(7));
t50 = cos(pkin(7));
t73 = t48 ^ 2 + t50 ^ 2;
t51 = cos(pkin(6));
t35 = -pkin(1) * t51 - pkin(2);
t72 = qJDD(1) * t35;
t26 = qJDD(3) + t72;
t47 = qJ(1) + pkin(6);
t41 = sin(t47);
t43 = cos(t47);
t65 = g(1) * t41 - g(2) * t43;
t77 = -t26 + t65;
t52 = sin(qJ(4));
t54 = cos(qJ(4));
t25 = t48 * t54 + t50 * t52;
t17 = t25 * qJD(1);
t49 = sin(pkin(6));
t33 = pkin(1) * t49 + qJ(3);
t22 = qJD(1) * qJD(3) + qJDD(1) * t33;
t76 = pkin(5) + t33;
t75 = t52 * t48;
t74 = t54 * t50;
t10 = t48 * qJDD(2) + t50 * t22;
t71 = t48 * qJDD(1);
t70 = t50 * qJDD(1);
t68 = qJD(4) * qJD(1) * t74 + t52 * t70 + t54 * t71;
t67 = qJD(1) * t75;
t66 = g(1) * t43 + g(2) * t41;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t64 = g(1) * t53 - g(2) * t55;
t63 = t52 * t71 - t54 * t70;
t37 = t50 * qJDD(2);
t9 = -t22 * t48 + t37;
t62 = t10 * t50 - t9 * t48;
t61 = t73 * t33 * qJD(1);
t20 = t76 * t48;
t21 = t76 * t50;
t60 = -t20 * t54 - t21 * t52;
t59 = -t20 * t52 + t21 * t54;
t24 = -t74 + t75;
t27 = -pkin(3) * t50 + t35;
t19 = t25 * qJD(4);
t58 = -t72 + t77;
t46 = pkin(7) + qJ(4);
t42 = cos(t46);
t40 = sin(t46);
t18 = t24 * qJD(4);
t16 = t24 * qJD(1);
t15 = qJD(1) * t27 + qJD(3);
t13 = qJDD(1) * t27 + qJDD(3);
t6 = pkin(5) * t70 + t10;
t5 = t37 + (-pkin(5) * qJDD(1) - t22) * t48;
t4 = -qJD(4) * t19 - qJDD(4) * t24;
t3 = -qJD(4) * t18 + qJDD(4) * t25;
t2 = qJD(1) * t19 + t63;
t1 = -qJD(4) * t67 + t68;
t7 = [qJDD(1), t64, g(1) * t55 + g(2) * t53, (t64 + (t49 ^ 2 + t51 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t58 * t50, -t58 * t48, t22 * t73 + t62 - t66, t26 * t35 - g(1) * (-pkin(1) * t53 - pkin(2) * t41 + qJ(3) * t43) - g(2) * (pkin(1) * t55 + pkin(2) * t43 + qJ(3) * t41) + t62 * t33 + t61 * qJD(3), t1 * t25 - t17 * t18, -t1 * t24 + t16 * t18 - t17 * t19 - t2 * t25, t3, t4, 0, t27 * t2 + t13 * t24 + t15 * t19 + t60 * qJDD(4) + t65 * t42 + (-qJD(3) * t25 - qJD(4) * t59) * qJD(4), t27 * t1 + t13 * t25 - t15 * t18 - t59 * qJDD(4) - t65 * t40 + (qJD(3) * t24 - qJD(4) * t60) * qJD(4); 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t10 * t48 + t50 * t9 - g(3), 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, -t70, t71, -t73 * qJD(1) ^ 2, -qJD(1) * t61 - t77, 0, 0, 0, 0, 0, 0.2e1 * t17 * qJD(4) + t63, (-t16 - t67) * qJD(4) + t68; 0, 0, 0, 0, 0, 0, 0, 0, t17 * t16, -t16 ^ 2 + t17 ^ 2, (t16 - t67) * qJD(4) + t68, -t63, qJDD(4), -g(3) * t42 - t15 * t17 + t40 * t66 + t54 * t5 - t52 * t6, g(3) * t40 + t15 * t16 + t42 * t66 - t52 * t5 - t54 * t6;];
tau_reg = t7;
