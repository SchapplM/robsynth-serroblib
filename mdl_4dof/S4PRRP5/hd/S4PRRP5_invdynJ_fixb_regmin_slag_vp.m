% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP5
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
% tau_reg [4x13]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:23
% EndTime: 2019-12-31 16:29:24
% DurationCPUTime: 0.37s
% Computational Cost: add. (275->98), mult. (631->141), div. (0->0), fcn. (378->6), ass. (0->61)
t64 = qJ(4) + pkin(5);
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t21 = sin(pkin(6));
t22 = cos(pkin(6));
t43 = g(1) * t22 + g(2) * t21;
t38 = t43 * t27;
t32 = g(3) * t25 + t38;
t29 = qJD(2) ^ 2;
t73 = qJDD(2) * t25 + t29 * t27;
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t44 = t25 * qJD(1) + t64 * qJD(2);
t5 = t44 * t24;
t59 = qJD(3) * pkin(3);
t4 = -t5 + t59;
t6 = t44 * t26;
t41 = t24 * t4 - t26 * t6;
t72 = t41 * qJD(2);
t48 = qJD(1) * qJD(2);
t17 = t25 * t48;
t28 = qJD(3) ^ 2;
t50 = t27 * qJDD(1);
t71 = 0.2e1 * qJDD(2) * pkin(2) - pkin(5) * t28 - g(3) * t27 + t25 * (t43 + t48) - t17 + t50;
t55 = qJDD(1) - g(3);
t70 = t43 * t25 + t55 * t27;
t69 = t5 + t4;
t19 = t24 ^ 2;
t20 = t26 ^ 2;
t63 = t19 - t20;
t62 = t19 + t20;
t61 = t28 + t29;
t60 = qJD(2) * pkin(2);
t18 = t26 * pkin(3) + pkin(2);
t56 = t27 * qJD(1);
t9 = -t18 * qJD(2) + qJD(4) - t56;
t58 = t9 * qJD(2);
t53 = qJDD(3) * t24;
t52 = t24 * qJDD(2);
t51 = t26 * qJDD(2);
t49 = t27 * qJDD(2);
t47 = qJD(2) * qJD(3);
t46 = t24 * t47;
t45 = qJD(3) * t64;
t42 = g(1) * t21 - g(2) * t22;
t39 = qJD(3) * t44;
t37 = t42 * t26;
t11 = qJDD(2) * pkin(5) + t25 * qJDD(1) + t27 * t48;
t36 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t11;
t33 = pkin(3) * t46 - t18 * qJDD(2) + qJDD(4) + t17;
t15 = -t56 - t60;
t31 = -pkin(5) * qJDD(3) + (t15 + t56 - t60) * qJD(3);
t30 = -t15 * qJD(2) - t11 + t32;
t13 = t64 * t26;
t12 = t64 * t24;
t8 = -t24 * qJD(4) - t26 * t45;
t7 = t26 * qJD(4) - t24 * t45;
t3 = t33 - t50;
t2 = -t24 * t39 + t36 * t26;
t1 = qJDD(3) * pkin(3) - t36 * t24 - t26 * t39;
t10 = [t55, 0, -t29 * t25 + t49, -t73, 0, 0, 0, 0, 0, (-0.2e1 * t46 + t51) * t27 + (-t61 * t26 - t53) * t25, (-qJDD(3) * t25 - 0.2e1 * t27 * t47) * t26 + (t61 * t25 - t49) * t24, t73 * t62, -g(3) + (-t3 - t72) * t27 + (t58 - t1 * t24 + t2 * t26 + (-t24 * t6 - t26 * t4) * qJD(3)) * t25; 0, qJDD(2), t70, -t55 * t25 + t38, t19 * qJDD(2) + 0.2e1 * t26 * t46, 0.2e1 * t24 * t51 - 0.2e1 * t63 * t47, t28 * t26 + t53, qJDD(3) * t26 - t28 * t24, 0, t31 * t24 + t71 * t26, -t71 * t24 + t31 * t26, (-qJD(3) * t4 + qJDD(2) * t13 + t2) * t26 + (-qJD(3) * t6 + qJDD(2) * t12 - t1) * t24 + (-t24 * t8 + t26 * t7 + (t12 * t26 - t13 * t24) * qJD(3) - t62 * t56) * qJD(2) - t32, t2 * t13 + t6 * t7 - t1 * t12 + t4 * t8 - t3 * t18 + t9 * t24 * t59 - g(3) * (t27 * t18 + t25 * t64) + (-t9 * t25 + t41 * t27) * qJD(1) + t43 * (t18 * t25 - t27 * t64); 0, 0, 0, 0, -t24 * t29 * t26, t63 * t29, t52, t51, qJDD(3), t30 * t24 - t37, t42 * t24 + t30 * t26, -pkin(3) * t52 + (-t59 + t69) * t26 * qJD(2), t69 * t6 + (t1 - t37 + (t32 - t58) * t24) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t29, t33 - t70 + t72;];
tau_reg = t10;
