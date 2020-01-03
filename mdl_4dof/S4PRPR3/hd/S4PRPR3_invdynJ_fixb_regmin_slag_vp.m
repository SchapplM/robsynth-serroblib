% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR3
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:56
% EndTime: 2019-12-31 16:20:57
% DurationCPUTime: 0.26s
% Computational Cost: add. (247->75), mult. (504->101), div. (0->0), fcn. (378->8), ass. (0->52)
t61 = qJD(2) * qJD(3);
t71 = qJ(3) * qJDD(2) + t61;
t45 = sin(pkin(7));
t46 = cos(pkin(7));
t67 = t45 ^ 2 + t46 ^ 2;
t47 = sin(qJ(4));
t48 = cos(qJ(4));
t16 = t48 * t45 + t47 * t46;
t12 = t16 * qJD(2);
t70 = t47 * t45;
t69 = t48 * t46;
t68 = pkin(5) + qJ(3);
t66 = qJDD(2) * pkin(2);
t64 = t45 * qJDD(2);
t63 = t46 * qJDD(2);
t60 = qJD(4) * qJD(2) * t69 + t47 * t63 + t48 * t64;
t8 = t45 * qJDD(1) + t71 * t46;
t59 = qJD(2) * t70;
t31 = -t46 * pkin(3) - pkin(2);
t44 = pkin(6) + qJ(2);
t38 = sin(t44);
t40 = cos(t44);
t58 = g(1) * t40 + g(2) * t38;
t57 = g(1) * t38 - g(2) * t40;
t56 = t47 * t64 - t48 * t63;
t55 = t67 * qJ(3) * qJD(2);
t21 = t68 * t45;
t22 = t68 * t46;
t54 = -t48 * t21 - t47 * t22;
t53 = -t47 * t21 + t48 * t22;
t15 = -t69 + t70;
t14 = t16 * qJD(4);
t52 = -t57 - t66;
t32 = qJDD(3) - t66;
t51 = -t32 - t52;
t34 = t46 * qJDD(1);
t7 = -t71 * t45 + t34;
t50 = -t7 * t45 + t8 * t46 - t58;
t43 = pkin(7) + qJ(4);
t39 = cos(t43);
t37 = sin(t43);
t20 = t31 * qJD(2) + qJD(3);
t19 = t31 * qJDD(2) + qJDD(3);
t13 = t15 * qJD(4);
t11 = t15 * qJD(2);
t6 = pkin(5) * t63 + t8;
t5 = t34 + (-t68 * qJDD(2) - t61) * t45;
t4 = -t14 * qJD(4) - t15 * qJDD(4);
t3 = -t13 * qJD(4) + t16 * qJDD(4);
t2 = qJD(2) * t14 + t56;
t1 = -qJD(4) * t59 + t60;
t9 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, t8 * t45 + t7 * t46 - g(3), 0, 0, 0, 0, 0, t4, -t3; 0, qJDD(2), t57, t58, t51 * t46, -t51 * t45, t71 * t67 + t50, t55 * qJD(3) + (-t32 + t57) * pkin(2) + t50 * qJ(3), t1 * t16 - t12 * t13, -t1 * t15 + t13 * t11 - t12 * t14 - t16 * t2, t3, t4, 0, t31 * t2 + t19 * t15 + t20 * t14 + t54 * qJDD(4) + t57 * t39 + (-t16 * qJD(3) - t53 * qJD(4)) * qJD(4), t31 * t1 + t19 * t16 - t20 * t13 - t53 * qJDD(4) - t57 * t37 + (t15 * qJD(3) - t54 * qJD(4)) * qJD(4); 0, 0, 0, 0, -t63, t64, -t67 * qJD(2) ^ 2, -t55 * qJD(2) + qJDD(3) + t52, 0, 0, 0, 0, 0, 0.2e1 * t12 * qJD(4) + t56, (-t11 - t59) * qJD(4) + t60; 0, 0, 0, 0, 0, 0, 0, 0, t12 * t11, -t11 ^ 2 + t12 ^ 2, (t11 - t59) * qJD(4) + t60, -t56, qJDD(4), -g(3) * t39 - t20 * t12 + t58 * t37 - t47 * t6 + t48 * t5, g(3) * t37 + t20 * t11 + t58 * t39 - t47 * t5 - t48 * t6;];
tau_reg = t9;
