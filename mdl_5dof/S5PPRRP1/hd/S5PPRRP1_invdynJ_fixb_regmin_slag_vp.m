% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:16
% EndTime: 2019-12-05 15:07:18
% DurationCPUTime: 0.53s
% Computational Cost: add. (517->109), mult. (1135->152), div. (0->0), fcn. (830->10), ass. (0->70)
t38 = sin(pkin(8));
t40 = cos(pkin(8));
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t99 = -t44 * t38 + t46 * t40;
t14 = t99 * qJD(1);
t19 = t46 * t38 + t44 * t40;
t15 = t19 * qJD(1);
t35 = pkin(8) + qJ(3);
t30 = sin(t35);
t31 = cos(t35);
t39 = sin(pkin(7));
t41 = cos(pkin(7));
t64 = g(1) * t41 + g(2) * t39;
t54 = g(3) * t31 - t64 * t30;
t17 = t19 * qJD(3);
t57 = -qJD(1) * t17 + t99 * qJDD(1);
t101 = t15 * qJD(3) - t54 + t57;
t82 = qJ(5) + pkin(6);
t47 = qJD(4) ^ 2;
t100 = (2 * qJDD(3) * pkin(3)) - pkin(6) * t47 + t101;
t97 = -g(1) * t39 + g(2) * t41;
t98 = qJDD(2) + t97;
t77 = qJD(3) * t99;
t96 = qJD(3) * t77 + t19 * qJDD(3);
t95 = t19 * qJDD(1);
t55 = g(3) * t30 + t64 * t31;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t65 = t82 * qJD(3) + t15;
t7 = t45 * qJD(2) - t65 * t43;
t78 = qJD(4) * pkin(4);
t6 = t7 + t78;
t93 = -t7 + t6;
t36 = t43 ^ 2;
t37 = t45 ^ 2;
t81 = t36 - t37;
t80 = t36 + t37;
t79 = qJD(3) * pkin(3);
t72 = t43 * qJDD(3);
t71 = t45 * qJDD(3);
t69 = qJD(3) * qJD(4);
t29 = t45 * pkin(4) + pkin(3);
t67 = t43 * t69;
t66 = qJD(4) * t82;
t8 = t43 * qJD(2) + t65 * t45;
t62 = t6 * t43 - t8 * t45;
t60 = t65 * qJD(4);
t59 = -t17 * qJD(3) + qJDD(3) * t99;
t58 = t97 * t45;
t56 = t19 * t47 - t59;
t53 = -0.2e1 * t77 * qJD(4) - qJDD(4) * t19;
t10 = -t14 - t79;
t52 = -pkin(6) * qJDD(4) + (t10 + t14 - t79) * qJD(4);
t4 = qJDD(3) * pkin(6) + qJD(1) * t77 + t95;
t51 = (qJ(5) * qJDD(3)) + qJD(4) * qJD(2) + qJD(3) * qJD(5) + t4;
t49 = -t10 * qJD(3) - t4 + t55;
t3 = pkin(4) * t67 - t29 * qJDD(3) + qJDD(5) - t57;
t48 = qJD(3) ^ 2;
t32 = t45 * qJDD(2);
t23 = t82 * t45;
t22 = t82 * t43;
t21 = qJDD(4) * t45 - t47 * t43;
t20 = qJDD(4) * t43 + t47 * t45;
t13 = -t43 * qJD(5) - t45 * t66;
t12 = t45 * qJD(5) - t43 * t66;
t9 = -t29 * qJD(3) + qJD(5) - t14;
t2 = (qJDD(2) - t60) * t43 + t51 * t45;
t1 = qJDD(4) * pkin(4) - t51 * t43 - t45 * t60 + t32;
t5 = [qJDD(1) - g(3), -g(3) + (t38 ^ 2 + t40 ^ 2) * qJDD(1), 0, t59, -t96, 0, 0, 0, 0, 0, t53 * t43 - t56 * t45, t56 * t43 + t53 * t45, t96 * t80, t9 * t17 - t3 * t99 - g(3) - t62 * t77 + (-t1 * t43 + t2 * t45 + (-t43 * t8 - t45 * t6) * qJD(4)) * t19; 0, t98, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t62 * qJD(4) + t1 * t45 + t2 * t43 + t97; 0, 0, qJDD(3), t101, -t95 + t55, t36 * qJDD(3) + 0.2e1 * t45 * t67, 0.2e1 * t43 * t71 - 0.2e1 * t81 * t69, t20, t21, 0, t100 * t45 + t52 * t43, -t100 * t43 + t52 * t45, (-qJD(4) * t6 + qJDD(3) * t23 + t2) * t45 + (-qJD(4) * t8 + qJDD(3) * t22 - t1) * t43 + (t12 * t45 - t13 * t43 - t80 * t14 + (t22 * t45 - t23 * t43) * qJD(4)) * qJD(3) - t55, t2 * t23 + t8 * t12 - t1 * t22 + t6 * t13 - t3 * t29 - g(3) * (t31 * t29 + t30 * t82) + (t43 * t78 - t15) * t9 + t62 * t14 + t64 * (t29 * t30 - t31 * t82); 0, 0, 0, 0, 0, -t43 * t48 * t45, t81 * t48, t72, t71, qJDD(4), t49 * t43 + t32 + t58, -t98 * t43 + t49 * t45, -pkin(4) * t72 + (-t78 + t93) * t45 * qJD(3), t93 * t8 + (t1 + t58 + (-t9 * qJD(3) + t55) * t43) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t48, t62 * qJD(3) + t3 + t54;];
tau_reg = t5;
