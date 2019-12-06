% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPPRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tau_reg [5x13]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:09
% EndTime: 2019-12-05 14:58:10
% DurationCPUTime: 0.30s
% Computational Cost: add. (285->65), mult. (640->119), div. (0->0), fcn. (589->12), ass. (0->56)
t37 = cos(pkin(8));
t25 = -t37 * qJDD(1) + qJDD(3);
t34 = sin(pkin(8));
t35 = sin(pkin(7));
t38 = cos(pkin(7));
t54 = (g(1) * t38 + g(2) * t35) * t34;
t82 = -g(3) * t37 - t25 + t54;
t33 = sin(pkin(9));
t36 = cos(pkin(9));
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t22 = t42 * t33 + t40 * t36;
t13 = t22 * t34;
t7 = t13 * qJD(4);
t44 = qJD(4) ^ 2;
t65 = qJDD(1) * t34;
t17 = t36 * qJDD(2) - t33 * t65;
t18 = t33 * qJDD(2) + t36 * t65;
t30 = pkin(9) + qJ(4);
t27 = sin(t30);
t28 = cos(t30);
t75 = t38 * t28;
t76 = t38 * t27;
t77 = t35 * t37;
t78 = g(3) * t34;
t81 = g(1) * (t35 * t27 + t37 * t75) + g(2) * (t28 * t77 - t76) + t28 * t78 - t40 * t17 - t42 * t18;
t80 = -g(1) * (t35 * t28 - t37 * t76) - g(2) * (-t27 * t77 - t75) + t27 * t78 + t42 * t17 - t40 * t18;
t21 = t40 * t33 - t42 * t36;
t69 = qJD(4) * t21;
t79 = t69 * qJD(4);
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t74 = t39 * t41;
t31 = t39 ^ 2;
t73 = -t41 ^ 2 + t31;
t72 = qJD(4) * pkin(4);
t14 = t21 * t34;
t68 = qJD(5) * t14;
t66 = qJDD(1) - g(3);
t64 = t14 * qJDD(5);
t63 = t41 * qJDD(4);
t62 = qJD(4) * qJD(5);
t61 = -g(1) * t35 + g(2) * t38;
t56 = -t13 * qJDD(4) + t34 * t79;
t55 = -t21 * qJDD(4) - t44 * t22;
t53 = qJD(5) * t37 + 0.2e1 * t7;
t43 = qJD(5) ^ 2;
t50 = t22 * t43 - t55;
t49 = t37 * qJDD(5) - t56;
t48 = 0.2e1 * t69 * qJD(5) - qJDD(5) * t22;
t47 = -pkin(6) * qJDD(5) - 0.2e1 * t72 * qJD(5);
t46 = -qJDD(4) * pkin(6) + qJD(4) * t72 + t81;
t45 = 0.2e1 * qJDD(4) * pkin(4) - pkin(6) * t43 + t80;
t24 = qJDD(5) * t41 - t43 * t39;
t23 = qJDD(5) * t39 + t43 * t41;
t1 = [t66, -g(3) + (t34 ^ 2 + t37 ^ 2) * qJDD(1), -t25 * t37 - g(3) + (-t17 * t33 + t18 * t36) * t34, 0, t56, t7 * qJD(4) + t14 * qJDD(4), 0, 0, 0, 0, 0, t39 * t64 - t49 * t41 + (t53 * t39 + t41 * t68) * qJD(5), t41 * t64 + t49 * t39 + (-t39 * t68 + t53 * t41) * qJD(5); 0, qJDD(2) + t61, t17 * t36 + t18 * t33 + t61, 0, t55, -t22 * qJDD(4) + t79, 0, 0, 0, 0, 0, t48 * t39 - t50 * t41, t50 * t39 + t48 * t41; 0, 0, -t66 * t37 + qJDD(3) - t54, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23; 0, 0, 0, qJDD(4), t80, t81, t31 * qJDD(4) + 0.2e1 * t62 * t74, 0.2e1 * t39 * t63 - 0.2e1 * t73 * t62, t23, t24, 0, t47 * t39 + t45 * t41, -t45 * t39 + t47 * t41; 0, 0, 0, 0, 0, 0, -t44 * t74, t73 * t44, t39 * qJDD(4), t63, qJDD(5), t46 * t39 - t82 * t41, t82 * t39 + t46 * t41;];
tau_reg = t1;
