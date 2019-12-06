% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:57
% EndTime: 2019-12-05 15:26:59
% DurationCPUTime: 0.41s
% Computational Cost: add. (345->102), mult. (663->130), div. (0->0), fcn. (493->10), ass. (0->69)
t33 = qJ(2) + pkin(8);
t30 = sin(t33);
t31 = cos(t33);
t43 = cos(qJ(2));
t32 = t43 * qJDD(1);
t41 = sin(qJ(2));
t64 = qJD(1) * qJD(2);
t15 = qJDD(2) * pkin(2) - t41 * t64 + t32;
t36 = sin(pkin(8));
t38 = cos(pkin(8));
t66 = t41 * qJDD(1);
t89 = t43 * t64 + t66;
t4 = t38 * t15 - t89 * t36;
t55 = qJDD(4) - t4;
t37 = sin(pkin(7));
t39 = cos(pkin(7));
t58 = g(1) * t39 + g(2) * t37;
t47 = g(3) * t31 - t58 * t30 + t55;
t73 = t43 * qJD(1);
t22 = qJD(2) * pkin(2) + t73;
t75 = qJD(1) * t41;
t10 = t36 * t22 + t38 * t75;
t8 = qJD(2) * qJ(4) + t10;
t91 = -t8 * qJD(2) + t47;
t90 = -(-pkin(3) - pkin(6)) * qJDD(2) - t91;
t87 = -t36 * t41 + t38 * t43;
t17 = t36 * t43 + t38 * t41;
t12 = t17 * qJD(1);
t29 = -t38 * pkin(2) - pkin(3);
t26 = -pkin(6) + t29;
t27 = t36 * pkin(2) + qJ(4);
t86 = qJDD(5) * t26 + (qJD(2) * t27 - t12 + t8) * qJD(5);
t11 = qJD(2) * t17;
t85 = 0.2e1 * qJD(5) * t11 - qJDD(5) * t87;
t82 = t36 * t15;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t79 = t40 * t42;
t35 = t42 ^ 2;
t78 = t40 ^ 2 - t35;
t44 = qJD(5) ^ 2;
t45 = qJD(2) ^ 2;
t77 = -t44 - t45;
t74 = qJDD(2) * pkin(3);
t24 = t36 * t75;
t62 = t38 * t73;
t14 = -t24 + t62;
t72 = qJD(4) - t14;
t71 = qJDD(1) - g(3);
t68 = qJDD(5) * t40;
t67 = qJDD(5) * t42;
t65 = t42 * qJDD(2);
t63 = qJD(2) * qJD(5);
t9 = t38 * t22 - t24;
t13 = t87 * qJD(2);
t57 = t13 * qJD(2) + t17 * qJDD(2);
t54 = g(1) * t37 - g(2) * t39 - qJDD(3);
t52 = t38 * t66 + t82;
t50 = t44 * t87 + t57;
t49 = -g(3) * t30 - t58 * t31;
t48 = -g(3) * t43 + t58 * t41;
t2 = qJDD(2) * qJ(4) + (qJD(4) + t62) * qJD(2) + t52;
t46 = t72 * qJD(2) + t27 * qJDD(2) - t26 * t44 + t2 + t49;
t20 = -t44 * t40 + t67;
t19 = -t44 * t42 - t68;
t7 = -qJD(2) * pkin(3) + qJD(4) - t9;
t5 = t89 * t38 + t82;
t3 = t55 - t74;
t1 = [t71, 0, t43 * qJDD(2) - t45 * t41, -qJDD(2) * t41 - t45 * t43, t10 * t13 - t9 * t11 + t5 * t17 + t4 * t87 - g(3), t11 * qJD(2) - qJDD(2) * t87, t57, t7 * t11 + t8 * t13 + t2 * t17 - t3 * t87 - g(3), 0, 0, 0, 0, 0, t50 * t40 + t85 * t42, -t85 * t40 + t50 * t42; 0, qJDD(2), t32 + t48, -t71 * t41 + t58 * t43, -t10 * t14 + t9 * t12 + (t36 * t5 + t38 * t4 + t48) * pkin(2), -t12 * qJD(2) + (-pkin(3) + t29) * qJDD(2) + t47, (qJ(4) + t27) * qJDD(2) + (0.2e1 * qJD(4) - t14 + t62) * qJD(2) + t49 + t52, t2 * t27 + t3 * t29 - t7 * t12 - g(3) * (t43 * pkin(2) + t31 * pkin(3) + t30 * qJ(4)) + t72 * t8 + t58 * (pkin(2) * t41 + pkin(3) * t30 - qJ(4) * t31), t35 * qJDD(2) - 0.2e1 * t63 * t79, -0.2e1 * t40 * t65 + 0.2e1 * t78 * t63, t20, t19, 0, t46 * t40 + t86 * t42, -t86 * t40 + t46 * t42; 0, 0, 0, 0, -t54, 0, 0, -t54, 0, 0, 0, 0, 0, t19, -t20; 0, 0, 0, 0, 0, qJDD(2), -t45, -t74 + t91, 0, 0, 0, 0, 0, t77 * t40 + t67, t77 * t42 - t68; 0, 0, 0, 0, 0, 0, 0, 0, t45 * t79, -t78 * t45, t65, -t40 * qJDD(2), qJDD(5), t54 * t40 - t90 * t42, t90 * t40 + t54 * t42;];
tau_reg = t1;
