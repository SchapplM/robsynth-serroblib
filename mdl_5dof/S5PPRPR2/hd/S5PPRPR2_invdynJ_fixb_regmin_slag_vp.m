% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:24
% EndTime: 2019-12-05 15:03:25
% DurationCPUTime: 0.34s
% Computational Cost: add. (285->82), mult. (576->103), div. (0->0), fcn. (448->10), ass. (0->65)
t27 = pkin(8) + qJ(3);
t25 = sin(t27);
t26 = cos(t27);
t33 = sin(pkin(7));
t35 = cos(pkin(7));
t56 = g(1) * t35 + g(2) * t33;
t95 = g(3) * t26 - t56 * t25;
t32 = sin(pkin(8));
t34 = cos(pkin(8));
t37 = sin(qJ(3));
t65 = qJD(1) * qJD(3);
t61 = t37 * t65;
t39 = cos(qJ(3));
t71 = qJDD(1) * t39;
t92 = qJDD(1) * t37 + t39 * t65;
t48 = (-t61 + t71) * t34 - t92 * t32;
t46 = qJDD(4) - t48;
t44 = t46 + t95;
t73 = qJD(3) * qJ(4);
t12 = t39 * t32 + t37 * t34;
t8 = t12 * qJD(1);
t6 = t8 + t73;
t94 = -t6 * qJD(3) + t44;
t45 = -g(3) * t25 - t56 * t26;
t81 = t37 * t32;
t62 = qJD(1) * t81;
t63 = t32 * t71 + t92 * t34;
t80 = t39 * t34;
t7 = -qJD(1) * t80 + t62;
t93 = -t45 + (-t7 + t62) * qJD(3) - t63;
t40 = -pkin(3) - pkin(6);
t91 = -t40 * qJDD(3) - t94;
t11 = -t80 + t81;
t87 = (t6 - t8 + t73) * qJD(5) + qJDD(5) * t40;
t10 = qJD(3) * t12;
t86 = 0.2e1 * qJD(5) * t10 + qJDD(5) * t11;
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t82 = t36 * t38;
t31 = t38 ^ 2;
t79 = t36 ^ 2 - t31;
t41 = qJD(5) ^ 2;
t42 = qJD(3) ^ 2;
t78 = -t41 - t42;
t76 = t8 * qJD(3);
t75 = qJD(4) + t7;
t74 = qJDD(3) * pkin(3);
t69 = qJDD(5) * t36;
t68 = qJDD(5) * t38;
t66 = t38 * qJDD(3);
t64 = qJD(3) * qJD(5);
t28 = qJDD(3) * qJ(4);
t9 = t11 * qJD(3);
t54 = t9 * qJD(3) - t12 * qJDD(3);
t53 = t10 * qJD(3) + t11 * qJDD(3);
t51 = g(1) * t33 - g(2) * t35 - qJDD(2);
t47 = -t11 * t41 - t54;
t29 = qJD(3) * qJD(4);
t2 = t32 * t61 - t28 - t29 - t63;
t43 = t75 * qJD(3) - t40 * t41 - t2 + t28 + t45;
t15 = -t41 * t36 + t68;
t14 = -t41 * t38 - t69;
t5 = -qJD(3) * pkin(3) + t75;
t3 = t46 - t74;
t1 = [qJDD(1) - g(3), -g(3) + (t32 ^ 2 + t34 ^ 2) * qJDD(1), 0, -t53, t54, t53, -t54, t5 * t10 + t3 * t11 - t2 * t12 - t6 * t9 - g(3), 0, 0, 0, 0, 0, t47 * t36 + t86 * t38, -t86 * t36 + t47 * t38; 0, -t51, 0, 0, 0, 0, 0, -t51, 0, 0, 0, 0, 0, t14, -t15; 0, 0, qJDD(3), t48 + t76 - t95, t93, t44 - 0.2e1 * t74 - t76, 0.2e1 * t28 + 0.2e1 * t29 - t93, -t2 * qJ(4) - t3 * pkin(3) - t5 * t8 - g(3) * (t26 * pkin(3) + t25 * qJ(4)) + t75 * t6 + t56 * (pkin(3) * t25 - qJ(4) * t26), t31 * qJDD(3) - 0.2e1 * t64 * t82, -0.2e1 * t36 * t66 + 0.2e1 * t79 * t64, t15, t14, 0, t43 * t36 + t87 * t38, -t87 * t36 + t43 * t38; 0, 0, 0, 0, 0, qJDD(3), -t42, -t74 + t94, 0, 0, 0, 0, 0, t78 * t36 + t68, t78 * t38 - t69; 0, 0, 0, 0, 0, 0, 0, 0, t42 * t82, -t79 * t42, t66, -t36 * qJDD(3), qJDD(5), t51 * t36 - t91 * t38, t91 * t36 + t51 * t38;];
tau_reg = t1;
