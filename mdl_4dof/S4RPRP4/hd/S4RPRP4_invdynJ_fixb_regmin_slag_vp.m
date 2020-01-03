% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP4
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:55
% EndTime: 2019-12-31 16:43:56
% DurationCPUTime: 0.33s
% Computational Cost: add. (384->98), mult. (732->133), div. (0->0), fcn. (391->8), ass. (0->68)
t86 = 2 * qJD(3);
t36 = sin(pkin(6));
t24 = t36 * pkin(1) + pkin(5);
t17 = t24 * qJDD(1);
t85 = (qJD(2) * qJD(3)) + t17;
t33 = qJ(1) + pkin(6);
t28 = sin(t33);
t29 = cos(t33);
t57 = g(1) * t29 + g(2) * t28;
t40 = cos(qJ(3));
t19 = t24 * qJD(1);
t38 = sin(qJ(3));
t77 = t38 * t19;
t7 = t40 * qJD(2) - t77;
t84 = qJD(4) - t7;
t4 = -(qJD(3) * pkin(3)) + t84;
t75 = t40 * t19;
t8 = t38 * qJD(2) + t75;
t5 = (qJD(3) * qJ(4)) + t8;
t54 = t40 * pkin(3) + t38 * qJ(4);
t52 = pkin(2) + t54;
t37 = cos(pkin(6));
t78 = t37 * pkin(1);
t10 = -t52 - t78;
t53 = pkin(3) * t38 - qJ(4) * t40;
t9 = t53 * qJD(3) - t38 * qJD(4);
t3 = qJD(1) * t9 + qJDD(1) * t10;
t71 = qJDD(3) * pkin(3);
t83 = qJDD(4) - t71;
t6 = qJD(1) * t10;
t68 = qJDD(3) * t24;
t82 = t6 * t86 - t68;
t81 = g(1) * t28;
t76 = t38 * t40;
t34 = t38 ^ 2;
t35 = t40 ^ 2;
t74 = t34 - t35;
t25 = -pkin(2) - t78;
t20 = qJD(1) * t25;
t72 = qJD(1) * t38;
t67 = t40 * qJDD(1);
t66 = qJD(1) * qJD(3);
t64 = qJDD(3) * qJ(4);
t43 = qJD(1) ^ 2;
t63 = t43 * t76;
t62 = t38 * qJDD(2) + t85 * t40;
t61 = t7 + t77;
t59 = qJD(3) * t75 - t40 * qJDD(2) + t85 * t38;
t39 = sin(qJ(1));
t41 = cos(qJ(1));
t56 = g(1) * t39 - g(2) * t41;
t42 = qJD(3) ^ 2;
t55 = g(2) * t29 + t24 * t42;
t51 = t56 * pkin(1);
t50 = -g(3) * t40 + t57 * t38 - t59;
t49 = -0.2e1 * qJDD(1) * t25 - t55;
t48 = t20 * t86 - t68;
t47 = t8 * qJD(3) + t50;
t46 = -0.2e1 * t3 - t55;
t1 = t64 + (qJD(4) - t77) * qJD(3) + t62;
t2 = t59 + t83;
t45 = t1 * t40 + t2 * t38 + (-t38 * t5 + t4 * t40) * qJD(3);
t31 = t38 * qJDD(1);
t22 = t40 * t81;
t16 = qJDD(3) * t40 - t42 * t38;
t15 = qJDD(3) * t38 + t42 * t40;
t14 = t53 * qJD(1);
t11 = [qJDD(1), t56, g(1) * t41 + g(2) * t39, (t36 ^ 2 + t37 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t51, t34 * qJDD(1) + 0.2e1 * t66 * t76, 0.2e1 * t38 * t67 - 0.2e1 * t74 * t66, t15, t16, 0, t48 * t38 + t49 * t40 + t22, t48 * t40 + (-t49 - t81) * t38, t82 * t38 + t46 * t40 + t22, (t34 + t35) * t17 + t45 - t57, -t82 * t40 + (t46 + t81) * t38, t3 * t10 + t6 * t9 + t51 + (-g(1) * pkin(5) - g(2) * t52) * t29 + (-g(2) * pkin(5) + g(1) * t52) * t28 + t45 * t24; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t16, -t15, t16, 0, t15, t1 * t38 - t2 * t40 - g(3) + (t38 * t4 + t40 * t5) * qJD(3); 0, 0, 0, 0, -t63, t74 * t43, t31, t67, qJDD(3), -t20 * t72 + t47, g(3) * t38 + t61 * qJD(3) + (-qJD(1) * t20 + t57) * t40 - t62, (2 * t71) - qJDD(4) + (t14 * t40 - t38 * t6) * qJD(1) + t47, -t53 * qJDD(1), (2 * t64) + (qJD(1) * t14 - g(3)) * t38 + (qJD(1) * t6 - t57) * t40 + (0.2e1 * qJD(4) - t61) * qJD(3) + t62, -t2 * pkin(3) - g(3) * t54 + t1 * qJ(4) - t6 * t14 - t4 * t8 + t84 * t5 + t57 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t63, t31, -t34 * t43 - t42, -t5 * qJD(3) + t6 * t72 - t50 + t83;];
tau_reg = t11;
