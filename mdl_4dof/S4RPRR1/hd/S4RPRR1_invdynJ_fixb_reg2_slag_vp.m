% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR1
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:52
% EndTime: 2019-03-08 18:31:52
% DurationCPUTime: 0.27s
% Computational Cost: add. (629->97), mult. (1256->123), div. (0->0), fcn. (726->14), ass. (0->69)
t53 = cos(pkin(7));
t40 = t53 * pkin(1) + pkin(2);
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t52 = sin(pkin(7));
t86 = pkin(1) * t52;
t21 = t58 * t40 - t55 * t86;
t50 = qJ(1) + pkin(7);
t46 = qJ(3) + t50;
t41 = qJ(4) + t46;
t35 = sin(t41);
t36 = cos(t41);
t89 = g(1) * t36 + g(2) * t35;
t88 = g(1) * t35 - g(2) * t36;
t38 = sin(t46);
t39 = cos(t46);
t87 = g(1) * t38 - g(2) * t39;
t79 = pkin(1) * qJDD(1);
t25 = t40 * qJDD(1);
t76 = qJD(1) * t86;
t69 = t55 * t76;
t75 = t52 * t79;
t27 = t40 * qJD(1);
t78 = qJD(3) * t27;
t11 = (t75 + t78) * t58 - qJD(3) * t69 + t55 * t25;
t16 = t58 * t27 - t69;
t49 = qJD(1) + qJD(3);
t14 = t49 * pkin(3) + t16;
t17 = t55 * t27 + t58 * t76;
t54 = sin(qJ(4));
t82 = t54 * t17;
t15 = qJD(4) * t82;
t57 = cos(qJ(4));
t48 = qJDD(1) + qJDD(3);
t62 = -t55 * t78 + t58 * t25 + (-qJD(1) * qJD(3) * t58 - qJDD(1) * t55) * t86;
t8 = t48 * pkin(3) + t62;
t1 = t54 * t8 - t15 + (qJD(4) * t14 + t11) * t57;
t42 = qJDD(4) + t48;
t85 = pkin(3) * t42;
t81 = t57 * t17;
t44 = cos(t50);
t59 = cos(qJ(1));
t80 = t59 * pkin(1) + pkin(2) * t44;
t45 = qJD(4) + t49;
t74 = -pkin(3) * t45 - t14;
t73 = -t54 * t11 + t57 * t8;
t43 = sin(t50);
t56 = sin(qJ(1));
t68 = -t56 * pkin(1) - pkin(2) * t43;
t67 = g(1) * t56 - g(2) * t59;
t6 = t54 * t14 + t81;
t20 = pkin(3) + t21;
t22 = t55 * t40 + t58 * t86;
t12 = t57 * t20 - t54 * t22;
t13 = t54 * t20 + t57 * t22;
t2 = -qJD(4) * t6 + t73;
t65 = -t1 + t89;
t64 = g(1) * t39 + g(2) * t38 - t11;
t63 = t2 + t88;
t61 = t62 + t87;
t51 = qJDD(2) - g(3);
t19 = t22 * qJD(3);
t18 = t21 * qJD(3);
t10 = t57 * t16 - t82;
t9 = -t54 * t16 - t81;
t5 = t57 * t14 - t82;
t4 = -t13 * qJD(4) - t54 * t18 - t57 * t19;
t3 = t12 * qJD(4) + t57 * t18 - t54 * t19;
t7 = [0, 0, 0, 0, 0, qJDD(1), t67, g(1) * t59 + g(2) * t56, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t43 - g(2) * t44 + 0.2e1 * t53 * t79, g(1) * t44 + g(2) * t43 - 0.2e1 * t75, 0 (t67 + (t52 ^ 2 + t53 ^ 2) * t79) * pkin(1), 0, 0, 0, 0, 0, t48, -t19 * t49 + t21 * t48 + t61, -t18 * t49 - t22 * t48 + t64, 0, -g(1) * t68 - g(2) * t80 + t11 * t22 - t16 * t19 + t17 * t18 + t62 * t21, 0, 0, 0, 0, 0, t42, t12 * t42 + t4 * t45 + t63, -t13 * t42 - t3 * t45 + t65, 0, t1 * t13 + t6 * t3 + t2 * t12 + t5 * t4 - g(1) * (-pkin(3) * t38 + t68) - g(2) * (pkin(3) * t39 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t17 * t49 + t61, t16 * t49 + t64, 0, 0, 0, 0, 0, 0, 0, t42, t57 * t85 - t9 * t45 + (t74 * t54 - t81) * qJD(4) + t73 + t88, t10 * t45 + t15 + (-t8 - t85) * t54 + (t74 * qJD(4) - t11) * t57 + t89, 0, -t6 * t10 - t5 * t9 + (t1 * t54 + t2 * t57 + (-t5 * t54 + t57 * t6) * qJD(4) + t87) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t6 * t45 + t63, t5 * t45 + t65, 0, 0;];
tau_reg  = t7;
