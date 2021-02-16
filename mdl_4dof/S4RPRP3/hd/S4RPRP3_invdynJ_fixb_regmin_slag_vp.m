% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP3
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
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:36
% EndTime: 2021-01-15 10:20:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (386->123), mult. (766->155), div. (0->0), fcn. (413->8), ass. (0->81)
t88 = qJDD(2) - g(3);
t44 = cos(qJ(3));
t39 = sin(pkin(6));
t27 = t39 * pkin(1) + pkin(5);
t18 = t27 * qJD(1);
t69 = qJ(4) * qJD(1);
t62 = t18 + t69;
t55 = t62 * t44;
t34 = t44 * qJD(2);
t42 = sin(qJ(3));
t5 = -t62 * t42 + t34;
t74 = qJD(3) * pkin(3);
t3 = t5 + t74;
t87 = t3 - t5;
t37 = t42 ^ 2;
t86 = pkin(3) * t37;
t36 = qJ(1) + pkin(6);
t30 = sin(t36);
t85 = g(1) * t30;
t31 = cos(t36);
t84 = g(1) * t31;
t83 = g(2) * t30;
t82 = g(2) * t31;
t81 = g(3) * t44;
t40 = cos(pkin(6));
t80 = t40 * pkin(1);
t79 = t44 * pkin(3);
t78 = t30 * t44;
t77 = t31 * t42;
t47 = qJD(1) ^ 2;
t76 = t44 * t47;
t38 = t44 ^ 2;
t75 = t37 - t38;
t73 = qJ(4) + t27;
t28 = -pkin(2) - t80;
t19 = qJD(1) * t28;
t72 = qJDD(3) * pkin(3);
t71 = t19 * qJD(1);
t70 = t42 * qJD(3);
t29 = pkin(2) + t79;
t13 = -t29 - t80;
t68 = qJDD(1) * t13;
t32 = t42 * qJDD(1);
t67 = t44 * qJDD(1);
t66 = qJD(1) * qJD(3);
t65 = t42 * t66;
t4 = pkin(3) * t65 + qJDD(4) + t68;
t64 = t4 + t68;
t63 = qJD(3) * t73;
t61 = 0.2e1 * t44 * t66;
t16 = t27 * qJDD(1);
t60 = -qJD(3) * qJD(2) - t16;
t59 = -t83 - t84;
t43 = sin(qJ(1));
t45 = cos(qJ(1));
t58 = g(1) * t43 - g(2) * t45;
t6 = t42 * qJD(2) + t55;
t57 = t3 * t42 - t44 * t6;
t33 = t44 * qJDD(2);
t56 = g(1) * t77 + t42 * t83 + t33 - t81;
t46 = qJD(3) ^ 2;
t54 = 0.2e1 * qJDD(1) * t28 + t27 * t46;
t12 = t18 * t70;
t53 = g(2) * t78 - t88 * t42 + t44 * t84 + t12;
t52 = -qJ(4) * qJDD(1) + t60;
t51 = 0.2e1 * qJD(3) * t19 - qJDD(3) * t27;
t50 = qJD(1) * qJD(4) - t52;
t9 = qJD(1) * t13 + qJD(4);
t49 = (-qJD(4) - t9) * qJD(1) + t52;
t41 = -qJ(4) - pkin(5);
t23 = g(1) * t78;
t22 = g(2) * t77;
t15 = qJDD(3) * t44 - t46 * t42;
t14 = qJDD(3) * t42 + t46 * t44;
t11 = t73 * t44;
t10 = t73 * t42;
t8 = -t42 * qJD(4) - t44 * t63;
t7 = t44 * qJD(4) - t42 * t63;
t2 = -t12 + (-qJ(4) * t66 + qJDD(2)) * t42 + t50 * t44;
t1 = -qJD(3) * t55 - t50 * t42 + t33 + t72;
t17 = [qJDD(1), t58, g(1) * t45 + g(2) * t43, (t58 + (t39 ^ 2 + t40 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t37 * qJDD(1) + t42 * t61, 0.2e1 * t42 * t67 - 0.2e1 * t75 * t66, t14, t15, 0, t23 + t51 * t42 + (-t54 - t82) * t44, t22 + t51 * t44 + (t54 - t85) * t42, -t10 * qJDD(3) + t23 + (-t64 - t82) * t44 + (t8 + (t9 + (t13 - t79) * qJD(1)) * t42) * qJD(3), -t11 * qJDD(3) + t22 + (t64 - t85) * t42 + (t44 * t9 - t7 + (t13 * t44 + t86) * qJD(1)) * qJD(3), (-qJD(3) * t3 + qJDD(1) * t11 + t2 + (qJD(3) * t10 + t7) * qJD(1)) * t44 + (-t6 * qJD(3) + qJDD(1) * t10 - t1 + (-qJD(3) * t11 - t8) * qJD(1)) * t42 + t59, t2 * t11 + t6 * t7 - t1 * t10 + t3 * t8 + t4 * t13 + t9 * pkin(3) * t70 - g(1) * (-t43 * pkin(1) - t30 * t29 - t31 * t41) - g(2) * (t45 * pkin(1) + t31 * t29 - t30 * t41); 0, 0, 0, t88, 0, 0, 0, 0, 0, t15, -t14, t15, -t14, 0, -t57 * qJD(3) + t1 * t44 + t2 * t42 - g(3); 0, 0, 0, 0, -t42 * t76, t75 * t47, t32, t67, qJDD(3), (-t16 - t71) * t42 + t56, (-t42 * t18 + t34) * qJD(3) + (t60 - t71) * t44 + t53, 0.2e1 * t72 + (t6 - t55) * qJD(3) + (pkin(3) * t76 + t49) * t42 + t56, -t47 * t86 + (t42 * t69 + t5) * qJD(3) + t49 * t44 + t53, -pkin(3) * t32 + (-t74 + t87) * t44 * qJD(1), t87 * t6 + (-t81 + t1 + (-qJD(1) * t9 - t59) * t42) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65 - t67, t32 + t61, (-t37 - t38) * t47, t57 * qJD(1) + t4 + t82 - t85;];
tau_reg = t17;
