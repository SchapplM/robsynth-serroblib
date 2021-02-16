% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tau_reg [4x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:40
% EndTime: 2021-01-15 10:27:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (351->113), mult. (633->137), div. (0->0), fcn. (294->4), ass. (0->77)
t39 = cos(qJ(3));
t41 = -pkin(1) - pkin(5);
t17 = t41 * qJD(1) + qJD(2);
t63 = qJ(4) * qJD(1);
t89 = t17 - t63;
t8 = t89 * t39;
t40 = cos(qJ(1));
t30 = g(2) * t40;
t38 = sin(qJ(1));
t31 = g(1) * t38;
t78 = t31 - t30;
t92 = qJDD(2) - t78;
t37 = sin(qJ(3));
t23 = t37 * pkin(3) + qJ(2);
t70 = qJD(1) * t23;
t13 = qJD(4) + t70;
t91 = -t13 * qJD(1) - t78;
t71 = pkin(1) * qJDD(1);
t90 = t71 - t92;
t33 = qJDD(1) * qJ(2);
t34 = qJD(1) * qJD(2);
t50 = g(1) * t40 + g(2) * t38;
t46 = -t50 + 0.2e1 * t34;
t88 = 0.2e1 * t33 + t46;
t59 = qJD(1) * qJD(3);
t55 = t39 * t59;
t61 = t37 * qJDD(1);
t3 = qJDD(4) + t33 + t34 + (t55 + t61) * pkin(3);
t87 = t3 - t50;
t86 = qJ(4) - t41;
t60 = qJ(4) * qJDD(1);
t85 = (qJD(4) + t13) * qJD(1) + t60;
t74 = qJD(3) * pkin(3);
t4 = t8 + t74;
t82 = t4 - t8;
t28 = g(3) * t37;
t81 = t37 * t17;
t43 = qJD(1) ^ 2;
t80 = t43 * t37;
t79 = g(3) * t39 + t37 * t31;
t35 = t37 ^ 2;
t36 = t39 ^ 2;
t77 = -t35 - t36;
t76 = t35 - t36;
t42 = qJD(3) ^ 2;
t75 = -t42 - t43;
t73 = t43 * qJ(2);
t15 = t86 * t39;
t69 = qJD(3) * t15;
t68 = qJD(3) * t37;
t67 = qJD(3) * t39;
t66 = qJDD(3) * pkin(3);
t62 = qJDD(3) * t37;
t26 = t39 * qJDD(1);
t16 = t41 * qJDD(1) + qJDD(2);
t11 = t39 * t16;
t58 = t39 * t30 + t11 + t28;
t57 = -t16 - t30;
t56 = t37 * t59;
t54 = t13 + t70;
t52 = -0.2e1 * t56;
t48 = -qJD(1) * qJD(4) - t60;
t47 = 0.2e1 * qJ(2) * t59 + qJDD(3) * t41;
t18 = pkin(3) * t67 + qJD(2);
t45 = qJD(1) * t18 + qJDD(1) * t23 + t87;
t44 = -t41 * t42 + t88;
t27 = qJDD(3) * t39;
t19 = qJ(4) * t56;
t14 = t86 * t37;
t10 = t75 * t37 + t27;
t9 = t75 * t39 - t62;
t7 = -t37 * t63 + t81;
t6 = -t37 * qJD(4) - t69;
t5 = -t39 * qJD(4) + t68 * t86;
t2 = t89 * t67 + (t16 + t48) * t37;
t1 = -t17 * t68 + t48 * t39 + t11 + t19 + t66;
t12 = [qJDD(1), t78, t50, -0.2e1 * t71 + t92, t88, t90 * pkin(1) + (t46 + t33) * qJ(2), t36 * qJDD(1) + t39 * t52, -0.2e1 * t37 * t26 + 0.2e1 * t76 * t59, -t42 * t37 + t27, -t42 * t39 - t62, 0, t44 * t37 + t47 * t39, -t47 * t37 + t44 * t39, -t15 * qJDD(3) + (t54 * t39 + t5) * qJD(3) + t45 * t37, t14 * qJDD(3) + (-t54 * t37 - t6) * qJD(3) + t45 * t39, (-t7 * qJD(3) + qJDD(1) * t15 - t1 + (qJD(3) * t14 - t5) * qJD(1)) * t39 + (qJD(3) * t4 + qJDD(1) * t14 - t2 + (-t6 - t69) * qJD(1)) * t37 + t78, -t2 * t14 + t7 * t6 - t1 * t15 + t4 * t5 + t3 * t23 + t13 * t18 - g(1) * (t23 * t40 - t38 * t86) - g(2) * (t23 * t38 + t40 * t86); 0, 0, 0, qJDD(1), -t43, -t73 - t90, 0, 0, 0, 0, 0, t10, t9, t10, t9, t77 * qJDD(1), t1 * t39 + t2 * t37 + (-t37 * t4 + t39 * t7) * qJD(3) + t91; 0, 0, 0, 0, 0, 0, t39 * t80, -t76 * t43, t26, -t61, qJDD(3), (-t73 - t31) * t39 + t58, (t57 + t73) * t37 + t79, 0.2e1 * t66 + t19 + (t7 - t81) * qJD(3) + (-pkin(3) * t80 - t31 - t85) * t39 + t58, -t36 * t43 * pkin(3) + (t57 + t85) * t37 + t79, -pkin(3) * t26 + (t74 - t82) * t37 * qJD(1), t82 * t7 + (t91 * t39 + t1 + t28) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55 + t61, t26 + t52, t77 * t43, (t37 * t7 + t39 * t4) * qJD(1) + t87;];
tau_reg = t12;
