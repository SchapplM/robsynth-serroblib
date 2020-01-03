% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR6
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:42
% EndTime: 2019-12-31 16:24:44
% DurationCPUTime: 0.43s
% Computational Cost: add. (324->107), mult. (746->152), div. (0->0), fcn. (553->10), ass. (0->73)
t42 = cos(qJ(2));
t82 = g(3) * t42;
t40 = sin(qJ(2));
t36 = sin(pkin(6));
t38 = cos(pkin(6));
t57 = g(1) * t38 + g(2) * t36;
t89 = t57 * t40;
t46 = t89 - t82;
t72 = qJDD(1) - g(3);
t91 = t72 * t42 + t89;
t35 = sin(pkin(7));
t37 = cos(pkin(7));
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t17 = t41 * t35 + t39 * t37;
t86 = t17 * qJD(4);
t90 = qJD(2) * t86;
t88 = t57 * t42;
t12 = t17 * qJD(2);
t73 = t42 * qJD(1);
t58 = qJD(3) - t73;
t66 = qJD(1) * qJD(2);
t75 = qJDD(2) * pkin(2);
t67 = t40 * t66 + qJDD(3);
t53 = -t42 * qJDD(1) + t67;
t9 = t53 - t75;
t85 = (t57 + t66) * t40 + t75 - t9 - t82;
t84 = qJD(4) ^ 2;
t83 = g(3) * t40;
t81 = t36 * t42;
t80 = t38 * t42;
t79 = t39 * t35;
t78 = t41 * t37;
t77 = pkin(5) + qJ(3);
t76 = t35 ^ 2 + t37 ^ 2;
t74 = t40 * qJD(1);
t16 = -t78 + t79;
t71 = t16 * qJDD(4);
t70 = t17 * qJDD(4);
t69 = t35 * qJDD(2);
t68 = t37 * qJDD(2);
t63 = qJD(2) * t78;
t65 = qJD(4) * t63 + t39 * t68 + t41 * t69;
t64 = qJD(2) * t79;
t28 = -t37 * pkin(3) - pkin(2);
t8 = qJDD(2) * qJ(3) + t40 * qJDD(1) + (qJD(3) + t73) * qJD(2);
t62 = t76 * t8;
t61 = t76 * t42;
t60 = pkin(5) * qJDD(2) + t8;
t59 = t76 * qJDD(2);
t56 = t39 * t69 - t41 * t68;
t18 = t77 * t35;
t19 = t77 * t37;
t55 = -t41 * t18 - t39 * t19;
t54 = -t39 * t18 + t41 * t19;
t43 = qJD(2) ^ 2;
t52 = t42 * qJDD(2) - t43 * t40;
t13 = t16 * qJD(4);
t45 = t58 * t76;
t44 = t62 - t83 - t88;
t34 = pkin(7) + qJ(4);
t31 = cos(t34);
t30 = sin(t34);
t22 = qJD(2) * qJ(3) + t74;
t20 = -qJD(2) * pkin(2) + t58;
t15 = t28 * qJD(2) + t58;
t10 = -t63 + t64;
t5 = t28 * qJDD(2) + t53;
t4 = t60 * t37;
t3 = t60 * t35;
t2 = t56 + t90;
t1 = -qJD(4) * t64 + t65;
t6 = [t72, 0, t52, -qJDD(2) * t40 - t43 * t42, t52 * t37, -t52 * t35, t40 * t59 + t43 * t61, -t9 * t42 - g(3) + t40 * t62 + (t20 * t40 + t22 * t61) * qJD(2), 0, 0, 0, 0, 0, (-t2 - t90) * t42 + (qJD(2) * t10 + t16 * t84 - t70) * t40, (qJD(2) * t13 - t1) * t42 + (qJD(2) * t12 + t17 * t84 + t71) * t40; 0, qJDD(2), t91, -t72 * t40 + t88, t85 * t37, -t85 * t35, qJ(3) * t59 + t45 * qJD(2) + t44, -t20 * t74 + (-t9 + t46) * pkin(2) + t44 * qJ(3) + t45 * t22, t1 * t17 - t12 * t13, -t1 * t16 + t13 * t10 - t12 * t86 - t17 * t2, -t13 * qJD(4) + t70, -qJD(4) * t86 - t71, 0, t55 * qJDD(4) + t28 * t2 + t5 * t16 + t15 * t86 + (-t17 * qJD(3) - t54 * qJD(4)) * qJD(4) + t46 * t31 + (-t40 * t10 + t42 * t86) * qJD(1), -t54 * qJDD(4) + t28 * t1 + t5 * t17 - t15 * t13 + (t16 * qJD(3) - t55 * qJD(4)) * qJD(4) - t46 * t30 + (-t40 * t12 - t42 * t13) * qJD(1); 0, 0, 0, 0, -t68, t69, -t76 * t43, -t76 * t22 * qJD(2) + t67 - t75 - t91, 0, 0, 0, 0, 0, 0.2e1 * t12 * qJD(4) + t56, (-t10 - t64) * qJD(4) + t65; 0, 0, 0, 0, 0, 0, 0, 0, t12 * t10, -t10 ^ 2 + t12 ^ 2, (t10 - t64) * qJD(4) + t65, -t56, qJDD(4), -t39 * t4 - t41 * t3 - t15 * t12 - g(1) * (-t30 * t80 + t36 * t31) - g(2) * (-t30 * t81 - t38 * t31) + t30 * t83, -t41 * t4 + t39 * t3 + t15 * t10 - g(1) * (-t36 * t30 - t31 * t80) - g(2) * (t38 * t30 - t31 * t81) + t31 * t83;];
tau_reg = t6;
