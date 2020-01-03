% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tau_reg [4x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:42
% EndTime: 2019-12-31 16:41:43
% DurationCPUTime: 0.29s
% Computational Cost: add. (299->93), mult. (569->122), div. (0->0), fcn. (372->8), ass. (0->60)
t44 = -pkin(1) - qJ(3);
t42 = sin(pkin(6));
t43 = cos(pkin(6));
t70 = t42 ^ 2 + t43 ^ 2;
t83 = (t44 * qJD(1) + qJD(2)) * t70;
t39 = qJDD(1) * qJ(2);
t40 = qJD(1) * qJD(2);
t80 = t39 + t40;
t21 = qJDD(3) + t80;
t46 = sin(qJ(1));
t48 = cos(qJ(1));
t59 = g(1) * t48 + g(2) * t46;
t82 = t21 - t59;
t81 = g(1) * t46 - g(2) * t48;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t73 = t47 * t43;
t79 = -t45 * t42 + t73;
t15 = t47 * t42 + t45 * t43;
t68 = qJD(1) * t42;
t64 = t45 * t68;
t51 = -qJD(4) * t64 + t15 * qJDD(1);
t78 = -qJD(1) * qJD(3) + qJDD(1) * t44;
t77 = 0.2e1 * t40;
t76 = -pkin(5) + t44;
t11 = t15 * qJD(4);
t75 = -t11 * qJD(4) + qJDD(4) * t79;
t72 = t48 * pkin(1) + t46 * qJ(2);
t69 = pkin(1) * qJDD(1);
t28 = qJD(1) * qJ(2) + qJD(3);
t67 = t42 * qJDD(1);
t66 = t43 * qJDD(1);
t63 = qJD(1) * t73;
t17 = qJDD(2) + t78;
t62 = t70 * t17;
t61 = -pkin(5) * qJDD(1) + t17;
t60 = qJDD(2) - t69;
t57 = -t45 * t67 + t47 * t66;
t18 = t76 * t42;
t19 = t76 * t43;
t56 = t47 * t18 + t45 * t19;
t55 = t45 * t18 - t47 * t19;
t12 = t79 * qJD(4);
t53 = -t12 * qJD(4) - t15 * qJDD(4);
t8 = t15 * qJD(1);
t50 = t80 + t82;
t49 = qJD(1) ^ 2;
t38 = pkin(6) + qJ(4);
t32 = t48 * qJ(2);
t30 = cos(t38);
t29 = sin(t38);
t27 = t42 * pkin(3) + qJ(2);
t20 = pkin(3) * t68 + t28;
t14 = pkin(3) * t67 + t21;
t10 = t63 - t64;
t4 = t61 * t43;
t3 = t61 * t42;
t2 = qJD(4) * t63 + t51;
t1 = -qJD(1) * t11 + t57;
t5 = [qJDD(1), t81, t59, qJDD(2) - 0.2e1 * t69 - t81, 0.2e1 * t39 + t77 - t59, -t60 * pkin(1) - g(1) * (-t46 * pkin(1) + t32) - g(2) * t72 + (t39 + t77) * qJ(2), t50 * t42, t50 * t43, t81 + t70 * (-t17 - t78), t21 * qJ(2) + t28 * qJD(2) - g(1) * (t44 * t46 + t32) - g(2) * (t48 * qJ(3) + t72) + t44 * t62 - qJD(3) * t83, t1 * t79 - t10 * t11, -t1 * t15 - t10 * t12 + t11 * t8 - t2 * t79, t75, t53, 0, qJD(2) * t8 + t27 * t2 + t14 * t15 + t20 * t12 - t55 * qJDD(4) - t59 * t29 + (-qJD(3) * t79 - t56 * qJD(4)) * qJD(4), qJD(2) * t10 + t27 * t1 + t14 * t79 - t20 * t11 - t56 * qJDD(4) - t59 * t30 + (t15 * qJD(3) + t55 * qJD(4)) * qJD(4); 0, 0, 0, qJDD(1), -t49, -t49 * qJ(2) + t60 - t81, -t49 * t42, -t49 * t43, -t70 * qJDD(1), -t28 * qJD(1) + t62 - t81, 0, 0, 0, 0, 0, -qJD(1) * t8 + t75, -qJD(1) * t10 + t53; 0, 0, 0, 0, 0, 0, t67, t66, -t70 * t49, qJD(1) * t83 + t82, 0, 0, 0, 0, 0, (t10 + t63) * qJD(4) + t51, -0.2e1 * t8 * qJD(4) + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t8, t10 ^ 2 - t8 ^ 2, t57, (t10 - t63) * qJD(4) - t51, qJDD(4), g(3) * t29 - t20 * t10 - t45 * t3 - t30 * t81 + t47 * t4, g(3) * t30 + t20 * t8 + t29 * t81 - t47 * t3 - t45 * t4;];
tau_reg = t5;
