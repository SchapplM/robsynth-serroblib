% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:40
% EndTime: 2019-12-31 17:56:41
% DurationCPUTime: 0.37s
% Computational Cost: add. (573->101), mult. (804->130), div. (0->0), fcn. (455->10), ass. (0->64)
t45 = sin(pkin(8));
t33 = pkin(1) * t45 + qJ(3);
t81 = t33 * qJDD(1);
t48 = sin(qJ(4));
t53 = qJD(5) ^ 2;
t39 = qJDD(1) - qJDD(4);
t51 = cos(qJ(4));
t83 = t51 * t39;
t74 = qJD(1) - qJD(4);
t92 = t74 ^ 2;
t93 = (t53 + t92) * t48 + t83;
t77 = qJD(5) * t74;
t46 = cos(pkin(8));
t34 = -pkin(1) * t46 - pkin(2);
t32 = -pkin(3) + t34;
t82 = t48 * t32 + t51 * t33;
t75 = qJ(1) + pkin(8);
t36 = cos(t75);
t73 = sin(t75);
t16 = -t36 * t51 - t73 * t48;
t17 = t36 * t48 - t73 * t51;
t20 = t32 * qJDD(1) + qJDD(3);
t22 = t32 * qJD(1) + qJD(3);
t76 = qJD(3) * qJD(1);
t23 = t76 + t81;
t27 = t33 * qJD(1);
t78 = qJD(4) * t51;
t79 = qJD(4) * t48;
t57 = g(1) * t16 + g(2) * t17 + t48 * t20 + t22 * t78 + t51 * t23 - t27 * t79;
t58 = -g(1) * t17 + g(2) * t16 - t51 * t20 + t22 * t79 + t48 * t23 + t27 * t78;
t67 = t34 * qJDD(1);
t89 = pkin(4) * t39;
t88 = pkin(4) * t74;
t87 = t74 * (t22 * t48 + t27 * t51);
t86 = (t48 * qJD(3) + t82 * qJD(4)) * t74;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t85 = t47 * t50;
t84 = t50 * t39;
t42 = t47 ^ 2;
t80 = -t50 ^ 2 + t42;
t44 = qJDD(2) - g(3);
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t69 = g(1) * t49 - g(2) * t52;
t5 = t22 * t51 - t27 * t48;
t68 = t32 * t51 - t33 * t48;
t65 = -t89 - t58;
t3 = -t5 + t88;
t63 = pkin(7) * t39 + t3 * t74 - t57;
t62 = -pkin(7) * qJDD(5) + (t3 + t5 + t88) * qJD(5);
t10 = pkin(4) - t68;
t11 = -pkin(7) + t82;
t7 = t51 * qJD(3) + t68 * qJD(4);
t61 = -qJDD(5) * t11 + (-t10 * t74 - t3 - t7) * qJD(5);
t60 = g(1) * t73 - g(2) * t36 - qJDD(3);
t59 = -qJDD(5) * t48 + 0.2e1 * t51 * t77;
t56 = pkin(7) * t53 - t65 + t87 + t89;
t55 = -t10 * t39 + t11 * t53 + t65 - t86;
t25 = qJDD(5) * t50 - t53 * t47;
t24 = qJDD(5) * t47 + t50 * t53;
t12 = -t39 * t42 - 0.2e1 * t77 * t85;
t9 = -t47 * t84 + t80 * t77;
t1 = [qJDD(1), t69, g(1) * t52 + g(2) * t49, (t69 + (t45 ^ 2 + t46 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t60 - 0.2e1 * t67, -g(1) * t36 - g(2) * t73 + 0.2e1 * t76 + 0.2e1 * t81, t23 * t33 + t27 * qJD(3) + (qJDD(3) + t67) * t34 - g(1) * (-t49 * pkin(1) - t73 * pkin(2) + t36 * qJ(3)) - g(2) * (t52 * pkin(1) + t36 * pkin(2) + t73 * qJ(3)), t39, -t68 * t39 + t58 + t86, t82 * t39 + t7 * t74 + t57, -t12, -0.2e1 * t9, -t24, -t25, 0, t61 * t47 - t55 * t50, t55 * t47 + t61 * t50; 0, 0, 0, t44, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24; 0, 0, 0, 0, -qJDD(1), -qJD(1) ^ 2, -t27 * qJD(1) - t60 + t67, 0, -t48 * t92 - t83, t48 * t39 - t51 * t92, 0, 0, 0, 0, 0, t59 * t47 - t93 * t50, t93 * t47 + t59 * t50; 0, 0, 0, 0, 0, 0, 0, -t39, -t58 - t87, -t5 * t74 - t57, t12, 0.2e1 * t9, t24, t25, 0, t62 * t47 - t56 * t50, t56 * t47 + t62 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 * t85, t80 * t92, -t47 * t39, -t84, qJDD(5), -t44 * t50 + t63 * t47, t44 * t47 + t63 * t50;];
tau_reg = t1;
