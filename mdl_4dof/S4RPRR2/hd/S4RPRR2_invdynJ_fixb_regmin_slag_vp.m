% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR2
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
% tau_reg [4x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:14
% EndTime: 2019-12-31 16:48:15
% DurationCPUTime: 0.26s
% Computational Cost: add. (319->69), mult. (573->97), div. (0->0), fcn. (328->10), ass. (0->59)
t39 = cos(pkin(7));
t31 = t39 * pkin(1) + pkin(2);
t21 = t31 * qJDD(1);
t32 = qJ(1) + pkin(7) + qJ(3);
t30 = cos(t32);
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t23 = t31 * qJD(1);
t68 = qJD(3) * t23;
t82 = -g(2) * t30 + t44 * t21 - t41 * t68;
t34 = qJDD(1) + qJDD(3);
t77 = t34 * pkin(3);
t62 = qJD(1) * qJD(3) * t44;
t38 = sin(pkin(7));
t79 = pkin(1) * t38;
t80 = (qJDD(1) * t41 + t62) * t79;
t81 = -t77 + t80 - t82;
t70 = t41 * t31 + t44 * t79;
t29 = sin(t32);
t64 = qJD(1) * t79;
t60 = t41 * t64;
t48 = g(1) * t30 + g(2) * t29 - (qJDD(1) * t79 + t68) * t44 + qJD(3) * t60 - t41 * t21;
t26 = g(1) * t29;
t35 = qJD(1) + qJD(3);
t76 = t35 * pkin(3);
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t9 = t44 * t23 - t60;
t7 = -t9 - t76;
t75 = t7 * qJD(4) * t40 + t43 * t26;
t74 = (t41 * t23 + t44 * t64) * t35;
t73 = t70 * qJD(3) * t35;
t72 = t43 * t34;
t71 = t44 * t31;
t36 = t40 ^ 2;
t69 = -t43 ^ 2 + t36;
t67 = t43 * qJD(4);
t66 = qJDD(2) - g(3);
t65 = t81 * t40 + t7 * t67;
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t59 = g(1) * t42 - g(2) * t45;
t57 = -t41 * t79 + t71;
t46 = qJD(4) ^ 2;
t54 = pkin(6) * t46 - t74 - t77;
t14 = -pkin(3) - t57;
t15 = pkin(6) + t70;
t53 = t14 * t34 + t15 * t46 + t73;
t52 = -t34 * pkin(6) - t7 * t35 + t48;
t51 = t26 + t82;
t50 = -pkin(6) * qJDD(4) + (t9 - t76) * qJD(4);
t11 = t57 * qJD(3);
t49 = -qJDD(4) * t15 + (t14 * t35 - t11) * qJD(4);
t33 = t35 ^ 2;
t18 = qJDD(4) * t43 - t46 * t40;
t17 = qJDD(4) * t40 + t46 * t43;
t13 = 0.2e1 * t40 * t35 * t67 + t36 * t34;
t6 = -0.2e1 * t69 * t35 * qJD(4) + 0.2e1 * t40 * t72;
t1 = [qJDD(1), t59, g(1) * t45 + g(2) * t42, (t59 + (t38 ^ 2 + t39 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t34, t34 * t71 - t73 + (-t62 + (-qJDD(1) - t34) * t41) * t79 + t51, -t11 * t35 - t70 * t34 + t48, t13, t6, t17, t18, 0, t49 * t40 + (-t53 - t81) * t43 + t75, t49 * t43 + (t53 - t26) * t40 + t65; 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17; 0, 0, 0, 0, t34, t51 + t74 - t80, t9 * t35 + t48, t13, t6, t17, t18, 0, t50 * t40 + (-t54 - t81) * t43 + t75, t50 * t43 + (t54 - t26) * t40 + t65; 0, 0, 0, 0, 0, 0, 0, -t40 * t33 * t43, t69 * t33, t40 * t34, t72, qJDD(4), t52 * t40 + t66 * t43, -t66 * t40 + t52 * t43;];
tau_reg = t1;
