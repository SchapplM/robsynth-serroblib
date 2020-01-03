% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR5
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tau_reg [4x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:34
% EndTime: 2019-12-31 17:03:35
% DurationCPUTime: 0.30s
% Computational Cost: add. (327->99), mult. (443->118), div. (0->0), fcn. (214->8), ass. (0->68)
t31 = qJD(1) + qJD(2);
t36 = sin(qJ(2));
t61 = qJDD(1) * t36;
t39 = cos(qJ(2));
t67 = qJD(1) * t39;
t34 = qJ(1) + qJ(2);
t27 = sin(t34);
t28 = cos(t34);
t74 = g(1) * t28 + g(2) * t27;
t84 = -((qJD(2) - t31) * t67 + t61) * pkin(1) + t74;
t83 = -g(1) * t27 + g(2) * t28;
t70 = pkin(1) * qJD(1);
t56 = t36 * t70;
t68 = t31 * qJ(3);
t10 = t56 + t68;
t82 = -t10 * t31 + t83;
t29 = t31 ^ 2;
t41 = -pkin(2) - pkin(6);
t35 = sin(qJ(4));
t65 = qJD(2) * t39;
t45 = (qJD(1) * t65 + t61) * pkin(1);
t63 = t31 * qJD(3);
t30 = qJDD(1) + qJDD(2);
t69 = t30 * qJ(3);
t5 = t45 + t63 + t69;
t38 = cos(qJ(4));
t64 = qJD(4) * t38;
t80 = t10 * t64 + t5 * t35;
t79 = t30 * pkin(2);
t78 = t39 * pkin(1);
t76 = t38 * t30;
t75 = t28 * pkin(2) + t27 * qJ(3);
t73 = -qJD(2) * t56 + qJDD(1) * t78;
t42 = qJD(4) ^ 2;
t72 = -t29 - t42;
t33 = t38 ^ 2;
t71 = t35 ^ 2 - t33;
t66 = qJD(2) * t36;
t24 = -pkin(2) - t78;
t15 = -pkin(6) + t24;
t60 = qJDD(4) * t15;
t59 = qJDD(4) * t35;
t58 = qJDD(4) * t41;
t57 = pkin(1) * t66;
t55 = t31 * t66;
t54 = qJDD(3) - t73;
t53 = -t27 * pkin(2) + t28 * qJ(3);
t52 = -t73 + t83;
t50 = -pkin(1) * t67 + qJD(3);
t6 = t54 - t79;
t18 = t36 * pkin(1) + qJ(3);
t49 = t18 * t31 + t57;
t48 = -t56 + t68;
t47 = -t41 * t30 - t54 - t82;
t46 = t31 * t56 - t52;
t13 = pkin(1) * t65 + qJD(3);
t44 = t13 * t31 - t15 * t42 + t18 * t30 - t74;
t43 = t50 * t31 - t41 * t42 + t69 - t74;
t40 = cos(qJ(1));
t37 = sin(qJ(1));
t26 = qJDD(4) * t38;
t12 = -t42 * t35 + t26;
t11 = -t42 * t38 - t59;
t9 = -t31 * pkin(2) + t50;
t7 = -0.2e1 * t31 * t35 * t64 + t33 * t30;
t3 = t5 * t38;
t1 = 0.2e1 * t71 * t31 * qJD(4) - 0.2e1 * t35 * t76;
t2 = [qJDD(1), g(1) * t37 - g(2) * t40, g(1) * t40 + g(2) * t37, t30, (t30 * t39 - t55) * pkin(1) - t52, ((-qJDD(1) - t30) * t36 + (-qJD(1) - t31) * t65) * pkin(1) + t74, pkin(1) * t55 + qJDD(3) + (-pkin(2) + t24) * t30 + t52, (qJD(3) + t13) * t31 + (qJ(3) + t18) * t30 + t45 - t74, t5 * t18 + t10 * t13 + t6 * t24 + t9 * t57 - g(1) * (-t37 * pkin(1) + t53) - g(2) * (t40 * pkin(1) + t75), t7, t1, t12, t11, 0, (t49 * qJD(4) + t60) * t38 + t44 * t35 + t80, t3 + (-t60 + (-t10 - t49) * qJD(4)) * t35 + t44 * t38; 0, 0, 0, t30, t46, t84, qJDD(3) - t46 - 0.2e1 * t79, 0.2e1 * t63 + 0.2e1 * t69 - t84, t5 * qJ(3) + t10 * qJD(3) - t6 * pkin(2) - g(1) * t53 - g(2) * t75 + (-t10 * t39 - t36 * t9) * t70, t7, t1, t12, t11, 0, (t48 * qJD(4) + t58) * t38 + t43 * t35 + t80, t3 + (-t58 + (-t10 - t48) * qJD(4)) * t35 + t43 * t38; 0, 0, 0, 0, 0, 0, t30, -t29, t6 + t82, 0, 0, 0, 0, 0, t72 * t35 + t26, t72 * t38 - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t29 * t35, -t71 * t29, t76, -t35 * t30, qJDD(4), g(3) * t35 - t47 * t38, g(3) * t38 + t47 * t35;];
tau_reg = t2;
