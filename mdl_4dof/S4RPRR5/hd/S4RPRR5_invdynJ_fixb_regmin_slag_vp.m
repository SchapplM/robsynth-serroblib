% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tau_reg [4x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:41
% EndTime: 2019-12-31 16:51:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (366->83), mult. (560->106), div. (0->0), fcn. (322->6), ass. (0->58)
t41 = cos(qJ(1));
t77 = sin(qJ(1));
t90 = g(1) * t77 - g(2) * t41;
t38 = sin(qJ(3));
t43 = qJD(4) ^ 2;
t33 = qJDD(1) - qJDD(3);
t40 = cos(qJ(3));
t74 = t40 * t33;
t65 = qJD(1) - qJD(3);
t88 = t65 ^ 2;
t89 = (t43 + t88) * t38 + t74;
t42 = -pkin(1) - pkin(2);
t68 = qJ(2) * qJD(1);
t87 = -qJD(3) * t68 + t42 * qJDD(1) + qJDD(2);
t86 = -qJDD(2) + t90;
t69 = qJD(4) * t65;
t20 = t42 * qJD(1) + qJD(2);
t66 = (qJD(1) * qJD(2));
t67 = (qJ(2) * qJDD(1));
t84 = qJD(3) * t20 + t66 + t67;
t73 = t40 * qJ(2) + t38 * t42;
t71 = pkin(1) * qJDD(1);
t83 = t71 + t86;
t12 = -t77 * t38 - t41 * t40;
t13 = t41 * t38 - t77 * t40;
t45 = g(1) * t12 + g(2) * t13 + t87 * t38 + t84 * t40;
t48 = -g(1) * t13 + g(2) * t12 + t84 * t38 - t87 * t40;
t81 = t33 * pkin(3);
t80 = t65 * pkin(3);
t79 = (t38 * qJD(2) + t73 * qJD(3)) * t65;
t78 = (t38 * t20 + t40 * t68) * t65;
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t76 = t37 * t39;
t75 = t39 * t33;
t35 = t37 ^ 2;
t72 = -t39 ^ 2 + t35;
t59 = -t38 * qJ(2) + t40 * t42;
t58 = -t81 - t48;
t57 = g(1) * t41 + g(2) * t77;
t8 = t40 * t20 - t38 * t68;
t4 = -t8 + t80;
t54 = t33 * pkin(6) + t4 * t65 - t45;
t52 = -pkin(6) * qJDD(4) + (t4 + t8 + t80) * qJD(4);
t51 = -t57 + (2 * t66);
t14 = pkin(3) - t59;
t15 = -pkin(6) + t73;
t6 = t40 * qJD(2) + t59 * qJD(3);
t50 = -qJDD(4) * t15 + (-t14 * t65 - t4 - t6) * qJD(4);
t49 = -qJDD(4) * t38 + 0.2e1 * t40 * t69;
t47 = pkin(6) * t43 - t58 + t78 + t81;
t46 = -t14 * t33 + t15 * t43 + t58 - t79;
t44 = qJD(1) ^ 2;
t18 = qJDD(4) * t39 - t43 * t37;
t17 = qJDD(4) * t37 + t43 * t39;
t5 = -t35 * t33 - 0.2e1 * t69 * t76;
t3 = -t37 * t75 + t72 * t69;
t1 = [qJDD(1), t90, t57, 0.2e1 * t71 + t86, t51 + (2 * t67), t83 * pkin(1) + (t51 + t67) * qJ(2), t33, -t59 * t33 + t48 + t79, t73 * t33 + t6 * t65 + t45, -t5, -0.2e1 * t3, -t17, -t18, 0, t50 * t37 - t46 * t39, t46 * t37 + t50 * t39; 0, 0, 0, -qJDD(1), -t44, -t44 * qJ(2) - t83, 0, -t38 * t88 - t74, t38 * t33 - t40 * t88, 0, 0, 0, 0, 0, t49 * t37 - t89 * t39, t89 * t37 + t49 * t39; 0, 0, 0, 0, 0, 0, -t33, -t48 - t78, -t65 * t8 - t45, t5, 0.2e1 * t3, t17, t18, 0, t52 * t37 - t47 * t39, t47 * t37 + t52 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t76, t72 * t88, -t37 * t33, -t75, qJDD(4), g(3) * t39 + t54 * t37, -g(3) * t37 + t54 * t39;];
tau_reg = t1;
