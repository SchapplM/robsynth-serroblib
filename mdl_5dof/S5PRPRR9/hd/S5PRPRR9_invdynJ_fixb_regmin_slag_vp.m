% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:47
% EndTime: 2019-12-31 17:39:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (454->86), mult. (606->106), div. (0->0), fcn. (344->6), ass. (0->60)
t69 = pkin(8) + qJ(2);
t32 = cos(t69);
t66 = sin(t69);
t94 = g(1) * t66 - g(2) * t32;
t42 = sin(qJ(4));
t46 = qJD(5) ^ 2;
t36 = qJDD(2) - qJDD(4);
t44 = cos(qJ(4));
t79 = t44 * t36;
t70 = qJD(2) - qJD(4);
t92 = t70 ^ 2;
t93 = (t46 + t92) * t42 + t79;
t45 = -pkin(2) - pkin(3);
t73 = qJ(3) * qJD(2);
t91 = -qJD(4) * t73 + t45 * qJDD(2) + qJDD(3);
t90 = -qJDD(3) + t94;
t74 = qJD(5) * t70;
t24 = t45 * qJD(2) + qJD(3);
t71 = (qJD(2) * qJD(3));
t72 = (qJ(3) * qJDD(2));
t88 = qJD(4) * t24 + t71 + t72;
t76 = pkin(2) * qJDD(2);
t87 = t76 + t90;
t78 = t44 * qJ(3) + t42 * t45;
t11 = -t32 * t44 - t66 * t42;
t12 = t32 * t42 - t66 * t44;
t48 = g(1) * t11 + g(2) * t12 + t91 * t42 + t88 * t44;
t52 = -g(1) * t12 + g(2) * t11 + t88 * t42 - t91 * t44;
t85 = t36 * pkin(4);
t84 = t70 * pkin(4);
t83 = (t42 * qJD(3) + t78 * qJD(4)) * t70;
t82 = (t42 * t24 + t44 * t73) * t70;
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t81 = t41 * t43;
t80 = t43 * t36;
t38 = t41 ^ 2;
t77 = -t43 ^ 2 + t38;
t40 = qJDD(1) - g(3);
t62 = -t42 * qJ(3) + t44 * t45;
t61 = -t85 - t52;
t9 = t44 * t24 - t42 * t73;
t4 = -t9 + t84;
t58 = t36 * pkin(7) + t4 * t70 - t48;
t57 = g(1) * t32 + g(2) * t66;
t55 = -pkin(7) * qJDD(5) + (t4 + t9 + t84) * qJD(5);
t15 = pkin(4) - t62;
t16 = -pkin(7) + t78;
t7 = t44 * qJD(3) + t62 * qJD(4);
t54 = -qJDD(5) * t16 + (-t15 * t70 - t4 - t7) * qJD(5);
t53 = -qJDD(5) * t42 + 0.2e1 * t44 * t74;
t51 = -t57 + (2 * t71);
t50 = pkin(7) * t46 - t61 + t82 + t85;
t49 = -t15 * t36 + t16 * t46 + t61 - t83;
t47 = qJD(2) ^ 2;
t19 = qJDD(5) * t43 - t46 * t41;
t18 = qJDD(5) * t41 + t46 * t43;
t6 = -t38 * t36 - 0.2e1 * t74 * t81;
t3 = -t41 * t80 + t77 * t74;
t1 = [t40, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18; 0, qJDD(2), t94, t57, 0.2e1 * t76 + t90, t51 + (2 * t72), t87 * pkin(2) + (t51 + t72) * qJ(3), t36, -t62 * t36 + t52 + t83, t78 * t36 + t7 * t70 + t48, -t6, -0.2e1 * t3, -t18, -t19, 0, t54 * t41 - t49 * t43, t49 * t41 + t54 * t43; 0, 0, 0, 0, -qJDD(2), -t47, -t47 * qJ(3) - t87, 0, -t42 * t92 - t79, t42 * t36 - t44 * t92, 0, 0, 0, 0, 0, t53 * t41 - t93 * t43, t93 * t41 + t53 * t43; 0, 0, 0, 0, 0, 0, 0, -t36, -t52 - t82, -t70 * t9 - t48, t6, 0.2e1 * t3, t18, t19, 0, t55 * t41 - t50 * t43, t50 * t41 + t55 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 * t81, t77 * t92, -t41 * t36, -t80, qJDD(5), -t40 * t43 + t58 * t41, t40 * t41 + t58 * t43;];
tau_reg = t1;
