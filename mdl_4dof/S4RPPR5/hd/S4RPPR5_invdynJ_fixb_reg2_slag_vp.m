% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR5
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:51
% EndTime: 2019-12-31 16:39:52
% DurationCPUTime: 0.38s
% Computational Cost: add. (507->125), mult. (872->157), div. (0->0), fcn. (495->6), ass. (0->77)
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t44 = cos(pkin(6));
t74 = qJD(1) * qJD(2);
t27 = t44 * t74;
t49 = -pkin(1) - pkin(2);
t23 = t49 * qJDD(1) + qJDD(2);
t43 = sin(pkin(6));
t75 = qJ(2) * qJDD(1);
t91 = t43 * t23 + t44 * t75;
t8 = t27 + t91;
t6 = -qJDD(1) * pkin(5) + t8;
t24 = t49 * qJD(1) + qJD(2);
t80 = qJ(2) * qJD(1);
t12 = t43 * t24 + t44 * t80;
t10 = -qJD(1) * pkin(5) + t12;
t3 = t47 * qJD(3) - t45 * t10;
t84 = t3 * qJD(4);
t1 = t45 * qJDD(3) + t47 * t6 + t84;
t34 = t47 * qJDD(3);
t92 = t47 * t10;
t4 = t45 * qJD(3) + t92;
t2 = -t4 * qJD(4) - t45 * t6 + t34;
t53 = -(t3 * t47 + t4 * t45) * qJD(4) + t1 * t47 - t2 * t45;
t94 = t44 * t24;
t51 = qJD(1) ^ 2;
t93 = t44 * t51;
t19 = t44 * qJ(2) + t43 * t49;
t46 = sin(qJ(1));
t48 = cos(qJ(1));
t90 = t48 * pkin(1) + t46 * qJ(2);
t89 = g(1) * t46 - g(2) * t48;
t41 = t45 ^ 2;
t42 = t47 ^ 2;
t88 = t41 - t42;
t87 = t41 + t42;
t50 = qJD(4) ^ 2;
t86 = t50 + t51;
t83 = t43 * qJ(2);
t9 = -t94 + (pkin(3) + t83) * qJD(1);
t85 = qJD(1) * t9;
t82 = pkin(1) * qJDD(1);
t81 = qJDD(3) + g(3);
t79 = qJDD(4) * t45;
t78 = qJDD(4) * t47;
t77 = t45 * qJDD(1);
t76 = t47 * qJDD(1);
t73 = qJD(1) * qJD(4);
t72 = t45 * t51 * t47;
t71 = t48 * pkin(2) + t90;
t70 = 0.2e1 * t74;
t69 = t47 * t73;
t68 = t43 * t74;
t67 = t44 * t23 - t43 * t75;
t66 = qJDD(1) * t87;
t65 = qJDD(2) - t82;
t64 = t45 * t69;
t36 = t48 * qJ(2);
t63 = t49 * t46 + t36;
t13 = -t46 * t43 - t48 * t44;
t14 = t48 * t43 - t46 * t44;
t62 = g(1) * t14 - g(2) * t13;
t61 = -g(1) * t13 - g(2) * t14;
t60 = g(1) * t48 + g(2) * t46;
t57 = t3 * t45 - t4 * t47;
t56 = (-t43 * t80 + t94) * t43 - t12 * t44;
t18 = t44 * t49 - t83;
t7 = t67 - t68;
t55 = -qJD(3) * qJD(4) - t6 + t61 + t85;
t15 = pkin(3) - t18;
t16 = -pkin(5) + t19;
t54 = -qJDD(4) * t16 + (-qJD(1) * t15 - qJD(2) * t44 - t9) * qJD(4);
t5 = qJDD(1) * pkin(3) - t7;
t52 = qJDD(1) * t15 - t16 * t50 + t5 - t62 + t68;
t22 = -t50 * t45 + t78;
t21 = -t50 * t47 - t79;
t11 = [0, 0, 0, 0, 0, qJDD(1), t89, t60, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t82 + t89, 0, -t60 + t70 + 0.2e1 * t75, -t65 * pkin(1) - g(1) * (-t46 * pkin(1) + t36) - g(2) * t90 + (t70 + t75) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -t18 * qJDD(1) - t62 - t67 + 0.2e1 * t68, t19 * qJDD(1) + 0.2e1 * t27 - t61 + t91, 0, -g(1) * t63 - g(2) * t71 - t56 * qJD(2) + t7 * t18 + t8 * t19, t41 * qJDD(1) + 0.2e1 * t64, 0.2e1 * t45 * t76 - 0.2e1 * t88 * t73, t21, t42 * qJDD(1) - 0.2e1 * t64, -t22, 0, t54 * t45 + t52 * t47, -t52 * t45 + t54 * t47, -t16 * t66 - t87 * t27 - t53 + t61, t5 * t15 - g(1) * (t14 * pkin(3) + t13 * pkin(5) + t63) - g(2) * (-t13 * pkin(3) + t14 * pkin(5) + t71) + (t9 * t43 - t57 * t44) * qJD(2) + t53 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t51, -t51 * qJ(2) + t65 - t89, 0, 0, 0, 0, 0, 0, -t44 * qJDD(1) - t43 * t51, t43 * qJDD(1) - t93, 0, t56 * qJD(1) + t8 * t43 + t7 * t44 - t89, 0, 0, 0, 0, 0, 0, (0.2e1 * t45 * t73 - t76) * t44 + (-t86 * t47 - t79) * t43, (0.2e1 * t69 + t77) * t44 + (t86 * t45 - t78) * t43, -t43 * t66 + t87 * t93, (t57 * qJD(1) - t5) * t44 + (t53 - t85) * t43 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, 0, 0, 0, 0, 0, t22, t21, 0, -t57 * qJD(4) + t1 * t45 + t2 * t47 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t88 * t51, -t77, t72, -t76, qJDD(4), g(3) * t47 + t34 + (t4 - t92) * qJD(4) + t55 * t45, t84 + (qJD(4) * t10 - t81) * t45 + t55 * t47, 0, 0;];
tau_reg = t11;
