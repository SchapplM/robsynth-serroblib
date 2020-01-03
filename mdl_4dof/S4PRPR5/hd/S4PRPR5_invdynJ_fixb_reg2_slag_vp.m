% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPR5
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:21
% EndTime: 2019-12-31 16:23:22
% DurationCPUTime: 0.45s
% Computational Cost: add. (480->110), mult. (1021->158), div. (0->0), fcn. (762->10), ass. (0->76)
t36 = qJ(2) + pkin(7);
t31 = sin(t36);
t32 = cos(t36);
t40 = sin(pkin(6));
t42 = cos(pkin(6));
t63 = g(1) * t42 + g(2) * t40;
t92 = -g(3) * t32 + t63 * t31;
t39 = sin(pkin(7));
t41 = cos(pkin(7));
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t21 = t39 * t46 + t41 * t44;
t16 = t21 * qJD(1);
t29 = t39 * pkin(2) + pkin(5);
t34 = t46 * qJDD(1);
t70 = qJD(1) * qJD(2);
t77 = qJDD(2) * pkin(2);
t19 = -t44 * t70 + t34 + t77;
t14 = t41 * t19;
t72 = t44 * qJDD(1);
t90 = t46 * t70 + t72;
t5 = -t90 * t39 + t14;
t3 = -qJDD(2) * pkin(3) - t5;
t30 = -t41 * pkin(2) - pkin(3);
t47 = qJD(4) ^ 2;
t91 = t16 * qJD(2) - qJDD(2) * t30 - t29 * t47 - t3 + t92;
t67 = -g(1) * t40 + g(2) * t42;
t55 = g(3) * t31 + t63 * t32;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t84 = t43 * t45;
t37 = t43 ^ 2;
t38 = t45 ^ 2;
t83 = t37 - t38;
t82 = t37 + t38;
t74 = t46 * qJD(1);
t26 = qJD(2) * pkin(2) + t74;
t79 = qJD(1) * t44;
t12 = t39 * t26 + t41 * t79;
t10 = qJD(2) * pkin(5) + t12;
t7 = t45 * qJD(3) - t43 * t10;
t81 = t7 * qJD(4);
t8 = t43 * qJD(3) + t45 * t10;
t80 = t8 * qJD(4);
t20 = t39 * t44 - t41 * t46;
t78 = qJD(2) * t20;
t76 = t78 * qJD(2);
t18 = t20 * qJD(1);
t75 = t18 * qJD(2);
t73 = qJDD(1) - g(3);
t71 = t45 * qJDD(2);
t69 = qJD(2) * qJD(4);
t48 = qJD(2) ^ 2;
t68 = t48 * t84;
t6 = t39 * t19 + t90 * t41;
t65 = qJDD(2) * t82;
t64 = t69 * t84;
t62 = t7 * t43 - t8 * t45;
t15 = t21 * qJD(2);
t61 = -t15 * qJD(2) - t20 * qJDD(2);
t11 = t41 * t26 - t39 * t79;
t57 = -qJD(4) * t10 + t67;
t56 = t21 * t47 - t61;
t54 = 0.2e1 * t78 * qJD(4) - qJDD(4) * t21;
t53 = -g(3) * t46 + t63 * t44;
t9 = -qJD(2) * pkin(3) - t11;
t52 = -qJDD(4) * t29 + (qJD(2) * t30 - t18 + t9) * qJD(4);
t4 = qJDD(2) * pkin(5) + t6;
t1 = t43 * qJDD(3) + t45 * t4 + t81;
t33 = t45 * qJDD(3);
t2 = -t43 * t4 + t33 - t80;
t50 = t1 * t45 - t2 * t43 + (-t43 * t8 - t45 * t7) * qJD(4);
t49 = -t9 * qJD(2) - qJD(4) * qJD(3) - t4 + t55;
t24 = qJDD(4) * t45 - t47 * t43;
t23 = qJDD(4) * t43 + t47 * t45;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, t46 * qJDD(2) - t48 * t44, -qJDD(2) * t44 - t48 * t46, 0, -g(3) + (t44 ^ 2 + t46 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t61, -t21 * qJDD(2) + t76, 0, -t11 * t15 - t12 * t78 - t5 * t20 + t6 * t21 - g(3), 0, 0, 0, 0, 0, 0, t54 * t43 - t56 * t45, t56 * t43 + t54 * t45, t21 * t65 - t82 * t76, t9 * t15 + t3 * t20 + t50 * t21 + t62 * t78 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t34 + t53, -t73 * t44 + t63 * t46, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t41 * t77 - t39 * t72 + t14 + (-t39 * t74 + t16) * qJD(2) + t92, -t39 * t77 + t55 - t6 - t75, 0, t11 * t16 + t12 * t18 + (t39 * t6 + t41 * t5 + t53) * pkin(2), t37 * qJDD(2) + 0.2e1 * t64, 0.2e1 * t43 * t71 - 0.2e1 * t83 * t69, t23, t38 * qJDD(2) - 0.2e1 * t64, t24, 0, t52 * t43 + t91 * t45, -t91 * t43 + t52 * t45, t29 * t65 + t82 * t75 + t50 - t55, t3 * t30 - t9 * t16 - g(3) * (t46 * pkin(2) + t32 * pkin(3) + t31 * pkin(5)) - t62 * t18 + t50 * t29 + t63 * (pkin(2) * t44 + pkin(3) * t31 - pkin(5) * t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + t67, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t62 * qJD(4) + t1 * t43 + t2 * t45 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t83 * t48, t43 * qJDD(2), t68, t71, qJDD(4), t49 * t43 + t57 * t45 + t33 + t80, t81 + (-qJDD(3) - t57) * t43 + t49 * t45, 0, 0;];
tau_reg = t13;
