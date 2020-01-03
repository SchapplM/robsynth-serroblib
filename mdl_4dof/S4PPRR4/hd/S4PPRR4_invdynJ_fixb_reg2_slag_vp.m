% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:43
% EndTime: 2019-12-31 16:18:43
% DurationCPUTime: 0.37s
% Computational Cost: add. (383->87), mult. (909->123), div. (0->0), fcn. (723->10), ass. (0->71)
t30 = sin(pkin(7));
t32 = cos(pkin(7));
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t14 = t37 * t30 + t35 * t32;
t12 = t14 * qJD(3);
t65 = qJDD(1) * t37;
t66 = qJDD(1) * t35;
t54 = t30 * t66 - t32 * t65;
t68 = qJDD(3) * pkin(3);
t4 = qJD(1) * t12 + t54 - t68;
t27 = pkin(7) + qJ(3);
t23 = sin(t27);
t24 = cos(t27);
t31 = sin(pkin(6));
t33 = cos(pkin(6));
t56 = g(1) * t33 + g(2) * t31;
t45 = -g(3) * t24 + t56 * t23;
t84 = -t4 + t45;
t38 = qJD(4) ^ 2;
t81 = t14 * qJD(1);
t82 = qJD(3) * t81;
t83 = -pkin(5) * t38 + t68 + t82 + t84;
t59 = -g(1) * t31 + g(2) * t33;
t75 = t37 * t32;
t76 = t35 * t30;
t13 = -t75 + t76;
t9 = t13 * qJD(1);
t47 = g(3) * t23 + t56 * t24;
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t77 = t34 * t36;
t28 = t34 ^ 2;
t29 = t36 ^ 2;
t74 = t28 - t29;
t73 = t28 + t29;
t72 = qJD(3) * pkin(3);
t8 = qJD(3) * pkin(5) + t81;
t5 = t36 * qJD(2) - t34 * t8;
t71 = t5 * qJD(4);
t6 = t34 * qJD(2) + t36 * t8;
t70 = t6 * qJD(4);
t69 = qJD(3) * t13;
t67 = t69 * qJD(3);
t64 = t36 * qJDD(3);
t63 = qJD(3) * qJD(4);
t39 = qJD(3) ^ 2;
t62 = t39 * t77;
t61 = -qJD(3) * qJD(1) * t75 - t30 * t65 - t32 * t66;
t60 = qJD(1) * t76;
t58 = t73 * qJDD(3);
t57 = t63 * t77;
t55 = t5 * t34 - t6 * t36;
t53 = -t12 * qJD(3) - t13 * qJDD(3);
t50 = -qJD(4) * t8 + t59;
t49 = -qJD(3) * t60 - t61;
t48 = t14 * t38 - t53;
t46 = 0.2e1 * t69 * qJD(4) - qJDD(4) * t14;
t7 = -t72 + t9;
t44 = -pkin(5) * qJDD(4) + (t7 - t9 - t72) * qJD(4);
t3 = qJDD(3) * pkin(5) + t49;
t1 = t34 * qJDD(2) + t36 * t3 + t71;
t25 = t36 * qJDD(2);
t2 = -t34 * t3 + t25 - t70;
t42 = t1 * t36 - t2 * t34 + (-t34 * t6 - t36 * t5) * qJD(4);
t41 = -qJD(4) * qJD(2) - t7 * qJD(3) - t3 + t47;
t40 = t42 - t47;
t17 = qJDD(4) * t36 - t38 * t34;
t16 = qJDD(4) * t34 + t38 * t36;
t15 = qJDD(2) + t59;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t30 ^ 2 + t32 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t53, -t14 * qJDD(3) + t67, 0, t49 * t14 - t81 * t69 - (-t54 - t82) * t13 + t9 * t12 - g(3), 0, 0, 0, 0, 0, 0, t46 * t34 - t48 * t36, t48 * t34 + t46 * t36, t14 * t58 - t73 * t67, t7 * t12 + t4 * t13 + t42 * t14 + t55 * t69 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t55 * qJD(4) + t1 * t34 + t2 * t36 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t45 - t54, (-t9 + t60) * qJD(3) + t47 + t61, 0, 0, t28 * qJDD(3) + 0.2e1 * t57, 0.2e1 * t34 * t64 - 0.2e1 * t74 * t63, t16, t29 * qJDD(3) - 0.2e1 * t57, t17, 0, t44 * t34 + t83 * t36, -t83 * t34 + t44 * t36, t73 * t9 * qJD(3) + pkin(5) * t58 + t40, t84 * pkin(3) + t40 * pkin(5) - t55 * t9 - t7 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t74 * t39, t34 * qJDD(3), t62, t64, qJDD(4), t41 * t34 + t50 * t36 + t25 + t70, t71 + (-qJDD(2) - t50) * t34 + t41 * t36, 0, 0;];
tau_reg = t10;
