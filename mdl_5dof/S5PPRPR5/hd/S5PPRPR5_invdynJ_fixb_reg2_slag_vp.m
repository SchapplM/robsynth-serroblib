% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:34
% EndTime: 2019-12-31 17:33:34
% DurationCPUTime: 0.45s
% Computational Cost: add. (436->111), mult. (738->131), div. (0->0), fcn. (477->6), ass. (0->74)
t35 = -pkin(3) - pkin(6);
t32 = sin(qJ(3));
t79 = t32 * qJD(2);
t23 = qJD(3) * t79;
t34 = cos(qJ(3));
t71 = t34 * qJDD(2);
t51 = qJDD(4) + t23 - t71;
t6 = t35 * qJDD(3) + t51;
t93 = qJD(5) * qJD(1) + t6;
t37 = qJD(3) ^ 2;
t92 = t34 * qJDD(3) - t37 * t32;
t77 = qJD(3) * qJ(4);
t19 = t77 + t79;
t91 = (t19 + t77 - t79) * qJD(5) + qJDD(5) * t35;
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t78 = t34 * qJD(2);
t61 = qJD(4) - t78;
t44 = -t35 * qJD(3) - t61;
t4 = t31 * qJD(1) - t33 * t44;
t84 = t4 * qJD(5);
t1 = -t33 * qJDD(1) + t31 * t6 + t84;
t41 = qJD(5) * t44;
t66 = t31 * qJDD(1) + t93 * t33;
t2 = t31 * t41 + t66;
t42 = t31 * t44;
t5 = t33 * qJD(1) + t42;
t39 = (t31 * t4 + t33 * t5) * qJD(5) - t1 * t31 - t2 * t33;
t90 = t19 * t34;
t89 = t31 * t33;
t28 = t31 ^ 2;
t29 = t33 ^ 2;
t87 = t28 - t29;
t86 = t28 + t29;
t36 = qJD(5) ^ 2;
t85 = t36 + t37;
t83 = cos(pkin(7));
t82 = sin(pkin(7));
t81 = qJDD(3) * pkin(3);
t80 = t19 * qJD(3);
t30 = qJDD(1) - g(3);
t76 = qJDD(3) * t32;
t75 = qJDD(5) * t31;
t73 = t32 * qJDD(2);
t72 = t33 * qJDD(3);
t69 = qJD(3) * qJD(5);
t67 = qJDD(3) * qJ(4);
t65 = t37 * t89;
t12 = -t82 * t32 - t83 * t34;
t13 = t83 * t32 - t82 * t34;
t64 = t12 * pkin(3) - t13 * qJ(4);
t63 = -t13 * pkin(3) - t12 * qJ(4);
t62 = t86 * qJDD(3);
t60 = t69 * t89;
t59 = g(1) * t13 - g(2) * t12;
t58 = g(1) * t12 + g(2) * t13;
t56 = t31 * t5 - t33 * t4;
t54 = -g(1) * t82 + g(2) * t83;
t53 = (-qJD(3) * pkin(3) + t61) * t32 + t90;
t7 = t67 + t73 + (qJD(4) + t78) * qJD(3);
t52 = t7 * qJ(4) + t19 * qJD(4);
t49 = t85 * t34 + t76;
t48 = -qJDD(5) * t34 + 0.2e1 * t32 * t69;
t47 = t59 + t80;
t46 = -t58 - t73;
t45 = t59 + t71;
t43 = qJDD(4) - t45;
t38 = t61 * qJD(3) - t35 * t36 + t58 + t67 + t7;
t26 = qJDD(5) * t33;
t16 = t37 * t34 + t76;
t15 = -t36 * t31 + t26;
t14 = t36 * t33 + t75;
t8 = t51 - t81;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, t14, t15, 0, -t56 * qJD(5) - t1 * t33 + t2 * t31 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t54, 0, 0, 0, 0, 0, 0, t92, -t16, 0, (t32 ^ 2 + t34 ^ 2) * qJDD(2) + t54, 0, 0, 0, 0, 0, 0, 0, -t92, t16, t53 * qJD(3) + t7 * t32 - t8 * t34 + t54, 0, 0, 0, 0, 0, 0, t49 * t31 + t48 * t33, -t48 * t31 + t49 * t33, t92 * t86, (-t56 * qJD(3) + t7) * t32 + (t39 + t80) * t34 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t45, t46, 0, 0, qJDD(3), 0, 0, 0, 0, 0, 0, t43 - 0.2e1 * t81, 0.2e1 * qJD(3) * qJD(4) - t46 + 0.2e1 * t67, -t8 * pkin(3) - g(1) * t63 - g(2) * t64 - t53 * qJD(2) + t52, t29 * qJDD(3) - 0.2e1 * t60, -0.2e1 * t31 * t72 + 0.2e1 * t87 * t69, t15, t28 * qJDD(3) + 0.2e1 * t60, -t14, 0, t38 * t31 + t91 * t33, -t91 * t31 + t38 * t33, t86 * t23 - t35 * t62 + t39 + t59, -g(1) * (-t13 * pkin(6) + t63) - g(2) * (t12 * pkin(6) + t64) + (t56 * t32 - t90) * qJD(2) - t39 * t35 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t37, t23 + t43 - t80 - t81, 0, 0, 0, 0, 0, 0, -t85 * t31 + t26, -t85 * t33 - t75, -t62, -t39 - t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t87 * t37, t72, -t65, -t31 * qJDD(3), qJDD(5), -g(3) * t31 + (t42 - t5) * qJD(5) - t47 * t33 + t66, t84 + (t41 + t30) * t33 + (t47 - t93) * t31, 0, 0;];
tau_reg = t3;
