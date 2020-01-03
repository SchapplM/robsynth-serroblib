% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:37
% EndTime: 2019-12-31 17:38:39
% DurationCPUTime: 0.42s
% Computational Cost: add. (596->92), mult. (1163->142), div. (0->0), fcn. (716->6), ass. (0->70)
t34 = sin(pkin(8));
t35 = cos(pkin(8));
t37 = sin(qJ(2));
t39 = cos(qJ(2));
t17 = t37 * t34 + t39 * t35;
t13 = t17 * qJD(1);
t86 = -t35 * qJD(3) + t13;
t89 = qJD(2) * t86;
t62 = t39 * qJD(1);
t63 = t37 * qJD(1);
t12 = t34 * t62 - t35 * t63;
t23 = (qJD(3) + t62) * qJD(2);
t57 = qJD(2) * t63;
t9 = t34 * t23 - t35 * t57;
t88 = (t34 * qJD(3) - t12) * qJD(2) + t9;
t41 = qJD(5) ^ 2;
t42 = qJD(2) ^ 2;
t87 = t34 * (t41 + t42);
t40 = -pkin(2) - pkin(3);
t73 = t35 * qJ(3) + t34 * t40;
t21 = -pkin(6) + t73;
t85 = -t21 * t41 + t88;
t84 = t39 * t34 - t37 * t35;
t10 = t35 * t23 + t34 * t57;
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t52 = qJD(3) - t62;
t19 = t40 * qJD(2) + t52;
t27 = qJD(2) * qJ(3) + t63;
t8 = t34 * t19 + t35 * t27;
t6 = -qJD(2) * pkin(6) + t8;
t3 = t38 * qJD(4) - t36 * t6;
t1 = t3 * qJD(5) + t38 * t10;
t4 = t36 * qJD(4) + t38 * t6;
t2 = -t4 * qJD(5) - t36 * t10;
t43 = -(t3 * t38 + t36 * t4) * qJD(5) + t1 * t38 - t2 * t36;
t83 = t9 * t17;
t82 = t27 * t39;
t81 = t35 * t42;
t80 = t36 * t38;
t77 = t41 * t36;
t76 = t41 * t38;
t75 = t42 * t37;
t74 = t42 * t39;
t32 = t36 ^ 2;
t33 = t38 ^ 2;
t72 = t32 - t33;
t71 = t32 + t33;
t69 = qJD(2) * pkin(2);
t7 = t35 * t19 - t34 * t27;
t5 = qJD(2) * pkin(4) - t7;
t68 = qJD(2) * t5;
t67 = qJD(5) * t36;
t66 = qJD(5) * t38;
t15 = t17 * qJD(2);
t65 = t15 * qJD(2);
t61 = qJD(2) * qJD(5);
t60 = t42 * t80;
t59 = 0.2e1 * t61;
t56 = -t10 + t68;
t53 = t35 * t59;
t51 = t61 * t80;
t48 = t3 * t36 - t38 * t4;
t47 = t7 * t34 - t8 * t35;
t46 = -t34 * qJ(3) + t35 * t40;
t20 = pkin(4) - t46;
t44 = qJD(5) * (-qJD(2) * t20 - t5 + t86);
t24 = t52 - t69;
t14 = t84 * qJD(2);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, 0, 0, 0, 0, 0, 0, 0, -t75, 0, t74, t23 * t37 + (t82 + (t24 - t62) * t37) * qJD(2), 0, 0, 0, 0, 0, 0, t14 * qJD(2), t65, 0, -t10 * t84 - t7 * t14 + t8 * t15 + t83, 0, 0, 0, 0, 0, 0, -t15 * t67 + t84 * t76 + (t14 * t38 - t17 * t67) * qJD(2), -t15 * t66 - t84 * t77 + (-t14 * t36 - t17 * t66) * qJD(2), -t71 * t65, t5 * t14 - t48 * t15 - t43 * t84 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t23 * qJ(3) + t27 * qJD(3) + (-t82 + (-t24 - t69) * t37) * qJD(1), 0, 0, 0, 0, 0, 0, t88, t10 - t89, 0, -t47 * qJD(3) + t10 * t73 + t7 * t12 - t8 * t13 - t9 * t46, 0.2e1 * t51, -t72 * t59, -t76, -0.2e1 * t51, t77, 0, t36 * t44 + t85 * t38, -t85 * t36 + t38 * t44, t71 * t89 - t43, -t5 * t12 + t9 * t20 + t48 * t13 + (t34 * t5 - t48 * t35) * qJD(3) + t43 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, (-t27 + t63) * qJD(2), 0, 0, 0, 0, 0, 0, -t34 * t42, -t81, 0, t47 * qJD(2) + t10 * t34 - t9 * t35, 0, 0, 0, 0, 0, 0, t36 * t53 - t38 * t87, t36 * t87 + t38 * t53, t71 * t81, (t48 * qJD(2) - t9) * t35 + (t43 - t68) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, -t48 * qJD(5) + t1 * t36 + t2 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t72 * t42, 0, t60, 0, 0, t56 * t36, t56 * t38, 0, 0;];
tauc_reg = t11;
