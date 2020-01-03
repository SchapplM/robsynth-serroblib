% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:57
% EndTime: 2019-12-31 19:55:00
% DurationCPUTime: 0.76s
% Computational Cost: add. (1421->108), mult. (3213->192), div. (0->0), fcn. (2973->6), ass. (0->67)
t52 = sin(qJ(2));
t53 = cos(qJ(2));
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t62 = t84 * t52 - t85 * t53;
t40 = t85 * t52 + t84 * t53;
t54 = 2 * qJD(5);
t49 = t85 * pkin(2) + pkin(3);
t77 = t84 * pkin(2);
t89 = sin(qJ(4));
t90 = cos(qJ(4));
t38 = t89 * t49 + t90 * t77;
t30 = t38 * qJD(4);
t87 = -qJ(3) - pkin(6);
t24 = t62 * t87;
t20 = -t62 * pkin(7) + t24;
t23 = t40 * t87;
t59 = -t40 * pkin(7) + t23;
t58 = t90 * t59;
t9 = t89 * t20 - t58;
t91 = t9 * t30;
t22 = t90 * t40 - t89 * t62;
t88 = t30 * t22;
t68 = qJD(2) * t84;
t69 = qJD(2) * t85;
t86 = t52 * t69 + t53 * t68;
t83 = t52 * qJD(2);
t82 = t53 * qJD(2);
t65 = -t52 * t68 + t53 * t69;
t12 = t22 * qJD(4) + t89 * t65 + t90 * t86;
t60 = t90 * t62;
t21 = t89 * t40 + t60;
t81 = 0.2e1 * t21 * t12;
t80 = -0.2e1 * pkin(1) * qJD(2);
t51 = pkin(2) * t83;
t57 = t89 * t59;
t10 = t90 * t20 + t57;
t70 = qJD(2) * t87;
t63 = t53 * qJD(3) + t52 * t70;
t64 = -t52 * qJD(3) + t53 * t70;
t18 = -t84 * t63 + t85 * t64;
t55 = t65 * pkin(7) - t18;
t19 = t85 * t63 + t84 * t64;
t56 = -t86 * pkin(7) + t19;
t75 = qJD(4) * t89;
t3 = -qJD(4) * t58 + t20 * t75 + t89 * t55 - t90 * t56;
t76 = qJD(4) * t90;
t4 = qJD(4) * t57 + t20 * t76 + t90 * t55 + t89 * t56;
t79 = -t10 * t3 + t9 * t4;
t78 = t52 * t82;
t50 = -t53 * pkin(2) - pkin(1);
t25 = t86 * pkin(3) + t51;
t11 = qJD(4) * t60 + t40 * t75 - t90 * t65 + t89 * t86;
t67 = t11 * t21 - t22 * t12;
t66 = t89 * t77;
t61 = -0.2e1 * t10 * t12 - 0.2e1 * t9 * t11 + 0.2e1 * t3 * t21 + 0.2e1 * t4 * t22;
t29 = qJD(4) * t66 - t49 * t76;
t37 = t90 * t49 - t66;
t28 = t62 * pkin(3) + t50;
t35 = -pkin(4) - t37;
t34 = qJ(5) + t38;
t27 = qJD(5) - t29;
t26 = 0.2e1 * t30;
t8 = t21 * pkin(4) - t22 * qJ(5) + t28;
t7 = -0.2e1 * t22 * t11;
t5 = t12 * pkin(4) + t11 * qJ(5) - t22 * qJD(5) + t25;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, 0.2e1 * (-t52 ^ 2 + t53 ^ 2) * qJD(2), 0, -0.2e1 * t78, 0, 0, t52 * t80, t53 * t80, 0, 0, 0.2e1 * t40 * t65, -0.2e1 * t40 * t86 - 0.2e1 * t65 * t62, 0, 0.2e1 * t62 * t86, 0, 0, 0.2e1 * t50 * t86 + 0.2e1 * t62 * t51, 0.2e1 * t40 * t51 + 0.2e1 * t50 * t65, -0.2e1 * t18 * t40 - 0.2e1 * t19 * t62 - 0.2e1 * t23 * t65 - 0.2e1 * t24 * t86, 0.2e1 * t23 * t18 + 0.2e1 * t24 * t19 + 0.2e1 * t50 * t51, t7, 0.2e1 * t67, 0, t81, 0, 0, 0.2e1 * t28 * t12 + 0.2e1 * t25 * t21, -0.2e1 * t28 * t11 + 0.2e1 * t25 * t22, t61, 0.2e1 * t28 * t25 + 0.2e1 * t79, t7, 0, -0.2e1 * t67, 0, 0, t81, 0.2e1 * t8 * t12 + 0.2e1 * t5 * t21, t61, 0.2e1 * t8 * t11 - 0.2e1 * t5 * t22, 0.2e1 * t8 * t5 + 0.2e1 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, -t83, 0, -pkin(6) * t82, pkin(6) * t83, 0, 0, 0, 0, t65, 0, -t86, 0, t18, -t19, (-t85 * t65 - t84 * t86) * pkin(2), (t85 * t18 + t84 * t19) * pkin(2), 0, 0, -t11, 0, -t12, 0, -t4, t3, t37 * t11 - t38 * t12 + t29 * t21 + t88, -t10 * t29 - t3 * t38 - t4 * t37 + t91, 0, -t11, 0, 0, t12, 0, -t4, -t35 * t11 - t34 * t12 - t27 * t21 + t88, -t3, t10 * t27 - t3 * t34 + t4 * t35 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0.2e1 * t29, 0, -0.2e1 * t38 * t29 - 0.2e1 * t37 * t30, 0, 0, 0, 0, 0, 0, -t26, 0, 0.2e1 * t27, 0.2e1 * t34 * t27 + 0.2e1 * t35 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t65, 0, t51, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t25, 0, 0, 0, 0, 0, 0, t12, 0, t11, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, -t4, t3, 0, 0, 0, -t11, 0, 0, t12, 0, -t4, pkin(4) * t11 - t12 * qJ(5) - t21 * qJD(5), -t3, -t4 * pkin(4) - t3 * qJ(5) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, t54 - t29, -t30 * pkin(4) + t27 * qJ(5) + t34 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, qJ(5) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
