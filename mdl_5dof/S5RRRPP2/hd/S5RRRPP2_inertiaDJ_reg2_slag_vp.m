% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:59
% EndTime: 2019-12-31 20:52:01
% DurationCPUTime: 0.75s
% Computational Cost: add. (385->108), mult. (970->149), div. (0->0), fcn. (574->4), ass. (0->69)
t62 = sin(qJ(3));
t60 = t62 ^ 2;
t64 = cos(qJ(3));
t61 = t64 ^ 2;
t65 = cos(qJ(2));
t83 = pkin(1) * qJD(2);
t75 = t65 * t83;
t95 = (t60 + t61) * t75;
t56 = t62 * qJ(4);
t93 = t64 * pkin(3) + t56;
t92 = 0.2e1 * (-t60 + t61) * qJD(3);
t67 = 2 * qJD(4);
t50 = -t65 * pkin(1) - pkin(2);
t25 = t50 - t93;
t57 = t64 * pkin(4);
t17 = -t25 + t57;
t63 = sin(qJ(2));
t76 = t63 * t83;
t54 = t62 * qJD(3);
t55 = t64 * qJD(3);
t23 = pkin(3) * t54 - qJ(4) * t55 - t62 * qJD(4);
t9 = -pkin(4) * t54 - t23;
t5 = t9 - t76;
t91 = t17 * t55 + t5 * t62;
t79 = pkin(2) + t93;
t26 = t57 + t79;
t90 = t26 * t55 + t9 * t62;
t89 = pkin(7) - qJ(5);
t10 = t76 + t23;
t88 = -t10 - t23;
t49 = t63 * pkin(1) + pkin(7);
t87 = t95 * t49;
t86 = t95 * pkin(7);
t85 = t50 * t55 + t62 * t76;
t80 = qJ(5) * qJD(3);
t84 = -t62 * qJD(5) - t64 * t80;
t82 = -qJ(5) + t49;
t81 = t64 * qJD(4);
t78 = pkin(2) * t55;
t77 = pkin(2) * t54;
t74 = pkin(7) * t54;
t73 = pkin(7) * t55;
t72 = t62 * t55;
t71 = -t64 * qJD(5) + t62 * t80;
t12 = t49 * t54 - t64 * t75;
t68 = t50 * t54 - t64 * t76;
t13 = t49 * t55 + t62 * t75;
t20 = -t93 * qJD(3) + t81;
t66 = -pkin(3) - pkin(4);
t59 = qJ(4) * t67;
t45 = -0.2e1 * t72;
t44 = 0.2e1 * t72;
t39 = t89 * t64;
t38 = t89 * t62;
t29 = t39 * t54;
t28 = t82 * t64;
t27 = t82 * t62;
t24 = t79 * t54;
t22 = t73 + t84;
t21 = t71 - t74;
t18 = -t81 + (-t64 * t66 + t56) * qJD(3);
t16 = t28 * t54;
t14 = t25 * t54;
t11 = 0.2e1 * t95;
t7 = t9 * t64;
t4 = t13 + t84;
t3 = -t12 + t71;
t2 = t5 * t64;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76, -0.2e1 * t75, 0, 0, t44, t92, 0, t45, 0, 0, 0.2e1 * t68, 0.2e1 * t85, t11, 0.2e1 * t50 * t76 + 0.2e1 * t87, t44, 0, -t92, 0, 0, t45, -0.2e1 * t10 * t64 + 0.2e1 * t14, t11, -0.2e1 * t10 * t62 - 0.2e1 * t25 * t55, 0.2e1 * t25 * t10 + 0.2e1 * t87, t44, -t92, 0, t45, 0, 0, -0.2e1 * t17 * t54 + 0.2e1 * t2, 0.2e1 * t91, -0.2e1 * t4 * t62 + 0.2e1 * t16 + 0.2e1 * (-qJD(3) * t27 - t3) * t64, 0.2e1 * t17 * t5 + 0.2e1 * t27 * t4 + 0.2e1 * t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0, t44, t92, 0, t45, 0, 0, t68 - t77, -t78 + t85, t95, -pkin(2) * t76 + t86, t44, 0, -t92, 0, 0, t45, t88 * t64 + t14 - t24, t95, t88 * t62 + (-t25 + t79) * t55, -t10 * t79 + t25 * t23 + t86, t44, -t92, 0, t45, 0, 0, t2 + t7 + (-t17 - t26) * t54, t90 + t91, t16 + t29 + (-t22 - t4) * t62 + (-t21 - t3 + (-t27 - t38) * qJD(3)) * t64, t17 * t9 + t28 * t21 + t27 * t22 + t5 * t26 + t3 * t39 + t4 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t92, 0, t45, 0, 0, -0.2e1 * t77, -0.2e1 * t78, 0, 0, t44, 0, -t92, 0, 0, t45, -0.2e1 * t23 * t64 - 0.2e1 * t24, 0, -0.2e1 * t23 * t62 + 0.2e1 * t55 * t79, -0.2e1 * t79 * t23, t44, -t92, 0, t45, 0, 0, -0.2e1 * t26 * t54 + 0.2e1 * t7, 0.2e1 * t90, -0.2e1 * t22 * t62 + 0.2e1 * t29 + 0.2e1 * (-qJD(3) * t38 - t21) * t64, 0.2e1 * t39 * t21 + 0.2e1 * t38 * t22 + 0.2e1 * t26 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t54, 0, -t13, t12, 0, 0, 0, t55, 0, 0, t54, 0, -t13, t20, -t12, (-pkin(3) * t62 + qJ(4) * t64) * t75 + t20 * t49, 0, 0, -t55, 0, -t54, 0, -t4, t3, t18, t3 * qJ(4) + t28 * qJD(4) + t4 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t54, 0, -t73, t74, 0, 0, 0, t55, 0, 0, t54, 0, -t73, t20, -t74, t20 * pkin(7), 0, 0, -t55, 0, -t54, 0, -t22, t21, t18, t21 * qJ(4) + t39 * qJD(4) + t22 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t59, 0, 0, 0, 0, 0, 0, 0, t67, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t55, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t55, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
