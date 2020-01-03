% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (377->93), mult. (925->143), div. (0->0), fcn. (668->6), ass. (0->74)
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t95 = t70 * pkin(3) + t67 * qJ(4);
t94 = 2 * qJD(4);
t93 = pkin(7) - pkin(8);
t68 = sin(qJ(2));
t54 = pkin(1) * t68 + pkin(7);
t92 = -pkin(8) + t54;
t60 = t70 * qJD(3);
t66 = sin(qJ(5));
t69 = cos(qJ(5));
t81 = t67 * qJD(3);
t82 = qJD(5) * t70;
t83 = qJD(5) * t69;
t15 = t67 * t83 - t69 * t81 + (t60 - t82) * t66;
t30 = pkin(3) * t81 - qJ(4) * t60 - t67 * qJD(4);
t20 = -pkin(4) * t81 - t30;
t85 = pkin(1) * qJD(2);
t77 = t68 * t85;
t17 = t20 - t77;
t71 = cos(qJ(2));
t55 = -pkin(1) * t71 - pkin(2);
t32 = t55 - t95;
t62 = t70 * pkin(4);
t25 = -t32 + t62;
t37 = t66 * t67 + t69 * t70;
t91 = t25 * t15 + t17 * t37;
t84 = qJD(5) * t66;
t16 = qJD(3) * t37 - t67 * t84 - t69 * t82;
t38 = -t66 * t70 + t67 * t69;
t90 = t25 * t16 + t17 * t38;
t80 = pkin(2) + t95;
t33 = t62 + t80;
t89 = t33 * t15 + t20 * t37;
t88 = t33 * t16 + t20 * t38;
t21 = t77 + t30;
t87 = -t21 - t30;
t86 = t55 * t60 + t67 * t77;
t79 = pkin(2) * t81;
t78 = pkin(2) * t60;
t76 = t71 * t85;
t75 = pkin(7) * t81;
t35 = t92 * t70;
t74 = t67 * t76;
t22 = t54 * t81 - t70 * t76;
t64 = t67 ^ 2;
t65 = t70 ^ 2;
t28 = (t64 + t65) * t76;
t73 = t55 * t81 - t70 * t77;
t29 = -qJD(3) * t95 + t70 * qJD(4);
t72 = -pkin(3) - pkin(4);
t57 = pkin(8) * t81;
t56 = pkin(7) * t60;
t48 = 0.2e1 * t67 * t60;
t45 = t93 * t70;
t44 = t93 * t67;
t40 = -pkin(8) * t60 + t56;
t39 = t57 - t75;
t36 = 0.2e1 * (-t64 + t65) * qJD(3);
t34 = t92 * t67;
t31 = t80 * t81;
t27 = t66 * qJD(4) + (qJ(4) * t69 + t66 * t72) * qJD(5);
t26 = -t69 * qJD(4) + (qJ(4) * t66 - t69 * t72) * qJD(5);
t24 = t32 * t81;
t23 = t54 * t60 + t74;
t19 = qJD(3) * t35 + t74;
t18 = -t22 + t57;
t10 = 0.2e1 * t38 * t16;
t5 = t66 * t39 - t69 * t40 + (t44 * t66 + t45 * t69) * qJD(5);
t4 = -t69 * t39 - t66 * t40 + (-t44 * t69 + t45 * t66) * qJD(5);
t3 = -0.2e1 * t15 * t38 - 0.2e1 * t16 * t37;
t2 = t66 * t18 - t69 * t19 + (t34 * t66 + t35 * t69) * qJD(5);
t1 = -t69 * t18 - t66 * t19 + (-t34 * t69 + t35 * t66) * qJD(5);
t6 = [0, 0, 0, 0, -0.2e1 * t77, -0.2e1 * t76, t48, t36, 0, 0, 0, 0.2e1 * t73, 0.2e1 * t86, -0.2e1 * t21 * t70 + 0.2e1 * t24, 0.2e1 * t28, -0.2e1 * t21 * t67 - 0.2e1 * t32 * t60, 0.2e1 * t32 * t21 + 0.2e1 * t28 * t54, t10, t3, 0, 0, 0, 0.2e1 * t91, 0.2e1 * t90; 0, 0, 0, 0, -t77, -t76, t48, t36, 0, 0, 0, t73 - t79, -t78 + t86, t70 * t87 + t24 - t31, t28, t87 * t67 + (-t32 + t80) * t60, pkin(7) * t28 - t21 * t80 + t32 * t30, t10, t3, 0, 0, 0, t89 + t91, t88 + t90; 0, 0, 0, 0, 0, 0, t48, t36, 0, 0, 0, -0.2e1 * t79, -0.2e1 * t78, -0.2e1 * t30 * t70 - 0.2e1 * t31, 0, -0.2e1 * t30 * t67 + 0.2e1 * t60 * t80, -0.2e1 * t80 * t30, t10, t3, 0, 0, 0, 0.2e1 * t89, 0.2e1 * t88; 0, 0, 0, 0, 0, 0, 0, 0, t60, -t81, 0, -t23, t22, -t23, t29, -t22, (-pkin(3) * t67 + qJ(4) * t70) * t76 + t29 * t54, 0, 0, -t16, t15, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, t60, -t81, 0, -t56, t75, -t56, t29, -t75, t29 * pkin(7), 0, 0, -t16, t15, 0, t5, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, qJ(4) * t94, 0, 0, 0, 0, 0, 0.2e1 * t27, -0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, -t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
