% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP9
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
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:31
% DurationCPUTime: 0.87s
% Computational Cost: add. (1186->158), mult. (2967->300), div. (0->0), fcn. (2537->6), ass. (0->79)
t102 = cos(qJ(4));
t65 = sin(qJ(2));
t66 = cos(qJ(2));
t73 = -t66 * pkin(2) - t65 * qJ(3);
t45 = -pkin(1) + t73;
t63 = cos(pkin(8));
t41 = t63 * t45;
t62 = sin(pkin(8));
t99 = t63 * t65;
t22 = -pkin(7) * t99 + t41 + (-pkin(6) * t62 - pkin(3)) * t66;
t101 = t62 * t65;
t98 = t63 * t66;
t55 = pkin(6) * t98;
t30 = t62 * t45 + t55;
t26 = -pkin(7) * t101 + t30;
t64 = sin(qJ(4));
t107 = t102 * t26 + t64 * t22;
t77 = qJD(4) * t102;
t93 = qJD(4) * t64;
t106 = -t62 * t93 + t63 * t77;
t81 = t102 * t63;
t105 = -t64 * t62 + t81;
t35 = -t65 * qJD(3) + (pkin(2) * t65 - qJ(3) * t66) * qJD(2);
t91 = t65 * qJD(2);
t85 = pkin(6) * t91;
t24 = t63 * t35 + t62 * t85;
t14 = (pkin(3) * t65 - pkin(7) * t98) * qJD(2) + t24;
t100 = t62 * t66;
t31 = t62 * t35;
t18 = t31 + (-pkin(6) * t99 - pkin(7) * t100) * qJD(2);
t4 = -qJD(4) * t107 + t102 * t14 - t64 * t18;
t104 = 0.2e1 * t106;
t103 = 2 * qJD(5);
t96 = pkin(7) + qJ(3);
t90 = t66 * qJD(2);
t58 = pkin(6) * t90;
t83 = t62 * t90;
t38 = pkin(3) * t83 + t58;
t44 = pkin(3) * t101 + t65 * pkin(6);
t94 = qJD(3) * t62;
t92 = t63 * qJD(3);
t89 = t66 * qJD(5);
t88 = pkin(6) * t100;
t87 = -0.2e1 * pkin(1) * qJD(2);
t86 = pkin(4) * t91;
t82 = t65 * t90;
t57 = -t63 * pkin(3) - pkin(2);
t80 = qJ(5) * t91;
t79 = t96 * t62;
t76 = t102 * qJD(3);
t75 = 0.2e1 * (t62 ^ 2 + t63 ^ 2) * qJD(3);
t25 = -t63 * t85 + t31;
t72 = -t24 * t62 + t25 * t63;
t71 = t102 * t79;
t48 = t96 * t63;
t12 = qJD(4) * t71 - t63 * t76 + (qJD(4) * t48 + t94) * t64;
t28 = t102 * t48 - t64 * t79;
t70 = -t12 * t66 - t28 * t91;
t13 = t48 * t77 + t64 * t92 + (-t96 * t93 + t76) * t62;
t27 = t64 * t48 + t71;
t69 = t13 * t66 - t27 * t91;
t68 = t102 * t22 - t64 * t26;
t43 = t102 * t62 + t64 * t63;
t3 = -t102 * t18 - t64 * t14 - t22 * t77 + t26 * t93;
t37 = t43 * qJD(4);
t34 = t105 * t65;
t33 = t43 * t65;
t29 = t41 - t88;
t21 = -pkin(4) * t105 - t43 * qJ(5) + t57;
t16 = t106 * t65 + t43 * t90;
t15 = t65 * t37 + t64 * t83 - t81 * t90;
t9 = t37 * pkin(4) - qJ(5) * t106 - t43 * qJD(5);
t8 = t33 * pkin(4) - t34 * qJ(5) + t44;
t7 = t66 * pkin(4) - t68;
t6 = -t66 * qJ(5) + t107;
t5 = t16 * pkin(4) + t15 * qJ(5) - t34 * qJD(5) + t38;
t2 = -t4 - t86;
t1 = -t3 + t80 - t89;
t10 = [0, 0, 0, 0.2e1 * t82, 0.2e1 * (-t65 ^ 2 + t66 ^ 2) * qJD(2), 0, 0, 0, t65 * t87, t66 * t87, -0.2e1 * t24 * t66 + 0.2e1 * (t29 + 0.2e1 * t88) * t91, 0.2e1 * t25 * t66 + 0.2e1 * (-t30 + 0.2e1 * t55) * t91, 0.2e1 * (-t24 * t63 - t25 * t62) * t65 + 0.2e1 * (-t29 * t63 - t30 * t62) * t90, 0.2e1 * pkin(6) ^ 2 * t82 + 0.2e1 * t29 * t24 + 0.2e1 * t30 * t25, -0.2e1 * t34 * t15, 0.2e1 * t15 * t33 - 0.2e1 * t34 * t16, 0.2e1 * t15 * t66 + 0.2e1 * t34 * t91, 0.2e1 * t16 * t66 - 0.2e1 * t33 * t91, -0.2e1 * t82, 0.2e1 * t44 * t16 + 0.2e1 * t38 * t33 - 0.2e1 * t4 * t66 + 0.2e1 * t68 * t91, -0.2e1 * t107 * t91 - 0.2e1 * t44 * t15 - 0.2e1 * t3 * t66 + 0.2e1 * t38 * t34, 0.2e1 * t8 * t16 + 0.2e1 * t2 * t66 + 0.2e1 * t5 * t33 - 0.2e1 * t7 * t91, -0.2e1 * t1 * t33 - 0.2e1 * t7 * t15 - 0.2e1 * t6 * t16 + 0.2e1 * t2 * t34, -0.2e1 * t1 * t66 + 0.2e1 * t8 * t15 - 0.2e1 * t5 * t34 + 0.2e1 * t6 * t91, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, t90, -t91, 0, -t58, t85, t66 * t94 + (t73 * t62 - t55) * qJD(2), t66 * t92 + (t73 * t63 + t88) * qJD(2), t72, -pkin(2) * t58 + (-t29 * t62 + t30 * t63) * qJD(3) + t72 * qJ(3), t106 * t34 - t15 * t43, -t105 * t15 - t106 * t33 - t43 * t16 - t34 * t37, -t106 * t66 + t43 * t91, t105 * t91 + t37 * t66, 0, -t105 * t38 + t57 * t16 + t44 * t37 + t69, t106 * t44 - t57 * t15 + t38 * t43 + t70, -t105 * t5 + t21 * t16 + t9 * t33 + t8 * t37 + t69, t1 * t105 + t106 * t7 + t12 * t33 + t13 * t34 - t27 * t15 - t28 * t16 + t2 * t43 - t6 * t37, -t106 * t8 + t21 * t15 - t9 * t34 - t5 * t43 - t70, t1 * t28 - t6 * t12 + t7 * t13 + t2 * t27 + t5 * t21 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, qJ(3) * t75, t43 * t104, 0.2e1 * t105 * t106 - 0.2e1 * t43 * t37, 0, 0, 0, 0.2e1 * t57 * t37, t57 * t104, -0.2e1 * t105 * t9 + 0.2e1 * t21 * t37, -0.2e1 * t105 * t12 + 0.2e1 * t106 * t27 + 0.2e1 * t13 * t43 - 0.2e1 * t28 * t37, -0.2e1 * t106 * t21 - 0.2e1 * t9 * t43, -0.2e1 * t28 * t12 + 0.2e1 * t27 * t13 + 0.2e1 * t21 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t63 * t90, 0, t58, 0, 0, 0, 0, 0, t16, -t15, t16, 0, t15, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t106, t37, 0, -t106, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, t91, t4, t3, t4 + 0.2e1 * t86, pkin(4) * t15 - t16 * qJ(5) - t33 * qJD(5), -t3 + 0.2e1 * t80 - 0.2e1 * t89, -t2 * pkin(4) + t1 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t37, 0, -t13, t12, -t13, -pkin(4) * t106 - t37 * qJ(5) + qJD(5) * t105, -t12, -t13 * pkin(4) - t12 * qJ(5) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, qJ(5) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
