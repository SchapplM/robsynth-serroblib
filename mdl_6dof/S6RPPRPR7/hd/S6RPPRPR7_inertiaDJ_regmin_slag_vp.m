% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:54
% EndTime: 2019-03-09 01:53:56
% DurationCPUTime: 0.90s
% Computational Cost: add. (1240->121), mult. (2734->235), div. (0->0), fcn. (2695->8), ass. (0->80)
t63 = sin(pkin(10));
t65 = cos(pkin(10));
t68 = sin(qJ(6));
t69 = cos(qJ(6));
t41 = t68 * t63 - t69 * t65;
t36 = t41 * qJD(6);
t103 = sin(qJ(4));
t104 = cos(qJ(4));
t64 = sin(pkin(9));
t66 = cos(pkin(9));
t43 = t103 * t66 + t104 * t64;
t84 = qJD(4) * t103;
t85 = qJD(4) * t104;
t39 = -t64 * t84 + t66 * t85;
t44 = t69 * t63 + t68 * t65;
t97 = t44 * t39;
t111 = t36 * t43 - t97;
t38 = -t64 * t85 - t66 * t84;
t42 = t103 * t64 - t104 * t66;
t12 = t42 * t36 + t38 * t44;
t47 = (t64 ^ 2 + t66 ^ 2) * qJD(3);
t108 = -0.2e1 * t36;
t107 = 2 * qJD(2);
t106 = t38 * pkin(4);
t67 = -pkin(1) - qJ(3);
t105 = -pkin(7) + t67;
t101 = t41 * t39;
t100 = t42 * t38;
t99 = t42 * t63;
t98 = t43 * t39;
t96 = t63 * t38;
t95 = t65 * t38;
t92 = pkin(8) + qJ(5);
t19 = t39 * pkin(4) - t38 * qJ(5) + t42 * qJD(5) + qJD(2);
t45 = t105 * t64;
t46 = t105 * t66;
t21 = t43 * qJD(3) + t45 * t84 - t46 * t85;
t6 = t63 * t19 - t65 * t21;
t54 = t64 * pkin(3) + qJ(2);
t30 = t43 * pkin(4) + t42 * qJ(5) + t54;
t33 = t103 * t46 + t104 * t45;
t14 = t63 * t30 + t65 * t33;
t91 = t63 ^ 2 + t65 ^ 2;
t89 = qJ(2) * qJD(2);
t86 = t91 * t39;
t5 = t65 * t19 + t63 * t21;
t13 = t65 * t30 - t63 * t33;
t83 = t91 * qJD(5);
t82 = 0.2e1 * t83;
t81 = t5 * t65 + t6 * t63;
t80 = -t5 * t63 + t6 * t65;
t7 = t65 * t42 * pkin(8) + t43 * pkin(5) + t13;
t8 = pkin(8) * t99 + t14;
t79 = t68 * t8 - t69 * t7;
t78 = t68 * t7 + t69 * t8;
t32 = t103 * t45 - t104 * t46;
t77 = -t13 * t63 + t14 * t65;
t22 = -qJD(3) * t42 + qJD(4) * t33;
t76 = -t22 * t42 + t32 * t38;
t37 = t44 * qJD(6);
t10 = t42 * t37 - t38 * t41;
t16 = -t37 * t43 - t101;
t48 = t92 * t63;
t49 = t92 * t65;
t74 = -t69 * t48 - t68 * t49;
t73 = -t68 * t48 + t69 * t49;
t71 = -0.2e1 * t98 + 0.2e1 * t100;
t70 = -qJ(5) * t39 - qJD(5) * t43 - t106;
t57 = -t65 * pkin(5) - pkin(4);
t29 = t41 * t42;
t28 = t44 * t42;
t25 = -t44 * qJD(5) - qJD(6) * t73;
t24 = t41 * qJD(5) - t74 * qJD(6);
t23 = -pkin(5) * t99 + t32;
t15 = pkin(5) * t96 + t22;
t4 = -pkin(8) * t96 + t6;
t3 = t39 * pkin(5) - pkin(8) * t95 + t5;
t2 = -t78 * qJD(6) + t69 * t3 - t68 * t4;
t1 = t79 * qJD(6) - t68 * t3 - t69 * t4;
t9 = [0, 0, 0, 0, t107, 0.2e1 * t89, t64 * t107, t66 * t107, 0.2e1 * t47, -0.2e1 * t67 * t47 + 0.2e1 * t89, -0.2e1 * t100, -0.2e1 * t38 * t43 + 0.2e1 * t42 * t39, 0, 0, 0, 0.2e1 * qJD(2) * t43 + 0.2e1 * t54 * t39, -0.2e1 * qJD(2) * t42 + 0.2e1 * t54 * t38, 0.2e1 * t13 * t39 + 0.2e1 * t5 * t43 + 0.2e1 * t76 * t63, -0.2e1 * t14 * t39 - 0.2e1 * t6 * t43 + 0.2e1 * t76 * t65, 0.2e1 * t81 * t42 + 0.2e1 * (-t13 * t65 - t14 * t63) * t38, 0.2e1 * t13 * t5 + 0.2e1 * t14 * t6 + 0.2e1 * t32 * t22, 0.2e1 * t29 * t10, 0.2e1 * t10 * t28 - 0.2e1 * t29 * t12, 0.2e1 * t10 * t43 + 0.2e1 * t29 * t39, -0.2e1 * t12 * t43 + 0.2e1 * t28 * t39, 0.2e1 * t98, 0.2e1 * t23 * t12 - 0.2e1 * t15 * t28 + 0.2e1 * t2 * t43 - 0.2e1 * t79 * t39, 0.2e1 * t1 * t43 + 0.2e1 * t23 * t10 + 0.2e1 * t15 * t29 - 0.2e1 * t78 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, 0, 0, 0, 0, 0, 0, t63 * t71, t65 * t71, 0, t77 * t39 + t80 * t43 - t76, 0, 0, 0, 0, 0, t42 * t12 + t38 * t28 + (t111 - t97) * t43, t42 * t10 - t38 * t29 + (-t16 + t101) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43 * t86 - 0.2e1 * t100, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t39, t38, t65 * t39, -t63 * t39, -t91 * t38, t81, 0, 0, 0, 0, 0, t16, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, -t22, t21, -t22 * t65 + t63 * t70, t22 * t63 + t65 * t70, t80, -t22 * pkin(4) + t80 * qJ(5) + t77 * qJD(5), t10 * t44 - t29 * t36, -t10 * t41 - t44 * t12 - t36 * t28 - t29 * t37, -t111, t16, 0, t57 * t12 + t15 * t41 + t23 * t37 + t25 * t43 + t74 * t39, t57 * t10 + t15 * t44 - t23 * t36 + t24 * t43 - t39 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, t95, -t96, t86, qJ(5) * t86 + t43 * t83 + t106, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, qJ(5) * t82, t44 * t108, 0.2e1 * t36 * t41 - 0.2e1 * t44 * t37, 0, 0, 0, 0.2e1 * t57 * t37, t57 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t95, 0, t22, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t12, t39, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
