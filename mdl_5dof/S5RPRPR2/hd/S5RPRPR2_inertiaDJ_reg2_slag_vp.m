% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:09
% DurationCPUTime: 0.60s
% Computational Cost: add. (638->72), mult. (1404->128), div. (0->0), fcn. (1182->8), ass. (0->63)
t57 = sin(pkin(9));
t59 = cos(pkin(9));
t78 = t57 ^ 2 + t59 ^ 2;
t100 = t78 * qJD(4);
t101 = 0.2e1 * t100;
t50 = cos(pkin(8)) * pkin(1) + pkin(2);
t61 = sin(qJ(3));
t62 = cos(qJ(3));
t97 = sin(pkin(8)) * pkin(1);
t98 = -t62 * t50 + t61 * t97;
t34 = t98 * qJD(3);
t33 = qJD(4) - t34;
t99 = t78 * t33;
t80 = t61 * t50 + t62 * t97;
t96 = t59 * pkin(4);
t36 = qJ(4) + t80;
t95 = -pkin(7) - t36;
t94 = cos(qJ(5));
t35 = t80 * qJD(3);
t93 = t35 * t57;
t92 = t35 * t59;
t74 = t94 * t57;
t60 = sin(qJ(5));
t87 = t60 * t59;
t41 = t74 + t87;
t38 = t41 * qJD(5);
t73 = t94 * t59;
t40 = t60 * t57 - t73;
t91 = t40 * t38;
t69 = qJD(5) * t94;
t77 = qJD(5) * t60;
t37 = t57 * t77 - t59 * t69;
t90 = t41 * t37;
t51 = -pkin(3) - t96;
t89 = t51 * t37;
t88 = t51 * t38;
t85 = -pkin(7) - qJ(4);
t65 = -pkin(3) + t98;
t32 = t65 - t96;
t84 = t32 * t38 + t35 * t40;
t83 = -t32 * t37 + t35 * t41;
t54 = t59 * pkin(7);
t29 = t59 * t36 + t54;
t75 = t60 * t95;
t10 = t94 * t29 + t57 * t75;
t64 = t95 * t74;
t4 = -qJD(5) * t64 - t33 * t73 + (qJD(5) * t29 + t33 * t57) * t60;
t5 = -t29 * t69 - t33 * t87 + (-qJD(5) * t75 - t94 * t33) * t57;
t9 = -t60 * t29 + t64;
t76 = -t10 * t38 + t9 * t37 + t4 * t40 - t5 * t41;
t72 = t85 * t57;
t68 = t94 * qJD(4);
t43 = t59 * qJ(4) + t54;
t63 = t94 * t72;
t13 = -qJD(5) * t63 - t59 * t68 + (qJD(4) * t57 + qJD(5) * t43) * t60;
t14 = -t43 * t69 - qJD(4) * t87 + (-t85 * t77 - t68) * t57;
t27 = -t60 * t43 + t63;
t28 = t94 * t43 + t60 * t72;
t66 = t13 * t40 - t14 * t41 + t27 * t37 - t28 * t38;
t23 = -0.2e1 * t90;
t22 = 0.2e1 * t91;
t8 = 0.2e1 * t40 * t37 - 0.2e1 * t41 * t38;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t35, 0.2e1 * t34, 0, -0.2e1 * t80 * t34 + 0.2e1 * t35 * t98, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, 0.2e1 * t93, 0.2e1 * t99, 0.2e1 * t65 * t35 + 0.2e1 * t36 * t99, t23, t8, 0, t22, 0, 0, 0.2e1 * t84, 0.2e1 * t83, 0.2e1 * t76, -0.2e1 * t10 * t4 + 0.2e1 * t32 * t35 + 0.2e1 * t9 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t37 - t9 * t38 - t4 * t41 - t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t90 + 0.2e1 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t93, t100 + t99, -t35 * pkin(3) + qJ(4) * t99 + t100 * t36, t23, t8, 0, t22, 0, 0, t84 + t88, t83 - t89, t66 + t76, -t10 * t13 + t9 * t14 + t5 * t27 - t4 * t28 + t35 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 * t13 - t40 * t14 - t38 * t27 - t37 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, qJ(4) * t101, t23, t8, 0, t22, 0, 0, 0.2e1 * t88, -0.2e1 * t89, 0.2e1 * t66, -0.2e1 * t28 * t13 + 0.2e1 * t27 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, t38, -t37, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t38, 0, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t38, 0, t14, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
