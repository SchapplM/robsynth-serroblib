% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:46
% DurationCPUTime: 0.66s
% Computational Cost: add. (645->80), mult. (1484->140), div. (0->0), fcn. (1250->8), ass. (0->70)
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t83 = t58 ^ 2 + t60 ^ 2;
t104 = qJD(4) * t83;
t105 = 0.2e1 * t104;
t59 = sin(pkin(8));
t61 = cos(pkin(8));
t64 = cos(qJ(2));
t82 = pkin(1) * qJD(2);
t80 = t64 * t82;
t63 = sin(qJ(2));
t81 = t63 * t82;
t37 = -t59 * t81 + t61 * t80;
t34 = qJD(4) + t37;
t103 = t34 * t83;
t102 = t60 * pkin(4);
t52 = t64 * pkin(1) + pkin(2);
t92 = t61 * t63;
t85 = pkin(1) * t92 + t59 * t52;
t35 = qJ(4) + t85;
t101 = -pkin(7) - t35;
t50 = t59 * pkin(2) + qJ(4);
t100 = -pkin(7) - t50;
t99 = cos(qJ(5));
t36 = (t59 * t64 + t92) * t82;
t98 = t36 * t58;
t97 = t36 * t60;
t75 = t99 * t58;
t62 = sin(qJ(5));
t90 = t62 * t60;
t42 = t75 + t90;
t39 = t42 * qJD(5);
t74 = t99 * t60;
t91 = t62 * t58;
t41 = -t74 + t91;
t96 = t41 * t39;
t70 = qJD(5) * t99;
t38 = qJD(5) * t91 - t60 * t70;
t95 = t42 * t38;
t78 = -t61 * pkin(2) - pkin(3);
t43 = t78 - t102;
t94 = t43 * t38;
t93 = t43 * t39;
t71 = -t59 * t63 * pkin(1) + t61 * t52;
t67 = -pkin(3) - t71;
t32 = t67 - t102;
t89 = t32 * t39 + t36 * t41;
t88 = -t32 * t38 + t36 * t42;
t55 = t60 * pkin(7);
t28 = t60 * t35 + t55;
t66 = t101 * t75;
t11 = -t62 * t28 + t66;
t77 = t62 * t101;
t12 = t99 * t28 + t58 * t77;
t4 = -qJD(5) * t66 - t34 * t74 + (qJD(5) * t28 + t34 * t58) * t62;
t5 = -t28 * t70 - t34 * t90 + (-qJD(5) * t77 - t99 * t34) * t58;
t79 = t11 * t38 - t12 * t39 + t4 * t41 - t5 * t42;
t76 = t62 * t100;
t40 = t60 * t50 + t55;
t65 = t100 * t75;
t69 = t99 * qJD(4);
t13 = -qJD(5) * t65 - t60 * t69 + (qJD(4) * t58 + qJD(5) * t40) * t62;
t14 = -t40 * t70 - qJD(4) * t90 + (-qJD(5) * t76 - t69) * t58;
t19 = -t62 * t40 + t65;
t20 = t99 * t40 + t58 * t76;
t72 = t13 * t41 - t14 * t42 + t19 * t38 - t20 * t39;
t25 = -0.2e1 * t95;
t24 = 0.2e1 * t96;
t10 = 0.2e1 * t41 * t38 - 0.2e1 * t42 * t39;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t81, -0.2e1 * t80, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t36, -0.2e1 * t37, 0, -0.2e1 * t71 * t36 + 0.2e1 * t85 * t37, 0, 0, 0, 0, 0, 0, -0.2e1 * t97, 0.2e1 * t98, 0.2e1 * t103, 0.2e1 * t103 * t35 + 0.2e1 * t67 * t36, t25, t10, 0, t24, 0, 0, 0.2e1 * t89, 0.2e1 * t88, 0.2e1 * t79, 0.2e1 * t11 * t5 - 0.2e1 * t12 * t4 + 0.2e1 * t32 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t80, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, (-t36 * t61 + t37 * t59) * pkin(2), 0, 0, 0, 0, 0, 0, -t97, t98, t104 + t103, t103 * t50 + t104 * t35 + t36 * t78, t25, t10, 0, t24, 0, 0, t89 + t93, t88 - t94, t72 + t79, t11 * t14 - t12 * t13 + t5 * t19 - t4 * t20 + t36 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t50 * t105, t25, t10, 0, t24, 0, 0, 0.2e1 * t93, -0.2e1 * t94, 0.2e1 * t72, -0.2e1 * t20 * t13 + 0.2e1 * t19 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t39 - t12 * t38 - t4 * t42 - t5 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t42 - t14 * t41 - t19 * t39 - t20 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95 + 0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t39, -t38, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t39, 0, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t39, 0, t14, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
