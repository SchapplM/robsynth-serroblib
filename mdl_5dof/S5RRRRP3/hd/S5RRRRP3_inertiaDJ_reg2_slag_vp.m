% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:29
% DurationCPUTime: 0.98s
% Computational Cost: add. (640->110), mult. (1751->155), div. (0->0), fcn. (1099->6), ass. (0->79)
t65 = sin(qJ(4));
t63 = t65 ^ 2;
t68 = cos(qJ(4));
t64 = t68 ^ 2;
t108 = t63 + t64;
t105 = qJD(2) + qJD(3);
t67 = sin(qJ(2));
t101 = pkin(1) * t67;
t66 = sin(qJ(3));
t56 = t66 * t101;
t69 = cos(qJ(3));
t97 = cos(qJ(2));
t81 = t97 * pkin(1);
t75 = t81 + pkin(2);
t71 = t69 * t75;
t76 = qJD(2) * t81;
t16 = -qJD(3) * t71 + t105 * t56 - t69 * t76;
t3 = t108 * t16;
t89 = qJD(3) * pkin(2);
t84 = t69 * t89;
t31 = t108 * t84;
t96 = t67 * t69;
t107 = t105 * pkin(1) * (t97 * t66 + t96);
t104 = 0.2e1 * (-t63 + t64) * qJD(4);
t103 = 2 * qJD(5);
t35 = pkin(1) * t96 + t66 * t75;
t30 = pkin(8) + t35;
t102 = t3 * t30;
t100 = t69 * pkin(2);
t61 = t65 * qJD(4);
t62 = t68 * qJD(4);
t33 = pkin(4) * t61 - qJ(5) * t62 - t65 * qJD(5);
t85 = t66 * t89;
t22 = t33 + t85;
t17 = t85 + t107;
t8 = t17 + t33;
t99 = -t22 - t8;
t98 = -t33 - t8;
t95 = t3 * pkin(8);
t34 = -t56 + t71;
t29 = -pkin(3) - t34;
t94 = t17 * t65 + t29 * t62;
t93 = -t22 - t33;
t58 = t66 * pkin(2) + pkin(8);
t92 = t31 * t58;
t91 = t31 * pkin(8);
t59 = -pkin(3) - t100;
t90 = t59 * t62 + t65 * t85;
t88 = pkin(3) * t61;
t87 = pkin(3) * t62;
t86 = qJD(2) * t101;
t83 = pkin(8) * t61;
t82 = pkin(8) * t62;
t80 = t65 * t62;
t23 = t29 * t61;
t79 = -t17 * t68 + t23;
t78 = -0.2e1 * t85;
t77 = -t3 * t58 + t31 * t30;
t74 = -t68 * pkin(4) - t65 * qJ(5);
t73 = pkin(4) * t65 - qJ(5) * t68;
t42 = t59 * t61;
t72 = -t68 * t85 + t42;
t45 = -pkin(3) + t74;
t32 = t74 * qJD(4) + t68 * qJD(5);
t55 = -0.2e1 * t80;
t54 = 0.2e1 * t80;
t37 = t45 - t100;
t36 = t45 * t61;
t28 = t37 * t61;
t27 = t58 * t62 + t65 * t84;
t26 = t58 * t61 - t68 * t84;
t25 = 0.2e1 * t31;
t19 = t29 + t74;
t18 = t19 * t61;
t5 = -t65 * t16 + t30 * t62;
t4 = t68 * t16 + t30 * t61;
t2 = -0.2e1 * t3;
t1 = t31 - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86, -0.2e1 * t76, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17, 0.2e1 * t16, 0, -0.2e1 * t35 * t16 - 0.2e1 * t34 * t17, t54, t104, 0, t55, 0, 0, 0.2e1 * t79, 0.2e1 * t94, t2, 0.2e1 * t29 * t17 - 0.2e1 * t102, t54, 0, -t104, 0, 0, t55, -0.2e1 * t8 * t68 + 0.2e1 * t18, t2, -0.2e1 * t19 * t62 - 0.2e1 * t8 * t65, 0.2e1 * t19 * t8 - 0.2e1 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t76, 0, 0, 0, 0, 0, 0, 0, 0, t78 - t107, t16 - t84, 0, (-t16 * t66 - t17 * t69 + (-t34 * t66 + t35 * t69) * qJD(3)) * pkin(2), t54, t104, 0, t55, 0, 0, t23 + t42 + (-t17 - t85) * t68, t90 + t94, t1, t17 * t59 + t29 * t85 + t77, t54, 0, -t104, 0, 0, t55, t99 * t68 + t18 + t28, t1, t99 * t65 + (-t19 - t37) * t62, t19 * t22 + t8 * t37 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -0.2e1 * t84, 0, 0, t54, t104, 0, t55, 0, 0, 0.2e1 * t72, 0.2e1 * t90, t25, 0.2e1 * t59 * t85 + 0.2e1 * t92, t54, 0, -t104, 0, 0, t55, -0.2e1 * t22 * t68 + 0.2e1 * t28, t25, -0.2e1 * t22 * t65 - 0.2e1 * t37 * t62, 0.2e1 * t37 * t22 + 0.2e1 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, t54, t104, 0, t55, 0, 0, t79 - t88, -t87 + t94, -t3, -t17 * pkin(3) - t95, t54, 0, -t104, 0, 0, t55, t98 * t68 + t18 + t36, -t3, t98 * t65 + (-t19 - t45) * t62, t19 * t33 + t8 * t45 - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t84, 0, 0, t54, t104, 0, t55, 0, 0, t72 - t88, -t87 + t90, t31, -pkin(3) * t85 + t91, t54, 0, -t104, 0, 0, t55, t93 * t68 + t28 + t36, t31, t93 * t65 + (-t37 - t45) * t62, t22 * t45 + t37 * t33 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t104, 0, t55, 0, 0, -0.2e1 * t88, -0.2e1 * t87, 0, 0, t54, 0, -t104, 0, 0, t55, -0.2e1 * t33 * t68 + 0.2e1 * t36, 0, -0.2e1 * t33 * t65 - 0.2e1 * t45 * t62, 0.2e1 * t45 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t61, 0, -t5, t4, 0, 0, 0, t62, 0, 0, t61, 0, -t5, t32, -t4, t16 * t73 + t30 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t61, 0, -t27, t26, 0, 0, 0, t62, 0, 0, t61, 0, -t27, t32, -t26, t32 * t58 - t73 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t61, 0, -t82, t83, 0, 0, 0, t62, 0, 0, t61, 0, -t82, t32, -t83, t32 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, qJ(5) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
