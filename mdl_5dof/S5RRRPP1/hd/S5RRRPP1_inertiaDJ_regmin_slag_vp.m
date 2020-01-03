% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:49:59
% DurationCPUTime: 0.43s
% Computational Cost: add. (785->86), mult. (1759->140), div. (0->0), fcn. (1385->6), ass. (0->74)
t73 = cos(qJ(3));
t67 = t73 * qJD(4);
t71 = sin(qJ(3));
t74 = cos(qJ(2));
t96 = pkin(1) * qJD(2);
t90 = t74 * t96;
t79 = t73 * t90;
t72 = sin(qJ(2));
t62 = pkin(1) * t72 + pkin(7);
t95 = -qJ(4) - t62;
t80 = qJD(3) * t95;
t38 = t71 * t80 + t67 + t79;
t70 = sin(pkin(8));
t75 = (-qJD(4) - t90) * t71 + t73 * t80;
t94 = cos(pkin(8));
t12 = t70 * t38 - t75 * t94;
t13 = t38 * t94 + t70 * t75;
t69 = t73 * qJ(4);
t46 = t62 * t73 + t69;
t85 = t94 * t71;
t34 = t70 * t46 - t85 * t95;
t101 = t70 * t71;
t35 = t101 * t95 + t46 * t94;
t48 = t70 * t73 + t85;
t44 = t48 * qJD(3);
t84 = t94 * t73;
t93 = t71 * qJD(3);
t45 = qJD(3) * t84 - t70 * t93;
t47 = -t84 + t101;
t112 = t12 * t48 - t13 * t47 + t34 * t45 - t35 * t44;
t116 = 0.2e1 * t112;
t100 = -qJ(4) - pkin(7);
t83 = qJD(3) * t100;
t43 = t71 * t83 + t67;
t76 = -t71 * qJD(4) + t73 * t83;
t31 = t70 * t43 - t76 * t94;
t32 = t43 * t94 + t70 * t76;
t55 = pkin(7) * t73 + t69;
t39 = -t100 * t85 + t70 * t55;
t40 = t100 * t101 + t55 * t94;
t113 = t31 * t48 - t32 * t47 + t39 * t45 - t40 * t44;
t115 = 0.2e1 * t113;
t114 = t113 + t112;
t111 = 2 * qJD(5);
t110 = t74 * pkin(1);
t65 = pkin(3) * t93;
t16 = pkin(4) * t44 - qJ(5) * t45 - qJD(5) * t48 + t65;
t66 = t72 * t96;
t14 = t16 + t66;
t64 = -pkin(3) * t73 - pkin(2);
t36 = pkin(4) * t47 - qJ(5) * t48 + t64;
t33 = t36 - t110;
t108 = t14 * t47 + t33 * t44;
t107 = -t14 * t48 - t33 * t45;
t106 = t16 * t47 + t36 * t44;
t99 = -t16 * t48 - t36 * t45;
t63 = -pkin(2) - t110;
t68 = t73 * qJD(3);
t97 = t63 * t68 + t71 * t66;
t92 = pkin(2) * t93;
t91 = pkin(2) * t68;
t89 = t12 * t34 + t35 * t13;
t86 = t31 * t39 + t40 * t32;
t78 = t63 * t93 - t66 * t73;
t77 = t12 * t39 + t13 * t40 + t31 * t34 + t35 * t32;
t61 = -pkin(3) * t94 - pkin(4);
t59 = pkin(3) * t70 + qJ(5);
t57 = 0.2e1 * t71 * t68;
t54 = t64 - t110;
t51 = t66 + t65;
t50 = 0.2e1 * (-t71 ^ 2 + t73 ^ 2) * qJD(3);
t30 = (-t44 * t70 - t45 * t94) * pkin(3);
t15 = -qJD(5) * t47 - t44 * t59 + t45 * t61;
t1 = [0, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t90, t57, t50, 0, 0, 0, 0.2e1 * t78, 0.2e1 * t97, t116, 0.2e1 * t51 * t54 + 0.2e1 * t89, 0.2e1 * t108, t116, 0.2e1 * t107, 0.2e1 * t14 * t33 + 0.2e1 * t89; 0, 0, 0, 0, -t66, -t90, t57, t50, 0, 0, 0, t78 - t92, -t91 + t97, t114, t51 * t64 + t54 * t65 + t77, t106 + t108, t114, t99 + t107, t14 * t36 + t16 * t33 + t77; 0, 0, 0, 0, 0, 0, t57, t50, 0, 0, 0, -0.2e1 * t92, -0.2e1 * t91, t115, 0.2e1 * t64 * t65 + 0.2e1 * t86, 0.2e1 * t106, t115, 0.2e1 * t99, 0.2e1 * t16 * t36 + 0.2e1 * t86; 0, 0, 0, 0, 0, 0, 0, 0, t68, -t93, 0, -t62 * t68 - t71 * t90, t62 * t93 - t79, t30, (-t12 * t94 + t13 * t70) * pkin(3), -t12, t15, t13, qJD(5) * t35 + t12 * t61 + t13 * t59; 0, 0, 0, 0, 0, 0, 0, 0, t68, -t93, 0, -pkin(7) * t68, pkin(7) * t93, t30, (-t31 * t94 + t32 * t70) * pkin(3), -t31, t15, t32, qJD(5) * t40 + t31 * t61 + t32 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t59 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t44, 0, -t45, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t44, 0, -t45, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
