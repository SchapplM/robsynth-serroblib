% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:59
% EndTime: 2019-12-31 18:41:01
% DurationCPUTime: 0.64s
% Computational Cost: add. (811->124), mult. (1496->140), div. (0->0), fcn. (1172->6), ass. (0->105)
t118 = cos(qJ(3));
t121 = pkin(1) * sin(pkin(8));
t69 = sin(qJ(3));
t80 = cos(pkin(8)) * pkin(1) + pkin(2);
t44 = -t118 * t80 + t69 * t121;
t68 = sin(qJ(4));
t70 = cos(qJ(4));
t79 = -t70 * pkin(4) - t68 * qJ(5);
t52 = -pkin(3) + t79;
t25 = t44 + t52;
t22 = t25 * t68;
t51 = t52 * t68;
t108 = t22 / 0.2e1 + t51 / 0.2e1;
t104 = t70 * qJ(5);
t119 = t68 * pkin(4);
t54 = t104 - t119;
t111 = t54 * t70;
t124 = -t108 - t111;
t42 = -pkin(3) + t44;
t90 = pkin(3) / 0.2e1 - t42 / 0.2e1;
t123 = t90 * t68;
t65 = t68 ^ 2;
t66 = t70 ^ 2;
t57 = t66 - t65;
t93 = qJD(1) + qJD(3);
t122 = t93 * t57;
t120 = pkin(3) * t70;
t117 = t25 * t54;
t116 = t25 * t70;
t115 = t42 * t70;
t114 = t52 * t70;
t113 = t54 * t52;
t112 = t54 * t68;
t110 = t68 * t44;
t109 = t70 * t44;
t45 = t118 * t121 + t69 * t80;
t41 = t45 * qJD(3);
t34 = t68 * t41;
t62 = t65 * qJD(5);
t106 = t62 - t34;
t17 = (-t65 - t66) * t44;
t43 = pkin(7) + t45;
t1 = t43 * t17 + t25 * t45;
t105 = t1 * qJD(1);
t103 = qJD(1) * t68;
t102 = qJD(3) * t68;
t15 = -t112 + t116;
t101 = t15 * qJD(1);
t16 = -t22 - t111;
t100 = t16 * qJD(1);
t99 = t17 * qJD(1);
t98 = t44 * qJD(1);
t97 = t45 * qJD(1);
t96 = t68 * qJD(4);
t63 = t70 * qJD(4);
t95 = t70 * qJD(5);
t94 = qJD(4) * qJ(5);
t92 = pkin(7) * t96;
t91 = pkin(7) * t63;
t89 = qJD(1) * t117;
t88 = t25 * t103;
t87 = t42 * t103;
t86 = qJD(1) * t115;
t85 = t70 * t41;
t84 = t43 * t96;
t83 = t43 * t63;
t82 = t52 / 0.2e1 + t25 / 0.2e1;
t58 = t68 * t95;
t81 = t58 - t85;
t23 = -t112 + t114;
t32 = -t109 / 0.2e1;
t5 = t82 * t70 - t112 + t32;
t78 = t5 * qJD(1) + t23 * qJD(3);
t24 = -t51 - t111;
t31 = t110 / 0.2e1;
t4 = t31 + t124;
t77 = t4 * qJD(1) + t24 * qJD(3);
t76 = t54 * qJD(4) + t68 * qJD(5);
t10 = t31 + t123;
t75 = pkin(3) * t102 + t10 * qJD(1);
t33 = t109 / 0.2e1;
t11 = t90 * t70 + t33;
t74 = t11 * qJD(1) + qJD(3) * t120;
t71 = (-t104 / 0.2e1 + t119 / 0.2e1) * t44;
t2 = t82 * t54 + t71;
t73 = t2 * qJD(1) + qJD(3) * t113;
t30 = -t110 / 0.2e1;
t8 = t30 + t108;
t72 = t8 * qJD(1) + t52 * t102;
t46 = t79 * qJD(4) + t95;
t59 = t68 * t63;
t56 = t57 * qJD(4);
t53 = t93 * t65;
t47 = t93 * t70 * t68;
t40 = t44 * qJD(3);
t36 = t70 * t97;
t35 = t68 * t97;
t14 = t17 * qJD(3);
t13 = -t120 / 0.2e1 + t115 / 0.2e1 + t33;
t12 = t31 - t123;
t9 = t30 - t108;
t7 = t31 - t124;
t6 = -t116 / 0.2e1 + t112 + t32 - t114 / 0.2e1;
t3 = -t113 / 0.2e1 - t117 / 0.2e1 + t71;
t18 = [0, 0, 0, 0, 0, -t41, t40, t59, t56, 0, 0, 0, t42 * t96 - t85, t42 * t63 + t34, -t16 * qJD(4) + t81, t14, -t15 * qJD(4) + t106, t1 * qJD(3) - t76 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t41 - t97, t40 + t98, t59, t56, 0, 0, 0, t12 * qJD(4) - t36 - t85, t13 * qJD(4) + t34 + t35, t7 * qJD(4) - t36 + t81, t14 + t99, t6 * qJD(4) + t106 - t35, t105 + (pkin(7) * t17 + t45 * t52) * qJD(3) + t3 * qJD(4) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t47, t122, t63, -t96, 0, t12 * qJD(3) - t83 + t87, t13 * qJD(3) + t84 + t86, t7 * qJD(3) - t100 - t83, t46, t6 * qJD(3) - t101 - t84, t3 * qJD(3) + t46 * t43 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t63, t53, t9 * qJD(3) + t83 - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t63, -t96, 0, t63, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96; 0, 0, 0, 0, 0, t97, -t98, t59, t56, 0, 0, 0, -t10 * qJD(4) + t36, -t11 * qJD(4) - t35, -t4 * qJD(4) + t36 + t58, -t99, -t5 * qJD(4) + t35 + t62, -t2 * qJD(4) - t8 * qJD(5) - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t59, t56, 0, 0, 0, -pkin(3) * t96, -pkin(3) * t63, -t24 * qJD(4) + t58, 0, -t23 * qJD(4) + t62, -t76 * t52; 0, 0, 0, 0, 0, 0, 0, t47, t122, t63, -t96, 0, -t75 - t91, -t74 + t92, -t77 - t91, t46, -t78 - t92, t46 * pkin(7) - t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t63, t53, -t72 + t91; 0, 0, 0, 0, 0, 0, 0, -t47, -t122, 0, 0, 0, t10 * qJD(3) - t87, t11 * qJD(3) - t86, t4 * qJD(3) + t100, 0, t5 * qJD(3) + t101, t2 * qJD(3) + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t47, -t122, 0, 0, 0, t75, t74, t77, 0, t78, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, -t53, t8 * qJD(3) + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, -t53, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t18;
