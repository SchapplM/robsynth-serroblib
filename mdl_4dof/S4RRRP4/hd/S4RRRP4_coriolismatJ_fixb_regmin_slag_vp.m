% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:43
% EndTime: 2019-12-31 17:15:44
% DurationCPUTime: 0.58s
% Computational Cost: add. (860->90), mult. (1709->138), div. (0->0), fcn. (1694->4), ass. (0->83)
t115 = pkin(5) + pkin(6);
t63 = sin(qJ(2));
t55 = t115 * t63;
t62 = sin(qJ(3));
t108 = t62 * t55;
t109 = cos(qJ(3));
t64 = cos(qJ(2));
t56 = t115 * t64;
t54 = t109 * t56;
t119 = -t54 + t108;
t49 = -t109 * t64 + t62 * t63;
t19 = t49 * qJ(4) + t119;
t111 = t19 * pkin(3);
t83 = qJD(2) + qJD(3);
t78 = t109 * pkin(2);
t25 = t109 * t55 + t56 * t62;
t51 = t109 * t63 + t62 * t64;
t121 = -qJ(4) * t51 - t25;
t120 = t83 * t25;
t118 = t51 ^ 2;
t117 = -pkin(3) / 0.2e1;
t69 = -t54 / 0.2e1;
t59 = t78 + pkin(3);
t116 = -t59 / 0.2e1;
t114 = pkin(2) * t62;
t113 = pkin(3) * t49;
t112 = pkin(3) * t51;
t110 = t63 * pkin(2);
t60 = -pkin(2) * t64 - pkin(1);
t37 = t60 + t113;
t5 = t37 * (t110 + t112);
t103 = t5 * qJD(1);
t6 = t37 * t112;
t102 = t6 * qJD(1);
t7 = -t121 * t51 + t19 * t49;
t101 = t7 * qJD(1);
t74 = -t62 * t49 / 0.2e1;
t9 = (t116 + t117) * t51 + (t74 - t63 / 0.2e1) * pkin(2);
t100 = t9 * qJD(1);
t99 = qJD(1) * t51;
t98 = qJD(1) * t60;
t97 = qJD(1) * t64;
t96 = qJD(2) * t63;
t95 = qJD(2) * t64;
t67 = -t78 / 0.2e1 + t59 / 0.2e1;
t11 = (t117 + t67) * t49;
t94 = t11 * qJD(1);
t48 = t49 ^ 2;
t24 = t48 - t118;
t93 = t24 * qJD(1);
t27 = t110 * t49 + t51 * t60;
t90 = t27 * qJD(1);
t28 = t110 * t51 - t49 * t60;
t89 = t28 * qJD(1);
t31 = t48 + t118;
t88 = t31 * qJD(1);
t34 = t69 + t54 / 0.2e1;
t87 = t34 * qJD(1);
t86 = t49 * qJD(3);
t85 = t51 * qJD(3);
t57 = -t63 ^ 2 + t64 ^ 2;
t84 = t57 * qJD(1);
t82 = pkin(1) * t63 * qJD(1);
t81 = pkin(1) * t97;
t80 = pkin(3) * t99;
t79 = qJD(3) * t114;
t77 = t49 * t98;
t76 = t51 * t98;
t75 = t63 * t97;
t72 = t109 * qJD(2);
t71 = t109 * qJD(3);
t30 = t83 * t51;
t65 = (t78 / 0.2e1 + t116) * t19;
t2 = -t111 / 0.2e1 - t65;
t40 = (t78 - t59) * t114;
t68 = -qJD(1) * t2 - qJD(2) * t40;
t32 = t49 * t99;
t29 = t83 * t49;
t26 = 0.2e1 * t69 + t108;
t10 = t113 / 0.2e1 + t67 * t49;
t8 = pkin(2) * t74 + t51 * t116 + t110 / 0.2e1 + t112 / 0.2e1;
t1 = t111 / 0.2e1 - t65;
t3 = [0, 0, 0, t63 * t95, t57 * qJD(2), 0, 0, 0, -pkin(1) * t96, -pkin(1) * t95, -t49 * t30, t83 * t24, 0, 0, 0, qJD(2) * t27 + t60 * t85, qJD(2) * t28 - t60 * t86, t31 * qJD(4), qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t7; 0, 0, 0, t75, t84, t95, -t96, 0, -pkin(5) * t95 - t82, pkin(5) * t96 - t81, -t32, t93, -t29, -t30, 0, qJD(2) * t119 + qJD(3) * t26 + t90, t89 + t120, (-t114 * t51 + t59 * t49) * qJD(2) + t10 * qJD(3), t103 + (t114 * t121 + t19 * t59) * qJD(2) + t1 * qJD(3) + t8 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t93, -t29, -t30, 0, qJD(2) * t26 + qJD(3) * t119 + t76, -t77 + t120, pkin(3) * t86 + qJD(2) * t10, qJD(2) * t1 + qJD(3) * t111 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, qJD(2) * t8 + t101; 0, 0, 0, -t75, -t84, 0, 0, 0, t82, t81, t32, -t93, 0, 0, 0, qJD(3) * t34 - t90, -t89, t11 * qJD(3), qJD(3) * t2 + qJD(4) * t9 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -pkin(2) * t71, 0, t40 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114 * t83 + t87, (-t72 - t71) * pkin(2), t94, -pkin(3) * t79 - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t93, 0, 0, 0, -qJD(2) * t34 - t76, t77, -t11 * qJD(2), -qJD(2) * t2 - qJD(4) * t112 - t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t114 - t87, pkin(2) * t72, -t94, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, pkin(3) * t85 - qJD(2) * t9 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
