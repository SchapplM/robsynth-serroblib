% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:51
% EndTime: 2019-12-05 16:15:54
% DurationCPUTime: 0.53s
% Computational Cost: add. (1057->116), mult. (1886->163), div. (0->0), fcn. (1240->6), ass. (0->84)
t72 = qJD(2) + qJD(3);
t78 = cos(qJ(3));
t96 = pkin(2) * qJD(3);
t88 = qJD(2) * t96;
t49 = t72 * qJD(4) + t78 * t88;
t73 = sin(pkin(9));
t74 = cos(pkin(9));
t98 = t73 ^ 2 + t74 ^ 2;
t113 = t98 * t49;
t97 = pkin(2) * qJD(2);
t82 = -t78 * t97 + qJD(4);
t77 = cos(qJ(5));
t105 = t77 * t74;
t75 = sin(qJ(5));
t106 = t75 * t73;
t50 = -t105 + t106;
t112 = t50 * t49;
t51 = t77 * t73 + t75 * t74;
t36 = t51 * t72;
t111 = t72 * t98;
t76 = sin(qJ(3));
t110 = t74 * t76;
t109 = t36 ^ 2;
t68 = t74 * pkin(7);
t108 = t78 * pkin(2);
t95 = t72 * t106;
t34 = -t72 * t105 + t95;
t107 = t36 * t34;
t56 = (-pkin(7) - qJ(4)) * t73;
t57 = t74 * qJ(4) + t68;
t24 = t77 * t56 - t75 * t57;
t104 = t24 * qJD(5) - t82 * t50;
t25 = t75 * t56 + t77 * t57;
t103 = -t25 * qJD(5) - t82 * t51;
t65 = -t74 * pkin(4) - pkin(3);
t33 = t65 * t72 + t82;
t45 = t51 * qJD(5);
t62 = t76 * t88;
t102 = t33 * t45 + t50 * t62;
t29 = t72 * t45;
t90 = qJD(5) * t105;
t91 = qJD(5) * t106;
t44 = -t90 + t91;
t101 = -t51 * t29 + t44 * t34;
t100 = -t33 * t44 + t51 * t62;
t93 = t76 * t97;
t53 = t72 * qJ(4) + t93;
t31 = t73 * qJD(1) + t74 * t53;
t94 = t76 * t96;
t67 = t74 * qJD(1);
t22 = t67 + (-pkin(7) * t72 - t53) * t73;
t23 = t72 * t68 + t31;
t12 = t77 * t22 - t75 * t23;
t13 = t75 * t22 + t77 * t23;
t4 = t12 * qJD(5) - t112;
t79 = t51 * t49;
t5 = -t13 * qJD(5) - t79;
t89 = t12 * t44 - t13 * t45 - t4 * t50 - t5 * t51;
t86 = t72 * t94;
t85 = t72 * t93;
t84 = (-qJD(3) + t72) * t97;
t83 = (-qJD(2) - t72) * t96;
t54 = t72 * t90;
t28 = t72 * t91 - t54;
t81 = -t50 * t28 + t36 * t45;
t80 = (-t73 * t53 + t67) * t73 - t31 * t74;
t64 = t76 * pkin(2) + qJ(4);
t47 = (-pkin(7) - t64) * t73;
t48 = t74 * t64 + t68;
t18 = t77 * t47 - t75 * t48;
t19 = t75 * t47 + t77 * t48;
t60 = t78 * t96 + qJD(4);
t58 = t73 * t62;
t55 = t65 - t108;
t52 = -t72 * pkin(3) + t82;
t41 = t45 * qJD(5);
t40 = t44 * qJD(5);
t32 = t34 ^ 2;
t11 = -t19 * qJD(5) - t51 * t60;
t10 = t18 * qJD(5) - t50 * t60;
t9 = t29 * t50 + t34 * t45;
t8 = -t28 * t51 - t36 * t44;
t1 = -t81 + t101;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t40, t81 + t101, -t12 * t45 - t13 * t44 + t4 * t51 - t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 - t86, t78 * t83, 0, 0, 0, 0, 0, 0, 0, 0, t83 * t110, t73 * t86 + t58, t60 * t111 + t113, -t80 * t60 + t64 * t113 + (t52 + (-pkin(3) - t108) * qJD(2)) * t94, t8, t1, -t40, t9, -t41, 0, t11 * qJD(5) + t55 * t29 + t34 * t94 + t102, -t10 * qJD(5) - t55 * t28 + t36 * t94 + t100, -t10 * t34 - t11 * t36 + t18 * t28 - t19 * t29 + t89, t13 * t10 + t12 * t11 + t5 * t18 + t4 * t19 + (qJD(2) * t55 + t33) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 + t85, t78 * t84, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t110, -t73 * t85 + t58, t82 * t111 + t113, -t80 * qJD(4) + qJ(4) * t113 + (t80 * t78 + (-pkin(3) * qJD(3) - t52) * t76) * t97, t8, t1, -t40, t9, -t41, 0, qJD(5) * t103 + t65 * t29 - t34 * t93 + t102, -qJD(5) * t104 - t65 * t28 - t36 * t93 + t100, -t103 * t36 - t104 * t34 + t24 * t28 - t25 * t29 + t89, t5 * t24 + t4 * t25 + t104 * t13 + t103 * t12 + (qJD(3) * t65 - t33) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98 * t72 ^ 2, t72 * t80 + t62, 0, 0, 0, 0, 0, 0, 0.2e1 * t36 * qJD(5), t54 + (-t34 - t95) * qJD(5), -t32 - t109, t12 * t36 + t13 * t34 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t32 + t109, t54 + (t34 - t95) * qJD(5), -t107, 0, 0, -t33 * t36 - t79, t33 * t34 + t112, 0, 0;];
tauc_reg = t2;
