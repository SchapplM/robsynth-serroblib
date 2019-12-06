% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:41
% EndTime: 2019-12-05 17:49:43
% DurationCPUTime: 0.52s
% Computational Cost: add. (629->97), mult. (1294->120), div. (0->0), fcn. (1188->8), ass. (0->84)
t63 = sin(pkin(9));
t61 = t63 ^ 2;
t65 = cos(pkin(9));
t62 = t65 ^ 2;
t55 = t61 + t62;
t49 = t55 * qJ(4);
t67 = sin(qJ(3));
t76 = cos(pkin(8)) * pkin(1) + pkin(2);
t73 = t67 * t76;
t103 = cos(qJ(3));
t106 = pkin(1) * sin(pkin(8));
t77 = t103 * t106;
t39 = t77 + t73;
t37 = qJ(4) + t39;
t11 = t55 * t37;
t66 = sin(qJ(5));
t100 = t66 * t65;
t102 = cos(qJ(5));
t80 = t102 * t63;
t47 = t80 + t100;
t105 = t65 * pkin(4);
t38 = -t103 * t76 + t67 * t106;
t75 = -pkin(3) + t38;
t29 = t75 - t105;
t59 = -pkin(3) - t105;
t82 = t59 / 0.2e1 + t29 / 0.2e1;
t114 = t82 * t47;
t101 = t66 * t63;
t79 = t102 * t65;
t45 = -t79 + t101;
t15 = t45 ^ 2 - t47 ^ 2;
t88 = qJD(1) + qJD(3);
t113 = t88 * t15;
t112 = t88 * t45;
t78 = t88 * t47;
t111 = t88 * t55;
t14 = t55 * t38;
t52 = t55 * qJD(4);
t99 = -t14 * qJD(3) + t52;
t1 = -t37 * t14 + t75 * t39;
t98 = t1 * qJD(1);
t97 = t11 * qJD(1);
t96 = t14 * qJD(1);
t95 = t38 * qJD(1);
t94 = t39 * qJD(1);
t36 = t39 * qJD(3);
t93 = t45 * qJD(1);
t92 = t45 * qJD(3);
t41 = t45 * qJD(5);
t91 = t47 * qJD(1);
t90 = t47 * qJD(3);
t89 = t47 * qJD(5);
t87 = t29 * t93;
t86 = t29 * t91;
t85 = t39 * t93;
t84 = t39 * t91;
t83 = t63 * t94;
t68 = t77 / 0.2e1 + t73 / 0.2e1;
t8 = t68 + (qJ(4) + t37) * (-t62 / 0.2e1 - t61 / 0.2e1);
t74 = t8 * qJD(1) - t49 * qJD(3);
t70 = (t100 / 0.2e1 + t80 / 0.2e1) * t38;
t2 = t70 - t114;
t72 = -t2 * qJD(1) + t59 * t90;
t69 = (t79 / 0.2e1 - t101 / 0.2e1) * t38;
t3 = t82 * t45 + t69;
t71 = -t3 * qJD(1) - t59 * t92;
t60 = t65 * pkin(7);
t51 = t65 * qJ(4) + t60;
t50 = (-pkin(7) - qJ(4)) * t63;
t44 = t47 * qJD(4);
t40 = t45 * qJD(4);
t35 = t38 * qJD(3);
t30 = t63 * t36;
t22 = t65 * t37 + t60;
t21 = (-pkin(7) - t37) * t63;
t20 = t45 * t89;
t19 = t39 * t90;
t18 = t39 * t92;
t13 = t15 * qJD(5);
t10 = t45 * t78;
t9 = t68 + t11 / 0.2e1 + t49 / 0.2e1;
t5 = t70 + t114;
t4 = t69 - (t29 + t59) * t45 / 0.2e1;
t6 = [0, 0, 0, 0, 0, -t36, t35, -t65 * t36, t30, t99, t1 * qJD(3) + t11 * qJD(4), -t20, t13, 0, 0, 0, t29 * t89 + t18, -t29 * t41 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t36 - t94, t35 + t95, -t88 * t65 * t39, t30 + t83, -t96 + t99, t98 + (-t39 * pkin(3) - qJ(4) * t14) * qJD(3) + t9 * qJD(4), -t20, t13, 0, 0, 0, t5 * qJD(5) + t18 + t85, t4 * qJD(5) + t19 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t9 * qJD(3) + t97, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t113, -t41, -t89, 0, t86 + t5 * qJD(3) + (-t102 * t22 - t66 * t21) * qJD(5), -t87 + t4 * qJD(3) + (-t102 * t21 + t66 * t22) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t41; 0, 0, 0, 0, 0, t94, -t95, t65 * t94, -t83, t52 + t96, -t8 * qJD(4) - t98, -t20, t13, 0, 0, 0, -t2 * qJD(5) - t85, -t3 * qJD(5) - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t49 * qJD(4), -t20, t13, 0, 0, 0, t59 * t89, -t59 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t74, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t113, -t41, -t89, 0, (-t102 * t51 - t66 * t50) * qJD(5) + t72, (-t102 * t50 + t66 * t51) * qJD(5) + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t8 * qJD(3) - t97, 0, 0, 0, 0, 0, t89, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t74, 0, 0, 0, 0, 0, t89, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t113, 0, 0, 0, t2 * qJD(3) - t44 - t86, t3 * qJD(3) + t40 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t113, 0, 0, 0, -t44 - t72, t40 - t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
