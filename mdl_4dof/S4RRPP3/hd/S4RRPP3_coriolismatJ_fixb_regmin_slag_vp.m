% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:07
% EndTime: 2021-01-15 10:36:09
% DurationCPUTime: 0.46s
% Computational Cost: add. (638->88), mult. (1286->133), div. (0->0), fcn. (1310->4), ass. (0->74)
t59 = sin(pkin(6));
t60 = sin(qJ(2));
t61 = cos(qJ(2));
t88 = cos(pkin(6));
t45 = t59 * t61 + t88 * t60;
t42 = t45 ^ 2;
t43 = t59 * t60 - t88 * t61;
t108 = t43 ^ 2 + t42;
t111 = qJD(3) * t108;
t110 = t108 * qJD(1);
t109 = -qJ(3) - pkin(5);
t49 = t109 * t60;
t50 = t109 * t61;
t23 = -t88 * t49 - t59 * t50;
t100 = t59 * t49;
t48 = t88 * t50;
t68 = -t48 + t100;
t65 = t23 * t45 - t43 * t68;
t106 = qJD(3) * t65;
t104 = t65 * qJD(1);
t103 = -t43 / 0.2e1;
t67 = -t48 / 0.2e1;
t102 = t45 * pkin(3);
t101 = t60 * pkin(2);
t99 = qJD(2) * pkin(2);
t55 = -t61 * pkin(2) - pkin(1);
t16 = t43 * pkin(3) - t45 * qJ(4) + t55;
t94 = t43 * qJ(4);
t17 = t101 + t94 + t102;
t1 = t16 * t17;
t98 = t1 * qJD(1);
t4 = t55 * t101;
t95 = t4 * qJD(1);
t5 = t16 * t45 + t17 * t43;
t93 = t5 * qJD(1);
t6 = t16 * t43 - t17 * t45;
t92 = t6 * qJD(1);
t52 = t59 * pkin(2) + qJ(4);
t54 = -t88 * pkin(2) - pkin(3);
t56 = t101 / 0.2e1;
t9 = t56 + (pkin(3) / 0.2e1 - t54 / 0.2e1) * t45 + (qJ(4) / 0.2e1 + t52 / 0.2e1) * t43;
t89 = t9 * qJD(1);
t87 = qJD(1) * t61;
t63 = t59 * t103 - t88 * t45 / 0.2e1;
t12 = (-t60 / 0.2e1 + t63) * pkin(2);
t86 = t12 * qJD(1);
t14 = t43 * t101 + t55 * t45;
t85 = t14 * qJD(1);
t15 = t45 * t101 - t55 * t43;
t84 = t15 * qJD(1);
t81 = t23 * qJD(2);
t80 = t42 * qJD(1);
t79 = t43 * qJD(1);
t33 = t43 * qJD(2);
t78 = t43 * qJD(3);
t77 = t45 * qJD(1);
t76 = t45 * qJD(4);
t51 = -t60 ^ 2 + t61 ^ 2;
t75 = t51 * qJD(1);
t74 = t60 * qJD(2);
t73 = t61 * qJD(2);
t72 = qJD(1) * pkin(1) * t60;
t71 = pkin(1) * t87;
t70 = t43 * t77;
t69 = t60 * t87;
t21 = t67 + t48 / 0.2e1;
t64 = t21 * qJD(1) + t52 * qJD(2);
t37 = t45 * qJD(3);
t35 = t45 * qJD(2);
t18 = t68 * qJD(2);
t13 = 0.2e1 * t67 + t100;
t11 = t63 * pkin(2) + t56;
t10 = t52 * t103 + t54 * t45 / 0.2e1 + t56 + t94 / 0.2e1 + t102 / 0.2e1;
t2 = [0, 0, 0, t60 * t73, t51 * qJD(2), 0, 0, 0, -pkin(1) * t74, -pkin(1) * t73, t14 * qJD(2), t15 * qJD(2), t111, t4 * qJD(2) + t106, t5 * qJD(2) - t43 * t76, t111, t6 * qJD(2) + t42 * qJD(4), t1 * qJD(2) - t16 * t76 + t106; 0, 0, 0, t69, t75, t73, -t74, 0, -pkin(5) * t73 - t72, pkin(5) * t74 - t71, -t18 + t85, t81 + t84, (t88 * t43 - t45 * t59) * t99, t95 + (-t23 * t59 - t68 * t88) * t99 + t11 * qJD(3), -t18 + t93, (-t54 * t43 - t52 * t45) * qJD(2) - qJD(4) * t43, -t81 + t92, t98 + (-t23 * t52 + t54 * t68) * qJD(2) + t10 * qJD(3) + t13 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t11 * qJD(2) + t104, 0, t110, 0, t10 * qJD(2) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t33, t80, t13 * qJD(2) - t16 * t77; 0, 0, 0, -t69, -t75, 0, 0, 0, t72, t71, -t37 - t85, t78 - t84, 0, t12 * qJD(3) - t95, -t37 - t93, 0, -t78 - t92, -t9 * qJD(3) + t21 * qJD(4) - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t52 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t79, 0, t86, -t77, 0, -t79, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, -t110, -t12 * qJD(2) - t104, t35, -t110, t33, t9 * qJD(2) - t104 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t79, 0, -t86, t77, 0, t79, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t80, -t21 * qJD(2) + (qJD(1) * t16 + qJD(3)) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
