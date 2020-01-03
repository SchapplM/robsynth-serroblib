% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:54
% EndTime: 2019-12-31 16:57:55
% DurationCPUTime: 0.69s
% Computational Cost: add. (696->103), mult. (1430->136), div. (0->0), fcn. (1500->4), ass. (0->80)
t65 = sin(pkin(6));
t66 = sin(qJ(2));
t67 = cos(qJ(2));
t95 = cos(pkin(6));
t50 = t65 * t67 + t95 * t66;
t111 = t50 ^ 2;
t48 = t65 * t66 - t95 * t67;
t45 = t48 ^ 2;
t14 = t45 - t111;
t121 = t14 * qJD(1);
t120 = t14 * qJD(2);
t116 = t45 + t111;
t119 = qJD(3) * t116;
t118 = t116 * qJD(1);
t117 = qJ(3) + pkin(5);
t54 = t117 * t66;
t55 = t117 * t67;
t27 = t95 * t54 + t65 * t55;
t107 = t65 * t54;
t53 = t95 * t55;
t74 = t53 - t107;
t71 = t27 * t50 - t48 * t74;
t114 = qJD(3) * t71;
t112 = t71 * qJD(1);
t110 = -t48 / 0.2e1;
t73 = t53 / 0.2e1;
t109 = t50 * pkin(3);
t108 = t66 * pkin(2);
t106 = qJD(2) * pkin(2);
t61 = -t67 * pkin(2) - pkin(1);
t18 = t48 * pkin(3) - t50 * qJ(4) + t61;
t101 = t48 * qJ(4);
t19 = t101 + t108 + t109;
t1 = t18 * t19;
t105 = t1 * qJD(1);
t4 = t61 * t108;
t102 = t4 * qJD(1);
t5 = t18 * t50 + t19 * t48;
t100 = t5 * qJD(1);
t6 = t18 * t48 - t19 * t50;
t99 = t6 * qJD(1);
t58 = t65 * pkin(2) + qJ(4);
t60 = -t95 * pkin(2) - pkin(3);
t62 = t108 / 0.2e1;
t9 = t62 + (pkin(3) / 0.2e1 - t60 / 0.2e1) * t50 + (qJ(4) / 0.2e1 + t58 / 0.2e1) * t48;
t96 = t9 * qJD(1);
t94 = qJD(1) * t67;
t69 = t65 * t110 - t95 * t50 / 0.2e1;
t12 = (-t66 / 0.2e1 + t69) * pkin(2);
t93 = t12 * qJD(1);
t16 = t48 * t108 + t61 * t50;
t90 = t16 * qJD(1);
t17 = t50 * t108 - t61 * t48;
t89 = t17 * qJD(1);
t86 = t27 * qJD(2);
t85 = t111 * qJD(1);
t84 = t48 * qJD(1);
t37 = t48 * qJD(2);
t83 = t48 * qJD(3);
t82 = t50 * qJD(1);
t39 = t50 * qJD(2);
t81 = t50 * qJD(4);
t56 = -t66 ^ 2 + t67 ^ 2;
t80 = t56 * qJD(1);
t79 = t66 * qJD(2);
t78 = t67 * qJD(2);
t77 = pkin(1) * t66 * qJD(1);
t76 = pkin(1) * t94;
t24 = t48 * t82;
t23 = t48 * t39;
t75 = t66 * t78;
t25 = t73 - t53 / 0.2e1;
t70 = t25 * qJD(1) + t58 * qJD(2);
t57 = t66 * t94;
t41 = t50 * qJD(3);
t20 = t74 * qJD(2);
t15 = 0.2e1 * t73 - t107;
t11 = t69 * pkin(2) + t62;
t10 = t58 * t110 + t60 * t50 / 0.2e1 + t62 + t101 / 0.2e1 + t109 / 0.2e1;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t56 * qJD(2), 0, -t75, 0, 0, -pkin(1) * t79, -pkin(1) * t78, 0, 0, -t23, t120, 0, t23, 0, 0, t16 * qJD(2), t17 * qJD(2), t119, t4 * qJD(2) + t114, -t23, 0, -t120, 0, 0, t23, t5 * qJD(2) - t48 * t81, t119, t6 * qJD(2) + qJD(4) * t111, t1 * qJD(2) - t18 * t81 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t80, t78, -t57, -t79, 0, -pkin(5) * t78 - t77, pkin(5) * t79 - t76, 0, 0, -t24, t121, -t37, t24, -t39, 0, -t20 + t90, t86 + t89, (t95 * t48 - t50 * t65) * t106, t102 + (-t27 * t65 - t74 * t95) * t106 + t11 * qJD(3), -t24, -t37, -t121, 0, t39, t24, -t20 + t100, (-t60 * t48 - t58 * t50) * qJD(2) - qJD(4) * t48, -t86 + t99, t105 + (-t27 * t58 + t60 * t74) * qJD(2) + t10 * qJD(3) + t15 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t11 * qJD(2) + t112, 0, 0, 0, 0, 0, 0, 0, t118, 0, t10 * qJD(2) + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t37, t85, t15 * qJD(2) - t18 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t80, 0, t57, 0, 0, t77, t76, 0, 0, t24, -t121, 0, -t24, 0, 0, -t41 - t90, t83 - t89, 0, t12 * qJD(3) - t102, t24, 0, t121, 0, 0, -t24, -t41 - t100, 0, -t83 - t99, -t9 * qJD(3) + t25 * qJD(4) - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t58 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t84, 0, t93, 0, 0, 0, 0, 0, 0, -t82, 0, -t84, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, -t118, -t12 * qJD(2) - t112, 0, 0, 0, 0, 0, 0, t39, -t118, t37, t9 * qJD(2) - t112 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t84, 0, -t93, 0, 0, 0, 0, 0, 0, t82, 0, t84, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, -t85, -t25 * qJD(2) + (qJD(1) * t18 + qJD(3)) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
