% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:36
% EndTime: 2020-01-03 11:20:39
% DurationCPUTime: 0.49s
% Computational Cost: add. (409->80), mult. (1061->134), div. (0->0), fcn. (1052->8), ass. (0->85)
t60 = sin(pkin(9));
t64 = sin(qJ(5));
t100 = t64 * t60;
t62 = cos(pkin(9));
t65 = cos(qJ(5));
t97 = t65 * t62;
t67 = -t97 + t100;
t105 = t62 ^ 2;
t106 = t60 ^ 2;
t107 = t105 + t106;
t63 = cos(pkin(8));
t59 = t63 ^ 2;
t61 = sin(pkin(8));
t103 = pkin(6) * t61;
t56 = sin(pkin(7)) * pkin(1) + qJ(3);
t102 = t56 * t63;
t101 = t63 * t61;
t99 = t64 * t62;
t98 = t65 * t60;
t33 = t67 * t61;
t82 = t33 * qJD(1);
t27 = t63 * t82;
t29 = t33 * qJD(5);
t96 = -t27 + t29;
t40 = -cos(pkin(7)) * pkin(1) - pkin(2) - t61 * qJ(4) - t63 * pkin(3);
t21 = t62 * t102 + t60 * t40;
t58 = t61 ^ 2;
t52 = t58 + t59;
t35 = t62 * t40;
t13 = -t62 * t103 + t35 + (-t56 * t60 - pkin(4)) * t63;
t14 = -t60 * t103 + t21;
t3 = -t65 * t13 + t64 * t14;
t68 = -t98 - t99;
t31 = t68 * t61;
t69 = (pkin(4) * t60 + t56) * t61;
t1 = -t3 * t63 + t31 * t69;
t95 = t1 * qJD(1);
t4 = -t64 * t13 - t65 * t14;
t2 = t33 * t69 + t4 * t63;
t94 = t2 * qJD(1);
t20 = -t60 * t102 + t35;
t6 = (t20 * t62 + t21 * t60) * t61;
t93 = t6 * qJD(1);
t7 = t31 ^ 2 - t33 ^ 2;
t92 = t7 * qJD(1);
t10 = t61 * t31 + t68 * t59;
t91 = t10 * qJD(1);
t11 = -t61 * t33 - t67 * t59;
t90 = t11 * qJD(1);
t15 = t68 * t63;
t89 = t15 * qJD(1);
t16 = t67 * t63;
t88 = t16 * qJD(1);
t70 = t106 / 0.2e1 + t105 / 0.2e1;
t23 = (-0.1e1 / 0.2e1 + t70) * t101;
t87 = t23 * qJD(1);
t30 = t52 * t56;
t86 = t30 * qJD(1);
t85 = t31 * qJD(1);
t84 = t31 * qJD(4);
t83 = t31 * qJD(5);
t81 = t33 * qJD(4);
t36 = (0.1e1 / 0.2e1 + t70) * t61;
t80 = t36 * qJD(1);
t41 = t107 * t58;
t79 = t41 * qJD(1);
t42 = t52 * t60;
t78 = t42 * qJD(1);
t43 = t52 * t62;
t77 = t43 * qJD(1);
t76 = t52 * qJD(1);
t75 = qJD(1) * t101;
t74 = qJD(4) * t101;
t73 = t31 * t82;
t72 = t60 * t75;
t71 = t62 * t75;
t5 = t58 * t56 + (-t20 * t60 + t21 * t62) * t63;
t66 = -t5 * qJD(1) - t23 * qJD(2);
t37 = (0.1e1 / 0.2e1 - t107 / 0.2e1) * t61;
t26 = t63 * t85;
t22 = t23 * qJD(3);
t18 = -t15 / 0.2e1 + (-t99 / 0.2e1 - t98 / 0.2e1) * t63;
t17 = (-t67 / 0.2e1 - t97 / 0.2e1 + t100 / 0.2e1) * t63;
t12 = -t26 + t83;
t8 = [0, 0, 0, 0, 0, 0, t52 * qJD(3), t30 * qJD(3), t42 * qJD(3) + t62 * t74, t43 * qJD(3) - t60 * t74, t41 * qJD(4), t5 * qJD(3) - t6 * qJD(4), -t31 * t29, t7 * qJD(5), -t63 * t83, -t63 * t29, 0, -t10 * qJD(3) - t2 * qJD(5) - t63 * t81, t11 * qJD(3) + t1 * qJD(5) + t63 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t76, t86, t78, t77, 0, t37 * qJD(4) - t66, 0, 0, 0, 0, 0, t18 * qJD(5) - t91, t17 * qJD(5) + t90; 0, 0, 0, 0, 0, 0, 0, 0, t71, -t72, t79, t37 * qJD(3) - t93, 0, 0, 0, 0, 0, -t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t92, t12, t96, 0, t18 * qJD(3) + t4 * qJD(5) - t94, t17 * qJD(3) + t3 * qJD(5) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t83; 0, 0, 0, 0, 0, 0, -t76, -t86, -t78, -t77, 0, -t36 * qJD(4) + t66, 0, 0, 0, 0, 0, -t15 * qJD(5) + t91, -t16 * qJD(5) - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * qJD(5) - t89, t67 * qJD(5) - t88; 0, 0, 0, 0, 0, 0, 0, 0, -t71, t72, -t79, t36 * qJD(3) + t93, 0, 0, 0, 0, 0, -t96, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t92, t26, t27, 0, t15 * qJD(3) + t81 + t94, t16 * qJD(3) - t84 - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
