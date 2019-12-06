% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:20
% EndTime: 2019-12-05 15:01:23
% DurationCPUTime: 0.48s
% Computational Cost: add. (767->108), mult. (2065->159), div. (0->0), fcn. (1669->8), ass. (0->77)
t62 = sin(pkin(8));
t66 = sin(qJ(3));
t91 = t66 * t62;
t52 = qJD(1) * t91;
t64 = cos(pkin(8));
t68 = cos(qJ(3));
t89 = t68 * t64;
t34 = qJD(1) * t89 - t52;
t71 = qJD(4) - t34;
t80 = qJD(1) * qJD(3);
t75 = t68 * t80;
t51 = t64 * t75;
t24 = t51 + (qJD(4) - t52) * qJD(3);
t63 = cos(pkin(9));
t67 = cos(qJ(5));
t61 = sin(pkin(9));
t65 = sin(qJ(5));
t92 = t65 * t61;
t42 = -t67 * t63 + t92;
t97 = t42 * t24;
t45 = t68 * t62 + t66 * t64;
t12 = t42 * t45;
t76 = t66 * t80;
t29 = t62 * t75 + t64 * t76;
t37 = t45 * qJD(1);
t74 = t37 * qJD(3) - t29;
t44 = t67 * t61 + t65 * t63;
t35 = t44 * qJD(3);
t96 = t35 ^ 2;
t88 = pkin(6) + qJ(4);
t46 = t88 * t61;
t47 = t88 * t63;
t21 = -t65 * t46 + t67 * t47;
t95 = -t21 * qJD(5) - t71 * t44;
t20 = -t67 * t46 - t65 * t47;
t94 = -t20 * qJD(5) + t71 * t42;
t43 = -t89 + t91;
t15 = t29 * t43;
t85 = qJD(3) * t63;
t78 = t67 * t85;
t79 = qJD(3) * t92;
t32 = -t78 + t79;
t93 = t35 * t32;
t39 = t44 * qJD(5);
t28 = qJD(3) * t39;
t38 = t42 * qJD(5);
t87 = -t44 * t28 + t38 * t32;
t26 = qJD(3) * qJ(4) + t37;
t19 = t61 * qJD(2) + t63 * t26;
t86 = t61 ^ 2 + t63 ^ 2;
t83 = t38 * qJD(5);
t40 = t43 * qJD(3);
t82 = t40 * qJD(3);
t41 = t45 * qJD(3);
t81 = t41 * qJD(3);
t56 = -t63 * pkin(4) - pkin(3);
t77 = t86 * t24;
t58 = t63 * qJD(2);
t13 = t58 + (-pkin(6) * qJD(3) - t26) * t61;
t14 = pkin(6) * t85 + t19;
t5 = t67 * t13 - t65 * t14;
t6 = t65 * t13 + t67 * t14;
t73 = (-t61 * t26 + t58) * t61 - t19 * t63;
t48 = qJD(5) * t78;
t27 = qJD(5) * t79 - t48;
t72 = -t42 * t27 + t35 * t39;
t70 = t44 * t24;
t11 = t44 * t45;
t31 = t32 ^ 2;
t30 = t39 * qJD(5);
t25 = -qJD(3) * pkin(3) + t71;
t22 = t56 * qJD(3) + t71;
t4 = qJD(5) * t12 + t44 * t40;
t3 = -qJD(5) * t11 + t42 * t40;
t2 = -t6 * qJD(5) - t70;
t1 = t5 * qJD(5) - t97;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t82, 0, (-t62 * t76 + t51) * t45 - t37 * t40 + t15 - t34 * t41, 0, 0, 0, 0, 0, 0, -t63 * t81, t61 * t81, -t86 * t82, t25 * t41 + t73 * t40 + t45 * t77 + t15, 0, 0, 0, 0, 0, 0, t4 * qJD(5) + t43 * t28 + t41 * t32, -t3 * qJD(5) - t43 * t27 + t41 * t35, -t11 * t27 + t12 * t28 - t3 * t32 - t4 * t35, -t1 * t12 - t2 * t11 + t22 * t41 + t6 * t3 + t5 * t4 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t83, t72 + t87, t1 * t44 - t2 * t42 - t6 * t38 - t5 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t51 + (t34 + t52) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, t74 * t63, -t74 * t61, t71 * qJD(3) * t86 + t77, -t29 * pkin(3) + qJ(4) * t77 - t25 * t37 - t71 * t73, -t27 * t44 - t35 * t38, -t72 + t87, -t83, t28 * t42 + t32 * t39, -t30, 0, t95 * qJD(5) + t22 * t39 + t56 * t28 + t29 * t42 - t37 * t32, t94 * qJD(5) - t22 * t38 - t56 * t27 + t29 * t44 - t37 * t35, -t1 * t42 - t2 * t44 + t20 * t27 - t21 * t28 + t94 * t32 - t95 * t35 + t5 * t38 - t6 * t39, t1 * t21 + t2 * t20 - t22 * t37 + t29 * t56 + t95 * t5 - t94 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 * qJD(3) ^ 2, t73 * qJD(3) + t29, 0, 0, 0, 0, 0, 0, 0.2e1 * t35 * qJD(5), t48 + (-t32 - t79) * qJD(5), -t31 - t96, t6 * t32 + t5 * t35 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t31 + t96, t48 + (t32 - t79) * qJD(5), -t93, 0, 0, -t22 * t35 - t70, t22 * t32 + t97, 0, 0;];
tauc_reg = t7;
