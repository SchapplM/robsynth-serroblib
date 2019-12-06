% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR3
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
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:32
% EndTime: 2019-12-05 17:51:35
% DurationCPUTime: 0.41s
% Computational Cost: add. (600->89), mult. (1295->149), div. (0->0), fcn. (719->8), ass. (0->78)
t42 = qJD(1) + qJD(3);
t37 = cos(pkin(8)) * pkin(1) + pkin(2);
t33 = t37 * qJD(1);
t50 = sin(qJ(3));
t96 = pkin(1) * sin(pkin(8));
t70 = qJD(3) * t96;
t63 = qJD(1) * t70;
t52 = cos(qJ(3));
t77 = qJD(3) * t52;
t61 = t33 * t77 - t50 * t63;
t13 = t42 * qJD(4) + t61;
t45 = sin(pkin(9));
t40 = t45 ^ 2;
t47 = cos(pkin(9));
t79 = t47 ^ 2 + t40;
t101 = t79 * t13;
t82 = t50 * t33;
t20 = qJD(3) * t82 + t52 * t63;
t55 = t50 * t37 + t52 * t96;
t24 = t55 * qJD(3);
t100 = -t24 * t42 - t20;
t71 = qJD(1) * t96;
t22 = t52 * t71 + t82;
t99 = t22 * t42 - t20;
t98 = -t52 * t37 + t50 * t96;
t21 = t52 * t33 - t50 * t71;
t58 = qJD(4) - t21;
t30 = -t47 * pkin(4) - t45 * pkin(7) - pkin(3);
t3 = t30 * t42 + t58;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t17 = t42 * qJ(4) + t22;
t6 = t45 * qJD(2) + t47 * t17;
t83 = t47 * t51;
t87 = t40 * t51;
t97 = (t13 * t83 + t49 * t20 + (t51 * t3 - t49 * t6) * qJD(5)) * t47 + t13 * t87;
t5 = -t47 * qJD(2) + t45 * t17;
t95 = t5 * t45;
t56 = t37 * t77 - t50 * t70;
t23 = qJD(4) + t56;
t92 = t23 * t42;
t39 = t42 ^ 2;
t90 = t40 * t39;
t89 = t40 * t42;
t88 = t40 * t49;
t86 = t42 * t45;
t85 = t47 * t42;
t84 = t47 * t49;
t34 = -qJD(5) + t85;
t81 = t51 * t34;
t78 = t49 ^ 2 - t51 ^ 2;
t76 = qJD(5) * t49;
t75 = qJD(5) * t51;
t74 = qJD(5) + t34;
t73 = t5 * t86;
t72 = t42 * t87;
t69 = qJD(5) * t89;
t68 = t45 * t76;
t67 = t45 * t75;
t60 = -t49 * t3 - t51 * t6;
t2 = t60 * qJD(5) - t13 * t84 + t51 * t20;
t65 = t13 * t88 - t2 * t47 + t5 * t67;
t64 = t34 * t68;
t62 = t74 * t86;
t59 = t6 * t47 + t95;
t57 = t34 * t47 + t89;
t54 = qJ(4) * t75 + t58 * t49;
t53 = -t34 ^ 2 - t90;
t28 = t68 * t85;
t26 = -0.2e1 * t51 * t49 * t69;
t25 = qJ(4) + t55;
t19 = 0.2e1 * t78 * t69;
t18 = t30 + t98;
t15 = t20 * t45;
t14 = -t42 * pkin(3) + t58;
t12 = (t34 + t85) * t67;
t11 = t28 + t64;
t1 = [0, 0, 0, 0, 0, t100, -t56 * t42 - t61, t100 * t47, t24 * t86 + t15, t79 * t92 + t101, t20 * (-pkin(3) + t98) + t14 * t24 + t59 * t23 + t25 * t101, t26, t19, t11, t12, 0, -(-t23 * t84 + t51 * t24) * t34 + t88 * t92 + (-(-t18 * t49 - t25 * t83) * t34 + t25 * t72) * qJD(5) + t65, (t23 * t83 + t49 * t24) * t34 + t23 * t72 + (t18 * t81 + (-t57 * t25 - t95) * t49) * qJD(5) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t34 - t85) * t67, t28 - t64; 0, 0, 0, 0, 0, t99, t21 * t42 - t61, t99 * t47, -t22 * t86 + t15, t58 * t42 * t79 + t101, -t20 * pkin(3) + qJ(4) * t101 - t14 * t22 + t58 * t59, t26, t19, t11, t12, 0, (t51 * t22 + t30 * t76 + t54 * t47) * t34 + t54 * t89 + t65, -t49 * t22 * t34 + t58 * t57 * t51 + (t30 * t81 + (-t57 * qJ(4) - t95) * t49) * qJD(5) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t39, -t59 * t42 + t20, 0, 0, 0, 0, 0, t53 * t49, t53 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t49 * t87, -t78 * t90, -t49 * t62, -t51 * t62, 0, t60 * t34 - t51 * t73 + t2, (-t47 * t13 - t74 * t3) * t51 + (t74 * t6 - t20 + t73) * t49;];
tauc_reg = t1;
