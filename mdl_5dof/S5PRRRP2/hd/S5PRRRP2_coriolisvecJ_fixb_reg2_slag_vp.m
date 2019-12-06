% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:58
% EndTime: 2019-12-05 16:42:02
% DurationCPUTime: 0.61s
% Computational Cost: add. (673->123), mult. (1317->160), div. (0->0), fcn. (621->4), ass. (0->87)
t44 = qJD(2) + qJD(3);
t47 = sin(qJ(4));
t45 = t47 ^ 2;
t49 = cos(qJ(4));
t46 = t49 ^ 2;
t108 = (t45 + t46) * t44;
t84 = t49 * qJD(1);
t48 = sin(qJ(3));
t90 = pkin(2) * qJD(2);
t80 = t48 * t90;
t30 = t44 * pkin(7) + t80;
t96 = t47 * t30;
t14 = t84 - t96;
t106 = qJD(5) - t14;
t50 = cos(qJ(3));
t89 = pkin(2) * qJD(3);
t74 = qJD(2) * t89;
t68 = t50 * t74;
t92 = qJD(4) * t84 + t49 * t68;
t2 = (qJD(5) - t96) * qJD(4) + t92;
t15 = t47 * qJD(1) + t49 * t30;
t62 = t47 * t68;
t7 = t15 * qJD(4) + t62;
t4 = t7 * t47;
t105 = t2 * t49 + t4;
t88 = qJD(4) * t47;
t6 = -t30 * t88 + t92;
t104 = t6 * t49 + t4 + (-t14 * t49 - t15 * t47) * qJD(4);
t63 = t14 * t47 - t15 * t49;
t103 = t63 * t50;
t11 = -qJD(4) * pkin(4) + t106;
t82 = qJD(4) * qJ(5);
t12 = t15 + t82;
t91 = -t45 + t46;
t102 = 0.2e1 * t91 * t44 * qJD(4);
t51 = qJD(4) ^ 2;
t101 = pkin(7) * t51;
t100 = t50 * pkin(2);
t99 = t7 * t49;
t32 = -t49 * pkin(4) - t47 * qJ(5) - pkin(3);
t98 = t32 * t44;
t97 = t44 * t49;
t41 = t51 * t47;
t42 = t51 * t49;
t79 = t50 * t90;
t31 = -t44 * pkin(3) - t79;
t69 = t48 * t74;
t87 = qJD(4) * t49;
t95 = t31 * t87 + t47 * t69;
t81 = t50 * t89;
t94 = t81 * t108;
t86 = qJD(4) * t50;
t76 = t47 * t86;
t93 = t76 * t90 + t80 * t97;
t85 = t47 * qJD(5);
t83 = -qJD(2) - t44;
t78 = t44 * t88;
t77 = t44 * t87;
t65 = pkin(4) * t47 - qJ(5) * t49;
t5 = t69 + (t65 * qJD(4) - t85) * t44;
t75 = -t5 - t101;
t73 = t14 + t96;
t70 = t47 * t77;
t67 = (-qJD(3) + t44) * t90;
t66 = t83 * t89;
t64 = t11 * t47 + t12 * t49;
t19 = pkin(4) * t88 - t49 * t82 - t85;
t13 = t48 * t89 + t19;
t38 = t48 * pkin(2) + pkin(7);
t61 = -t13 * t44 - t38 * t51 - t5;
t23 = t32 - t100;
t60 = t23 * t44 - t81;
t59 = t11 * t87 - t12 * t88 + t105;
t58 = t79 * t108;
t57 = -t44 * t47 * t48 + t49 * t86;
t53 = (t11 * t49 - t12 * t47) * qJD(4) + t105;
t43 = t44 ^ 2;
t39 = -pkin(3) - t100;
t36 = t47 * t43 * t49;
t29 = -0.2e1 * t70;
t28 = 0.2e1 * t70;
t24 = t91 * t43;
t21 = t31 * t88;
t20 = t65 * t44;
t10 = -t79 + t98;
t8 = t10 * t88;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, 0, -t63 * qJD(4) + t6 * t47 - t99, 0, 0, 0, 0, 0, 0, -t41, 0, t42, t64 * qJD(4) + t2 * t47 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t66, t50 * t66, 0, 0, t28, t102, t42, t29, -t41, 0, t39 * t78 - t38 * t42 + t21 + (t83 * t49 * t48 - t76) * t89, t38 * t41 + t39 * t77 - t57 * t89 + t95, t104 + t94, t104 * t38 + (-t103 + (qJD(2) * t39 + t31) * t48) * t89, t28, t42, -t102, 0, t41, t29, t49 * t61 + t60 * t88 + t8, t59 + t94, t61 * t47 + (-t10 - t60) * t87, t10 * t13 + t5 * t23 + t38 * t53 + t64 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t67, t50 * t67, 0, 0, t28, t102, t42, t29, -t41, 0, -pkin(3) * t78 + t21 + (-t69 - t101) * t49 + t93, -pkin(3) * t77 + pkin(7) * t41 + t57 * t90 + t95, -t58 + t104, t104 * pkin(7) + (t103 + (-pkin(3) * qJD(3) - t31) * t48) * t90, t28, t42, -t102, 0, t41, t29, t32 * t78 + t8 + (-t19 * t44 + t75) * t49 + t93, -t58 + t59, (-t10 - t79 - t98) * t87 + ((-t19 + t80) * t44 + t75) * t47, t10 * t19 + t5 * t32 + (-t10 * t48 - t50 * t64) * t90 + t53 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t24, 0, t36, 0, 0, (-t31 * t44 - t68) * t47, t73 * qJD(4) - t31 * t97 - t92, 0, 0, -t36, 0, t24, 0, 0, t36, -t62 + (-t10 * t47 + t20 * t49) * t44, 0, (t10 * t49 + t20 * t47) * t44 + (0.2e1 * qJD(5) - t73) * qJD(4) + t92, -t7 * pkin(4) + t2 * qJ(5) - t10 * t20 + t106 * t12 - t11 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, -t45 * t43 - t51, (t10 * t44 + t68) * t47 + (-t12 + t15) * qJD(4);];
tauc_reg = t1;
