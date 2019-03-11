% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:34
% EndTime: 2019-03-09 02:39:37
% DurationCPUTime: 0.87s
% Computational Cost: add. (1377->125), mult. (3117->243), div. (0->0), fcn. (2956->10), ass. (0->86)
t66 = sin(pkin(11));
t68 = cos(pkin(11));
t70 = sin(qJ(6));
t72 = cos(qJ(6));
t116 = -t70 * t66 + t72 * t68;
t115 = t116 * qJD(6);
t100 = t66 ^ 2 + t68 ^ 2;
t67 = sin(pkin(10));
t73 = cos(qJ(3));
t98 = cos(pkin(10));
t90 = t98 * t73;
t71 = sin(qJ(3));
t97 = t71 * qJD(3);
t44 = qJD(3) * t90 - t67 * t97;
t114 = t44 * t100;
t113 = 0.2e1 * t115;
t57 = t67 * pkin(3) + qJ(5);
t112 = pkin(8) + t57;
t60 = sin(pkin(9)) * pkin(1) + pkin(7);
t99 = qJ(4) + t60;
t88 = qJD(3) * t99;
t40 = t73 * qJD(4) - t71 * t88;
t74 = -t71 * qJD(4) - t73 * t88;
t27 = t67 * t40 - t98 * t74;
t49 = t99 * t73;
t91 = t98 * t71;
t36 = t67 * t49 + t99 * t91;
t111 = t36 * t27;
t52 = t67 * t73 + t91;
t43 = t52 * qJD(3);
t110 = t43 * t66;
t109 = t43 * t68;
t105 = t67 * t71;
t50 = -t90 + t105;
t108 = t50 * t43;
t107 = t52 * t66;
t106 = t66 * t44;
t104 = t68 * t44;
t63 = pkin(3) * t97;
t20 = t43 * pkin(4) - t44 * qJ(5) - t52 * qJD(5) + t63;
t28 = t98 * t40 + t67 * t74;
t7 = t66 * t20 + t68 * t28;
t62 = -cos(pkin(9)) * pkin(1) - pkin(2);
t76 = -t73 * pkin(3) + t62;
t32 = t50 * pkin(4) - t52 * qJ(5) + t76;
t37 = -t99 * t105 + t98 * t49;
t14 = t66 * t32 + t68 * t37;
t96 = t73 * qJD(3);
t95 = 0.2e1 * t96;
t6 = t68 * t20 - t66 * t28;
t13 = t68 * t32 - t66 * t37;
t89 = t100 * qJD(5);
t87 = 0.2e1 * t89;
t5 = -t68 * t52 * pkin(8) + t50 * pkin(5) + t13;
t8 = -pkin(8) * t107 + t14;
t86 = t72 * t5 - t70 * t8;
t85 = t70 * t5 + t72 * t8;
t84 = t6 * t68 + t7 * t66;
t83 = -t6 * t66 + t7 * t68;
t61 = -t98 * pkin(3) - pkin(4);
t82 = -t13 * t66 + t14 * t68;
t81 = t27 * t50 + t36 * t43;
t80 = t27 * t52 + t36 * t44;
t53 = t72 * t66 + t70 * t68;
t18 = t115 * t50 + t53 * t43;
t46 = t53 * qJD(6);
t79 = -t116 * t43 + t50 * t46;
t47 = t112 * t66;
t48 = t112 * t68;
t78 = -t72 * t47 - t70 * t48;
t77 = -t70 * t47 + t72 * t48;
t75 = -qJD(5) * t50 - t43 * t57 + t44 * t61;
t54 = -t68 * pkin(5) + t61;
t35 = t116 * t52;
t34 = t53 * t52;
t26 = pkin(5) * t107 + t36;
t24 = -t53 * qJD(5) - qJD(6) * t77;
t23 = -qJD(5) * t116 - t78 * qJD(6);
t15 = pkin(5) * t106 + t27;
t12 = t115 * t52 + t44 * t53;
t11 = -t116 * t44 + t46 * t52;
t4 = -pkin(8) * t106 + t7;
t3 = t43 * pkin(5) - pkin(8) * t104 + t6;
t2 = -t85 * qJD(6) + t72 * t3 - t70 * t4;
t1 = -t86 * qJD(6) - t70 * t3 - t72 * t4;
t9 = [0, 0, 0, 0, t71 * t95, 0.2e1 * (-t71 ^ 2 + t73 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t62 * t97, t62 * t95, -0.2e1 * t28 * t50 - 0.2e1 * t37 * t43 + 0.2e1 * t80, 0.2e1 * t37 * t28 + 0.2e1 * t76 * t63 + 0.2e1 * t111, 0.2e1 * t13 * t43 + 0.2e1 * t6 * t50 + 0.2e1 * t80 * t66, -0.2e1 * t14 * t43 - 0.2e1 * t7 * t50 + 0.2e1 * t80 * t68, -0.2e1 * t84 * t52 + 0.2e1 * (-t13 * t68 - t14 * t66) * t44, 0.2e1 * t13 * t6 + 0.2e1 * t14 * t7 + 0.2e1 * t111, -0.2e1 * t35 * t11, 0.2e1 * t11 * t34 - 0.2e1 * t35 * t12, -0.2e1 * t11 * t50 + 0.2e1 * t35 * t43, -0.2e1 * t50 * t12 - 0.2e1 * t43 * t34, 0.2e1 * t108, 0.2e1 * t26 * t12 + 0.2e1 * t15 * t34 + 0.2e1 * t2 * t50 + 0.2e1 * t86 * t43, 0.2e1 * t1 * t50 - 0.2e1 * t26 * t11 + 0.2e1 * t15 * t35 - 0.2e1 * t85 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t52 + t37 * t44 + t81, 0, 0, 0, t82 * t44 + t83 * t52 + t81, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52 * t44 + 0.2e1 * t108, 0, 0, 0, 0.2e1 * t114 * t52 + 0.2e1 * t108, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t96, -t97, 0, -t60 * t96, t60 * t97 (-t43 * t67 - t98 * t44) * pkin(3) (-t98 * t27 + t28 * t67) * pkin(3), -t27 * t68 + t66 * t75, t27 * t66 + t68 * t75, t83, t82 * qJD(5) + t27 * t61 + t83 * t57, -t11 * t53 + t115 * t35, -t11 * t116 - t115 * t34 - t53 * t12 - t35 * t46, t18, -t79, 0, -t116 * t15 + t54 * t12 + t24 * t50 + t26 * t46 + t78 * t43, -t54 * t11 + t115 * t26 + t15 * t53 + t23 * t50 - t43 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0 (-t98 * t43 + t44 * t67) * pkin(3), -t109, t110, t114, t114 * t57 + t43 * t61 + t52 * t89, 0, 0, 0, 0, 0, t79, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t57 * t87, t53 * t113, 0.2e1 * t115 * t116 - 0.2e1 * t53 * t46, 0, 0, 0, 0.2e1 * t54 * t46, t54 * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t109, -t110, -t114, t84, 0, 0, 0, 0, 0, -t79, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t104, 0, t27, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, t43, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t46, 0, t24, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
