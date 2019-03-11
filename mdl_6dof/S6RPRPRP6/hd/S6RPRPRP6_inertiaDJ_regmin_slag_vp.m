% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:45
% EndTime: 2019-03-09 03:19:46
% DurationCPUTime: 0.72s
% Computational Cost: add. (1336->142), mult. (2885->245), div. (0->0), fcn. (2741->6), ass. (0->85)
t102 = pkin(3) + pkin(8);
t101 = cos(qJ(3));
t59 = cos(pkin(9));
t80 = t101 * t59;
t58 = sin(pkin(9));
t61 = sin(qJ(3));
t97 = t61 * t58;
t42 = -t80 + t97;
t43 = t101 * t58 + t61 * t59;
t53 = -t59 * pkin(2) - pkin(1);
t71 = -t43 * qJ(4) + t53;
t18 = t102 * t42 + t71;
t95 = pkin(7) + qJ(2);
t46 = t95 * t58;
t47 = t95 * t59;
t30 = t101 * t46 + t61 * t47;
t24 = t43 * pkin(4) + t30;
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t94 = t62 * t18 + t60 * t24;
t31 = -t101 * t47 + t61 * t46;
t56 = t60 ^ 2;
t57 = t62 ^ 2;
t77 = (t56 - t57) * qJD(5);
t79 = qJD(3) * t101;
t20 = (qJD(2) * t58 + qJD(3) * t47) * t61 - qJD(2) * t80 + t46 * t79;
t35 = qJD(3) * t97 - t59 * t79;
t104 = -0.2e1 * t35;
t103 = 2 * qJD(4);
t36 = t43 * qJD(3);
t13 = -t36 * pkin(4) - t20;
t100 = t13 * t42;
t99 = t56 * t36;
t98 = t60 * t36;
t96 = t62 * t36;
t91 = qJ(6) * t42;
t90 = qJ(6) + t102;
t89 = qJD(5) * t60;
t88 = qJD(5) * t62;
t87 = qJD(5) * t102;
t86 = t62 * qJD(6);
t85 = qJ(4) * qJD(5);
t84 = t60 * t96;
t83 = pkin(5) * t89;
t82 = pkin(5) * t88;
t81 = t60 * t88;
t45 = t90 * t62;
t78 = -t18 - t91;
t76 = 0.2e1 * (t58 ^ 2 + t59 ^ 2) * qJD(2);
t23 = t62 * t24;
t5 = t43 * pkin(5) + t78 * t60 + t23;
t6 = t62 * t91 + t94;
t75 = -t5 * t60 + t6 * t62;
t74 = t35 * qJ(4) - t43 * qJD(4);
t73 = -qJ(4) * t36 - qJD(4) * t42;
t70 = -t60 * t35 + t43 * t88;
t69 = t62 * t35 + t43 * t89;
t68 = t42 * t88 + t98;
t67 = t42 * t89 - t96;
t10 = t102 * t36 + t74;
t21 = t43 * qJD(2) - t31 * qJD(3);
t14 = -t35 * pkin(4) + t21;
t3 = -t62 * t10 - t60 * t14 + t18 * t89 - t24 * t88;
t66 = t13 + (qJ(4) * t42 + t102 * t43) * qJD(5);
t25 = -t42 * pkin(4) - t31;
t65 = -qJD(5) * t25 - t102 * t35 - t73;
t12 = t62 * t14;
t1 = -t35 * pkin(5) + t12 + t78 * t88 + (-qJ(6) * t36 - qJD(5) * t24 - qJD(6) * t42 - t10) * t60;
t2 = -t67 * qJ(6) + t42 * t86 - t3;
t64 = -t1 * t60 + t2 * t62 + (-t5 * t62 - t6 * t60) * qJD(5);
t52 = t60 * pkin(5) + qJ(4);
t48 = qJD(4) + t82;
t44 = t90 * t60;
t38 = t42 ^ 2;
t34 = -qJD(5) * t45 - t60 * qJD(6);
t33 = t90 * t89 - t86;
t32 = t57 * t36;
t29 = t43 * t104;
t28 = t42 * pkin(3) + t71;
t17 = t36 * pkin(3) + t74;
t15 = (-pkin(5) * t62 - pkin(4)) * t42 - t31;
t8 = t33 * t62 + t34 * t60 + (-t44 * t62 + t45 * t60) * qJD(5);
t7 = t67 * pkin(5) + t13;
t4 = -t94 * qJD(5) - t60 * t10 + t12;
t9 = [0, 0, 0, 0, 0, t76, qJ(2) * t76, t29, 0.2e1 * t35 * t42 - 0.2e1 * t43 * t36, 0, 0, 0, 0.2e1 * t53 * t36, t53 * t104, 0.2e1 * t20 * t42 + 0.2e1 * t21 * t43 - 0.2e1 * t30 * t35 + 0.2e1 * t31 * t36, -0.2e1 * t17 * t42 - 0.2e1 * t28 * t36, -0.2e1 * t17 * t43 + 0.2e1 * t28 * t35, 0.2e1 * t28 * t17 + 0.2e1 * t31 * t20 + 0.2e1 * t30 * t21, 0.2e1 * t38 * t81 + 0.2e1 * t42 * t99, -0.2e1 * t38 * t77 + 0.4e1 * t42 * t84, 0.2e1 * t70 * t42 + 0.2e1 * t43 * t98, -0.2e1 * t69 * t42 + 0.2e1 * t43 * t96, t29, 0.2e1 * t4 * t43 - 0.2e1 * (-t60 * t18 + t23) * t35 - 0.2e1 * t62 * t100 + 0.2e1 * t67 * t25, 0.2e1 * t60 * t100 + 0.2e1 * t68 * t25 + 0.2e1 * t3 * t43 + 0.2e1 * t94 * t35, 0.2e1 * t75 * t36 + 0.2e1 * t64 * t42, 0.2e1 * t5 * t1 + 0.2e1 * t15 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, 0, -t36, t35, t17, 0, 0, 0, 0, 0, -t70, t69, t32 + t99, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t21, t20, pkin(3) * t35 + t73, t21, -t20, -t21 * pkin(3) - t20 * qJ(4) - t31 * qJD(4), -t42 * t77 + t84, -0.4e1 * t42 * t81 + t32 - t99, -t69, -t70, 0, t66 * t60 - t65 * t62, t65 * t60 + t66 * t62 (t34 * t42 - t36 * t44 - t1 + (t42 * t45 - t6) * qJD(5)) * t62 + (-t33 * t42 + t36 * t45 - t2 + (t42 * t44 + t5) * qJD(5)) * t60, -t1 * t45 + t15 * t48 - t2 * t44 + t5 * t33 + t6 * t34 + t7 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * t33 + t62 * t34 + (t44 * t60 + t45 * t62) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, qJ(4) * t103, -0.2e1 * t81, 0.2e1 * t77, 0, 0, 0, 0.2e1 * qJD(4) * t60 + 0.2e1 * t62 * t85, 0.2e1 * qJD(4) * t62 - 0.2e1 * t60 * t85, -0.2e1 * t8, -0.2e1 * t45 * t33 - 0.2e1 * t44 * t34 + 0.2e1 * t52 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, t21, 0, 0, 0, 0, 0, -t69, -t70, 0, t75 * qJD(5) + t1 * t62 + t2 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t67, -t35, t4, t3, -t68 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t89, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0, t60 * t87, t62 * t87, t83, t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
