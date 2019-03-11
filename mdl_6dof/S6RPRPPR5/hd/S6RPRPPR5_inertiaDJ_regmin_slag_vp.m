% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:37
% EndTime: 2019-03-09 02:51:40
% DurationCPUTime: 0.82s
% Computational Cost: add. (1345->137), mult. (3137->253), div. (0->0), fcn. (3059->8), ass. (0->86)
t70 = sin(pkin(10));
t72 = cos(pkin(10));
t75 = sin(qJ(6));
t77 = cos(qJ(6));
t112 = -t75 * t70 + t77 * t72;
t44 = t112 * qJD(6);
t107 = cos(qJ(3));
t71 = sin(pkin(9));
t99 = pkin(7) + qJ(2);
t58 = t99 * t71;
t73 = cos(pkin(9));
t59 = t99 * t73;
t76 = sin(qJ(3));
t40 = -t107 * t59 + t76 * t58;
t97 = t70 ^ 2 + t72 ^ 2;
t57 = t97 * qJD(5);
t93 = qJD(3) * t107;
t94 = t107 * t73;
t31 = (qJD(2) * t71 + qJD(3) * t59) * t76 - qJD(2) * t94 + t58 * t93;
t102 = t76 * t71;
t45 = qJD(3) * t102 - t73 * t93;
t111 = -0.2e1 * t45;
t110 = 2 * qJD(4);
t52 = -t94 + t102;
t109 = pkin(8) * t52;
t100 = pkin(3) + qJ(5);
t108 = -pkin(8) - t100;
t106 = t70 * t45;
t105 = t72 * t45;
t54 = t107 * t71 + t76 * t73;
t46 = t54 * qJD(3);
t104 = t72 * t46;
t83 = t45 * qJ(4) - t54 * qJD(4);
t15 = t52 * qJD(5) + t100 * t46 + t83;
t32 = t54 * qJD(2) - t40 * qJD(3);
t20 = -t45 * pkin(4) + t32;
t7 = t72 * t15 + t70 * t20;
t65 = -t73 * pkin(2) - pkin(1);
t80 = -t54 * qJ(4) + t65;
t26 = t100 * t52 + t80;
t39 = t107 * t58 + t76 * t59;
t33 = t54 * pkin(4) + t39;
t11 = t72 * t26 + t70 * t33;
t96 = qJ(4) * qJD(4);
t95 = -pkin(5) * t72 - pkin(4);
t92 = 0.2e1 * (t71 ^ 2 + t73 ^ 2) * qJD(2);
t18 = t72 * t20;
t6 = -t70 * t15 + t18;
t3 = t6 * t72 + t7 * t70;
t91 = -t6 * t70 + t7 * t72;
t30 = t72 * t33;
t8 = t54 * pkin(5) + t30 + (-t26 - t109) * t70;
t9 = t72 * t109 + t11;
t90 = t75 * t9 - t77 * t8;
t89 = t75 * t8 + t77 * t9;
t19 = -t46 * pkin(4) - t31;
t34 = -t52 * pkin(4) - t40;
t88 = t19 * t52 + t34 * t46;
t51 = t77 * t70 + t75 * t72;
t43 = t51 * qJD(6);
t87 = t112 * t45 + t43 * t54;
t22 = -t44 * t54 + t51 * t45;
t55 = t108 * t70;
t56 = t108 * t72;
t86 = t77 * t55 + t75 * t56;
t85 = t75 * t55 - t77 * t56;
t82 = -qJ(4) * t46 - qJD(4) * t52;
t78 = qJD(5) * t54 - t100 * t45 - t82;
t63 = t70 * pkin(5) + qJ(4);
t38 = t54 * t111;
t37 = t52 * pkin(3) + t80;
t36 = t51 * t52;
t35 = t112 * t52;
t28 = -qJD(5) * t112 - t86 * qJD(6);
t27 = t51 * qJD(5) + t85 * qJD(6);
t25 = t46 * pkin(3) + t83;
t23 = t95 * t52 - t40;
t16 = t95 * t46 - t31;
t14 = -t112 * t46 + t52 * t43;
t13 = t52 * t44 + t51 * t46;
t10 = -t70 * t26 + t30;
t5 = pkin(8) * t104 + t7;
t4 = -t45 * pkin(5) + t18 + (-pkin(8) * t46 - t15) * t70;
t2 = -t89 * qJD(6) + t77 * t4 - t75 * t5;
t1 = t90 * qJD(6) - t75 * t4 - t77 * t5;
t12 = [0, 0, 0, 0, 0, t92, qJ(2) * t92, t38, 0.2e1 * t45 * t52 - 0.2e1 * t54 * t46, 0, 0, 0, 0.2e1 * t65 * t46, t65 * t111, 0.2e1 * t31 * t52 + 0.2e1 * t32 * t54 - 0.2e1 * t39 * t45 + 0.2e1 * t40 * t46, -0.2e1 * t25 * t52 - 0.2e1 * t37 * t46, -0.2e1 * t25 * t54 + 0.2e1 * t37 * t45, 0.2e1 * t37 * t25 + 0.2e1 * t40 * t31 + 0.2e1 * t39 * t32, -0.2e1 * t10 * t45 + 0.2e1 * t6 * t54 - 0.2e1 * t88 * t72, 0.2e1 * t11 * t45 - 0.2e1 * t7 * t54 + 0.2e1 * t88 * t70, 0.2e1 * t91 * t52 + 0.2e1 * (-t10 * t70 + t11 * t72) * t46, 0.2e1 * t10 * t6 + 0.2e1 * t11 * t7 + 0.2e1 * t34 * t19, 0.2e1 * t36 * t13, 0.2e1 * t13 * t35 - 0.2e1 * t36 * t14, 0.2e1 * t13 * t54 - 0.2e1 * t36 * t45, -0.2e1 * t14 * t54 - 0.2e1 * t35 * t45, t38, 0.2e1 * t23 * t14 - 0.2e1 * t16 * t35 + 0.2e1 * t2 * t54 + 0.2e1 * t90 * t45, 0.2e1 * t1 * t54 + 0.2e1 * t23 * t13 + 0.2e1 * t16 * t36 + 0.2e1 * t89 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, -t46, t45, t25, t106, t105, t97 * t46, t91, 0, 0, 0, 0, 0, t22, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, -t32, t31, pkin(3) * t45 + t82, t32, -t31, -t32 * pkin(3) - t31 * qJ(4) - t40 * qJD(4), t19 * t70 - t78 * t72, t19 * t72 + t78 * t70, -t3, t19 * qJ(4) + t34 * qJD(4) - t3 * t100 + (-t10 * t72 - t11 * t70) * qJD(5), t112 * t13 - t36 * t43, -t112 * t14 - t13 * t51 - t43 * t35 - t36 * t44, -t87, t22, 0, -qJD(4) * t35 + t63 * t14 + t16 * t51 + t23 * t44 + t28 * t54 + t85 * t45, qJD(4) * t36 + t112 * t16 + t63 * t13 - t23 * t43 + t27 * t54 + t86 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0.2e1 * t96, t70 * t110, t72 * t110, 0.2e1 * t57, 0.2e1 * t100 * t57 + 0.2e1 * t96, -0.2e1 * t112 * t43, -0.2e1 * t112 * t44 + 0.2e1 * t43 * t51, 0, 0, 0, 0.2e1 * qJD(4) * t51 + 0.2e1 * t63 * t44, 0.2e1 * qJD(4) * t112 - 0.2e1 * t63 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, t32, -t105, t106, 0, t3, 0, 0, 0, 0, 0, -t87, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t70 * t46, 0, t19, 0, 0, 0, 0, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), 0, 0, 0, 0, 0, t44, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, -t45, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
