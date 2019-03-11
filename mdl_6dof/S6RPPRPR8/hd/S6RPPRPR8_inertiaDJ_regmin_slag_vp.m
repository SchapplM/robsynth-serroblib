% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:15
% EndTime: 2019-03-09 01:56:16
% DurationCPUTime: 0.54s
% Computational Cost: add. (645->97), mult. (1389->181), div. (0->0), fcn. (1297->6), ass. (0->69)
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t82 = sin(qJ(4));
t66 = qJD(4) * t82;
t83 = cos(qJ(4));
t67 = qJD(4) * t83;
t23 = -t42 * t67 - t43 * t66;
t27 = t82 * t42 - t83 * t43;
t80 = t27 * t23;
t24 = -t42 * t66 + t43 * t67;
t28 = t83 * t42 + t82 * t43;
t81 = t24 * t28;
t89 = 0.2e1 * t80 - 0.2e1 * t81;
t45 = sin(qJ(6));
t40 = t45 ^ 2;
t46 = cos(qJ(6));
t75 = -t46 ^ 2 + t40;
t65 = t75 * qJD(6);
t31 = (t42 ^ 2 + t43 ^ 2) * qJD(3);
t44 = -pkin(1) - qJ(3);
t84 = -pkin(7) + t44;
t29 = t84 * t42;
t30 = t84 * t43;
t8 = t28 * qJD(3) + t29 * t66 - t30 * t67;
t88 = 2 * qJD(2);
t87 = 2 * qJD(5);
t86 = pkin(4) + pkin(8);
t4 = -t24 * pkin(5) - t8;
t85 = t4 * t28;
t79 = t45 * t23;
t78 = t46 * t23;
t77 = t46 * t24;
t35 = t42 * pkin(3) + qJ(2);
t74 = qJD(6) * t45;
t73 = qJD(6) * t46;
t72 = qJD(6) * t86;
t71 = qJ(2) * qJD(2);
t70 = qJ(5) * qJD(6);
t69 = t45 * t77;
t68 = t45 * t73;
t26 = t28 ^ 2;
t64 = qJD(6) * (t27 ^ 2 + t26);
t63 = t27 * qJ(5) + t35;
t19 = t82 * t29 - t83 * t30;
t10 = -t27 * pkin(5) + t19;
t7 = t86 * t28 + t63;
t62 = t46 * t10 - t45 * t7;
t61 = t45 * t10 + t46 * t7;
t59 = qJ(5) * t24 + qJD(5) * t28;
t57 = t27 * t73 - t79;
t14 = t27 * t74 + t78;
t16 = t45 * t24 + t28 * t73;
t56 = t28 * t74 - t77;
t55 = -t23 * qJ(5) + t27 * qJD(5) + qJD(2);
t52 = t83 * t29 + t82 * t30;
t51 = pkin(4) * t23 + t59;
t50 = t4 + (qJ(5) * t28 - t27 * t86) * qJD(6);
t9 = -t27 * qJD(3) + t52 * qJD(4);
t49 = t19 * t23 - t24 * t52 - t9 * t27 + t8 * t28;
t11 = -t28 * pkin(5) + t52;
t48 = -qJD(6) * t11 + t23 * t86 + t59;
t21 = -0.2e1 * t80;
t18 = t28 * pkin(4) + t63;
t6 = t24 * pkin(4) + t55;
t5 = t23 * pkin(5) + t9;
t3 = t86 * t24 + t55;
t2 = -t61 * qJD(6) - t45 * t3 + t46 * t5;
t1 = -t62 * qJD(6) - t46 * t3 - t45 * t5;
t12 = [0, 0, 0, 0, t88, 0.2e1 * t71, t42 * t88, t43 * t88, 0.2e1 * t31, -0.2e1 * t44 * t31 + 0.2e1 * t71, t21, -0.2e1 * t23 * t28 + 0.2e1 * t27 * t24, 0, 0, 0, 0.2e1 * qJD(2) * t28 + 0.2e1 * t35 * t24, -0.2e1 * qJD(2) * t27 + 0.2e1 * t35 * t23, 0.2e1 * t49, -0.2e1 * t18 * t24 - 0.2e1 * t6 * t28, -0.2e1 * t18 * t23 + 0.2e1 * t6 * t27, 0.2e1 * t18 * t6 + 0.2e1 * t19 * t9 - 0.2e1 * t52 * t8, 0.2e1 * t26 * t68 + 0.2e1 * t40 * t81, -0.2e1 * t26 * t65 + 0.4e1 * t28 * t69, -0.2e1 * t16 * t27 + 0.2e1 * t28 * t79, 0.2e1 * t56 * t27 + 0.2e1 * t28 * t78, t21, 0.2e1 * t56 * t11 - 0.2e1 * t2 * t27 + 0.2e1 * t62 * t23 - 0.2e1 * t46 * t85, -0.2e1 * t1 * t27 + 0.2e1 * t16 * t11 - 0.2e1 * t61 * t23 + 0.2e1 * t45 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, -t49, 0, 0, 0, 0, 0, t45 * t64 + t46 * t89, -t45 * t89 + t46 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t24, t23, 0, -t24, -t23, t6, 0, 0, 0, 0, 0, t57, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, -t9, t8, -t51, t9, -t8, -t9 * pkin(4) - t8 * qJ(5) + qJD(5) * t52, -t28 * t65 + t69, -t75 * t24 - 0.4e1 * t28 * t68, t14, t57, 0, t50 * t45 - t48 * t46, t48 * t45 + t50 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, -t23, t24, t51, 0, 0, 0, 0, 0, t16, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, qJ(5) * t87, -0.2e1 * t68, 0.2e1 * t65, 0, 0, 0, 0.2e1 * qJD(5) * t45 + 0.2e1 * t46 * t70, 0.2e1 * qJD(5) * t46 - 0.2e1 * t45 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, t9, 0, 0, 0, 0, 0, t14, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t56, t23, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73, 0, t45 * t72, t46 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
