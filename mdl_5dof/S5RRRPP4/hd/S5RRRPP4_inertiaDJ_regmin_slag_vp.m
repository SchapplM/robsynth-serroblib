% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:47
% DurationCPUTime: 0.56s
% Computational Cost: add. (1174->86), mult. (2718->158), div. (0->0), fcn. (2377->6), ass. (0->59)
t75 = qJD(2) + qJD(3);
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t45 = sin(qJ(2));
t74 = pkin(6) + pkin(7);
t64 = t74 * t45;
t47 = cos(qJ(2));
t76 = t74 * t47;
t51 = t44 * t64 - t46 * t76;
t15 = t75 * t51;
t32 = t44 * t45 - t46 * t47;
t22 = t75 * t32;
t33 = t44 * t47 + t46 * t45;
t77 = t22 * qJ(4) - t33 * qJD(4) + t15;
t52 = t44 * t76 + t46 * t64;
t14 = t75 * t52;
t48 = 2 * qJD(5);
t17 = -t32 * qJ(4) - t51;
t43 = sin(pkin(8));
t49 = -t33 * qJ(4) - t52;
t70 = cos(pkin(8));
t10 = t43 * t17 - t70 * t49;
t60 = t70 * t44;
t71 = pkin(2) * qJD(3);
t27 = (t43 * t46 + t60) * t71;
t73 = t10 * t27;
t21 = -t43 * t32 + t70 * t33;
t72 = t27 * t21;
t40 = t46 * pkin(2) + pkin(3);
t30 = pkin(2) * t60 + t43 * t40;
t69 = t45 * qJD(2);
t68 = t47 * qJD(2);
t67 = -0.2e1 * pkin(1) * qJD(2);
t42 = pkin(2) * t69;
t66 = t44 * t71;
t65 = t46 * t71;
t41 = -t47 * pkin(2) - pkin(1);
t11 = t70 * t17 + t43 * t49;
t23 = t75 * t33;
t8 = -t23 * qJ(4) - t32 * qJD(4) - t14;
t3 = t43 * t8 - t70 * t77;
t4 = t77 * t43 + t70 * t8;
t61 = t10 * t3 + t11 * t4;
t18 = t23 * pkin(3) + t42;
t28 = -t43 * t66 + t70 * t65;
t55 = t32 * pkin(3) + t41;
t29 = -t43 * t44 * pkin(2) + t70 * t40;
t12 = -t43 * t22 + t70 * t23;
t13 = -t70 * t22 - t43 * t23;
t20 = t70 * t32 + t43 * t33;
t54 = 0.2e1 * t10 * t13 - 0.2e1 * t11 * t12 - 0.2e1 * t4 * t20 + 0.2e1 * t3 * t21;
t38 = -t70 * pkin(3) - pkin(4);
t37 = t43 * pkin(3) + qJ(5);
t26 = -pkin(4) - t29;
t25 = qJ(5) + t30;
t24 = qJD(5) + t28;
t9 = t20 * pkin(4) - t21 * qJ(5) + t55;
t5 = t12 * pkin(4) - t13 * qJ(5) - t21 * qJD(5) + t18;
t1 = [0, 0, 0, 0.2e1 * t45 * t68, 0.2e1 * (-t45 ^ 2 + t47 ^ 2) * qJD(2), 0, 0, 0, t45 * t67, t47 * t67, -0.2e1 * t33 * t22, 0.2e1 * t22 * t32 - 0.2e1 * t33 * t23, 0, 0, 0, 0.2e1 * t41 * t23 + 0.2e1 * t32 * t42, -0.2e1 * t41 * t22 + 0.2e1 * t33 * t42, t54, 0.2e1 * t18 * t55 + 0.2e1 * t61, 0.2e1 * t9 * t12 + 0.2e1 * t5 * t20, t54, -0.2e1 * t9 * t13 - 0.2e1 * t5 * t21, 0.2e1 * t9 * t5 + 0.2e1 * t61; 0, 0, 0, 0, 0, t68, -t69, 0, -pkin(6) * t68, pkin(6) * t69, 0, 0, -t22, -t23, 0, t15, t14, -t30 * t12 - t29 * t13 - t28 * t20 + t72, t11 * t28 - t3 * t29 + t4 * t30 + t73, -t3, -t25 * t12 + t26 * t13 - t24 * t20 + t72, t4, t11 * t24 + t4 * t25 + t3 * t26 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t65, 0, -0.2e1 * t29 * t27 + 0.2e1 * t30 * t28, -0.2e1 * t27, 0, 0.2e1 * t24, 0.2e1 * t25 * t24 + 0.2e1 * t26 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, t15, t14, (-t12 * t43 - t13 * t70) * pkin(3), (-t3 * t70 + t4 * t43) * pkin(3), -t3, -qJD(5) * t20 - t37 * t12 + t38 * t13, t4, t11 * qJD(5) + t3 * t38 + t4 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, (-t27 * t70 + t28 * t43) * pkin(3), -t27, 0, t28 + t48, t25 * qJD(5) + t24 * t37 + t27 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t37 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t12, 0, -t13, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
