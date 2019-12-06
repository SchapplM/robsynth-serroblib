% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:28
% EndTime: 2019-12-05 18:16:30
% DurationCPUTime: 0.30s
% Computational Cost: add. (333->63), mult. (803->96), div. (0->0), fcn. (652->8), ass. (0->62)
t42 = cos(pkin(9)) * pkin(1) + pkin(2);
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t75 = pkin(1) * sin(pkin(9));
t78 = -t52 * t42 + t49 * t75;
t77 = qJD(4) + qJD(5);
t76 = -pkin(7) - pkin(8);
t51 = cos(qJ(4));
t74 = t51 * pkin(4);
t53 = t49 * t42 + t52 * t75;
t27 = pkin(7) + t53;
t73 = -pkin(8) - t27;
t47 = sin(qJ(5));
t48 = sin(qJ(4));
t50 = cos(qJ(5));
t32 = t47 * t51 + t50 * t48;
t16 = t77 * t32;
t25 = t53 * qJD(3);
t62 = t48 * qJD(4);
t59 = pkin(4) * t62;
t19 = t25 + t59;
t26 = -pkin(3) + t78;
t23 = t26 - t74;
t70 = t47 * t48;
t31 = -t50 * t51 + t70;
t72 = t23 * t16 + t19 * t31;
t44 = t51 * qJD(4);
t63 = qJD(5) * t50;
t15 = -t50 * t44 - t51 * t63 + t77 * t70;
t71 = -t23 * t15 + t19 * t32;
t24 = t78 * qJD(3);
t69 = t48 * t24;
t68 = t51 * t24;
t66 = t25 * t48 + t26 * t44;
t43 = -pkin(3) - t74;
t65 = t43 * t16 + t31 * t59;
t64 = -t43 * t15 + t32 * t59;
t61 = pkin(3) * t62;
t60 = pkin(3) * t44;
t58 = qJD(5) * t47 * pkin(4);
t57 = pkin(4) * t63;
t56 = qJD(4) * t76;
t55 = -t25 * t51 + t26 * t62;
t54 = qJD(4) * t73;
t45 = t51 * pkin(8);
t38 = 0.2e1 * t48 * t44;
t37 = t51 * pkin(7) + t45;
t36 = t76 * t48;
t34 = t51 * t56;
t33 = t48 * t56;
t30 = 0.2e1 * (-t48 ^ 2 + t51 ^ 2) * qJD(4);
t18 = t51 * t27 + t45;
t17 = t73 * t48;
t10 = -0.2e1 * t32 * t15;
t9 = t51 * t54 + t69;
t8 = t48 * t54 - t68;
t5 = -t47 * t33 + t50 * t34 + (-t36 * t47 - t37 * t50) * qJD(5);
t4 = -t50 * t33 - t47 * t34 + (-t36 * t50 + t37 * t47) * qJD(5);
t3 = 0.2e1 * t15 * t31 - 0.2e1 * t32 * t16;
t2 = -t47 * t8 + t50 * t9 + (-t17 * t47 - t18 * t50) * qJD(5);
t1 = -t47 * t9 - t50 * t8 + (-t17 * t50 + t18 * t47) * qJD(5);
t6 = [0, 0, 0, 0, 0, -0.2e1 * t25, 0.2e1 * t24, t38, t30, 0, 0, 0, 0.2e1 * t55, 0.2e1 * t66, t10, t3, 0, 0, 0, 0.2e1 * t72, 0.2e1 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t25, t24, t38, t30, 0, 0, 0, t55 - t61, -t60 + t66, t10, t3, 0, 0, 0, t65 + t72, t64 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t38, t30, 0, 0, 0, -0.2e1 * t61, -0.2e1 * t60, t10, t3, 0, 0, 0, 0.2e1 * t65, 0.2e1 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t62, 0, -t27 * t44 + t69, t27 * t62 + t68, 0, 0, -t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t44, 0, 0, 0, 0, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t62, 0, -pkin(7) * t44, pkin(7) * t62, 0, 0, -t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t58, -0.2e1 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
