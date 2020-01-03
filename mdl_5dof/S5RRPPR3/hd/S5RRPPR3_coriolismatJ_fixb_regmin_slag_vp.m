% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:37
% EndTime: 2019-12-31 19:26:39
% DurationCPUTime: 0.46s
% Computational Cost: add. (319->91), mult. (679->105), div. (0->0), fcn. (525->6), ass. (0->69)
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t31 = t42 ^ 2 - t44 ^ 2;
t39 = qJD(1) + qJD(2);
t84 = t39 * t31;
t41 = cos(pkin(8));
t43 = sin(qJ(2));
t76 = t41 * t43;
t45 = cos(qJ(2));
t79 = t45 * pkin(1);
t36 = pkin(2) + t79;
t40 = sin(pkin(8));
t78 = t40 * t36;
t48 = pkin(1) * t76 + t78;
t18 = qJ(4) + t48;
t77 = t40 * t45;
t22 = (t76 + t77) * pkin(1);
t80 = t40 * pkin(2);
t34 = qJ(4) + t80;
t83 = t18 + t22 + t34;
t32 = t40 * t43 * pkin(1);
t23 = t41 * t79 - t32;
t37 = qJD(4) * t42;
t63 = t42 * qJD(2);
t75 = t23 * t63 + t37;
t38 = qJD(4) * t44;
t61 = t44 * qJD(2);
t74 = t23 * t61 + t38;
t73 = pkin(1) * qJD(1);
t72 = pkin(1) * qJD(2);
t53 = t41 * t36 - t32;
t51 = -pkin(3) - t53;
t1 = t18 * t23 + t51 * t22;
t71 = t1 * qJD(1);
t2 = -t53 * t22 + t48 * t23;
t70 = t2 * qJD(1);
t69 = qJD(5) * t42;
t68 = qJD(5) * t44;
t67 = t18 * qJD(1);
t66 = t22 * qJD(1);
t65 = t23 * qJD(1);
t64 = t42 * qJD(1);
t62 = t44 * qJD(1);
t60 = t23 * qJD(2) + qJD(4);
t59 = t18 * t64;
t58 = t18 * t62;
t57 = t23 * t64;
t56 = t23 * t62;
t55 = t42 * t68;
t54 = -t41 * pkin(2) - pkin(3);
t52 = pkin(1) * t39;
t28 = t39 * t44;
t50 = t22 / 0.2e1 - t34 / 0.2e1 - t18 / 0.2e1;
t8 = -qJ(4) + (t79 / 0.2e1 - pkin(2) / 0.2e1 - t36 / 0.2e1) * t40;
t49 = t8 * qJD(1) - t34 * qJD(2);
t3 = t50 * t42;
t47 = -t3 * qJD(1) + t34 * t63;
t4 = t50 * t44;
t46 = t4 * qJD(1) - t34 * t61;
t33 = -pkin(7) + t54;
t29 = t31 * qJD(5);
t27 = t39 * t42;
t21 = t42 * t28;
t19 = t22 * qJD(2);
t17 = -pkin(7) + t51;
t7 = t80 / 0.2e1 + qJ(4) + t78 / 0.2e1 + (t76 + t77 / 0.2e1) * pkin(1);
t6 = t83 * t44 / 0.2e1;
t5 = -t83 * t42 / 0.2e1;
t9 = [0, 0, 0, 0, -t43 * t72, -t45 * t72, t2 * qJD(2), t19, t60, t1 * qJD(2) + t18 * qJD(4), -t55, t29, 0, 0, 0, t18 * t68 + t75, -t18 * t69 + t74; 0, 0, 0, 0, -t43 * t52, -t45 * t52, t70 + (-t22 * t41 + t23 * t40) * qJD(2) * pkin(2), t19 + t66, t60 + t65, t71 + (t22 * t54 + t23 * t34) * qJD(2) + t7 * qJD(4), -t55, t29, 0, 0, 0, t6 * qJD(5) + t57 + t75, t5 * qJD(5) + t56 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t39, t7 * qJD(2) + t67, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t84, -t69, -t68, 0, t6 * qJD(2) - t17 * t69 + t58, t5 * qJD(2) - t17 * t68 - t59; 0, 0, 0, 0, t43 * t73, t45 * t73, -t70, -t66, qJD(4) - t65, -t8 * qJD(4) - t71, -t55, t29, 0, 0, 0, -t4 * qJD(5) + t37 - t57, t3 * qJD(5) + t38 - t56; 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t34 * qJD(4), -t55, t29, 0, 0, 0, t34 * t68 + t37, -t34 * t69 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t39, -t49, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t84, -t69, -t68, 0, -t33 * t69 - t46, -t33 * t68 - t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t69; 0, 0, 0, 0, 0, 0, 0, 0, -t39, t8 * qJD(2) - t67, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, -t39, t49, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t84, 0, 0, 0, t4 * qJD(2) - t58, -t3 * qJD(2) + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t84, 0, 0, 0, t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
