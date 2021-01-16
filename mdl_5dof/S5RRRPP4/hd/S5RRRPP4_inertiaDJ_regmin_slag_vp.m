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
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 22:24:49
% EndTime: 2021-01-15 22:24:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (1336->97), mult. (3116->169), div. (0->0), fcn. (2713->6), ass. (0->60)
t83 = qJD(2) + qJD(3);
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t46 = sin(qJ(2));
t79 = pkin(6) + pkin(7);
t66 = t79 * t46;
t48 = cos(qJ(2));
t80 = t79 * t48;
t53 = t45 * t66 - t47 * t80;
t15 = t83 * t53;
t34 = t45 * t46 - t47 * t48;
t22 = t83 * t34;
t58 = t45 * t48 + t47 * t46;
t82 = t22 * qJ(4) - t58 * qJD(4) + t15;
t14 = t83 * (t45 * t80 + t47 * t66);
t49 = 2 * qJD(5);
t17 = -t34 * qJ(4) - t53;
t44 = sin(pkin(8));
t50 = t58 * (-qJ(4) - t79);
t73 = cos(pkin(8));
t10 = t44 * t17 - t73 * t50;
t62 = t73 * t45;
t74 = pkin(2) * qJD(3);
t29 = (t44 * t47 + t62) * t74;
t78 = t10 * t29;
t21 = -t44 * t34 + t73 * t58;
t77 = t29 * t21;
t41 = t47 * pkin(2) + pkin(3);
t32 = pkin(2) * t62 + t44 * t41;
t72 = t46 * qJD(2);
t71 = t48 * qJD(2);
t70 = -0.2e1 * pkin(1) * qJD(2);
t43 = pkin(2) * t72;
t68 = t45 * t74;
t67 = t47 * t74;
t42 = -t48 * pkin(2) - pkin(1);
t11 = t73 * t17 + t44 * t50;
t57 = t58 * qJD(2);
t23 = t58 * qJD(3) + t57;
t8 = -t23 * qJ(4) - t34 * qJD(4) - t14;
t3 = t44 * t8 - t73 * t82;
t4 = t82 * t44 + t73 * t8;
t63 = t10 * t3 + t11 * t4;
t18 = t23 * pkin(3) + t43;
t25 = t34 * pkin(3) + t42;
t30 = -t44 * t68 + t73 * t67;
t31 = -t44 * t45 * pkin(2) + t73 * t41;
t12 = -t44 * t22 + t73 * t23;
t13 = -t73 * t22 - t44 * t23;
t20 = t73 * t34 + t44 * t58;
t56 = 0.2e1 * t10 * t13 - 0.2e1 * t11 * t12 - 0.2e1 * t4 * t20 + 0.2e1 * t3 * t21;
t39 = -t73 * pkin(3) - pkin(4);
t38 = t44 * pkin(3) + qJ(5);
t28 = -pkin(4) - t31;
t27 = qJ(5) + t32;
t26 = qJD(5) + t30;
t24 = 0.2e1 * t29;
t9 = t20 * pkin(4) - t21 * qJ(5) + t25;
t5 = t12 * pkin(4) - t13 * qJ(5) - t21 * qJD(5) + t18;
t1 = [0, 0, 0, 0.2e1 * t46 * t71, 0.2e1 * (-t46 ^ 2 + t48 ^ 2) * qJD(2), 0, 0, 0, t46 * t70, t48 * t70, -0.2e1 * t58 * t22, 0.2e1 * t22 * t34 - 0.2e1 * t58 * t23, 0, 0, 0, 0.2e1 * t42 * t23 + 0.2e1 * t34 * t43, 0.2e1 * t46 * pkin(2) * t57 - 0.2e1 * t42 * t22, 0.2e1 * t25 * t12 + 0.2e1 * t18 * t20, 0.2e1 * t25 * t13 + 0.2e1 * t18 * t21, t56, 0.2e1 * t25 * t18 + 0.2e1 * t63, 0.2e1 * t9 * t12 + 0.2e1 * t5 * t20, t56, -0.2e1 * t9 * t13 - 0.2e1 * t5 * t21, 0.2e1 * t9 * t5 + 0.2e1 * t63; 0, 0, 0, 0, 0, t71, -t72, 0, -pkin(6) * t71, pkin(6) * t72, 0, 0, -t22, -t23, 0, t15, t14, -t3, -t4, -t32 * t12 - t31 * t13 - t30 * t20 + t77, t11 * t30 - t3 * t31 + t4 * t32 + t78, -t3, -t27 * t12 + t28 * t13 - t26 * t20 + t77, t4, t11 * t26 + t4 * t27 + t3 * t28 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t68, -0.2e1 * t67, -t24, -0.2e1 * t30, 0, -0.2e1 * t31 * t29 + 0.2e1 * t32 * t30, -t24, 0, 0.2e1 * t26, 0.2e1 * t27 * t26 + 0.2e1 * t28 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, t15, t14, -t3, -t4, (-t12 * t44 - t73 * t13) * pkin(3), (-t73 * t3 + t4 * t44) * pkin(3), -t3, -qJD(5) * t20 - t38 * t12 + t39 * t13, t4, t11 * qJD(5) + t3 * t39 + t4 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, -t29, -t30, 0, (-t73 * t29 + t30 * t44) * pkin(3), -t29, 0, t49 + t30, t27 * qJD(5) + t26 * t38 + t29 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t38 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, t18, t12, 0, -t13, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
