% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:31
% DurationCPUTime: 0.35s
% Computational Cost: add. (333->50), mult. (640->109), div. (0->0), fcn. (593->6), ass. (0->42)
t28 = sin(pkin(7));
t29 = cos(pkin(8));
t49 = cos(qJ(5));
t39 = t49 * t29;
t27 = sin(pkin(8));
t31 = sin(qJ(5));
t46 = t31 * t27;
t33 = t39 - t46;
t10 = t33 * t28;
t30 = cos(pkin(7));
t42 = t30 * qJD(2);
t19 = -qJD(4) + t42;
t38 = (t27 ^ 2 + t29 ^ 2) * t19;
t22 = t28 * qJD(2);
t52 = 0.2e1 * t22;
t51 = 0.2e1 * qJD(2);
t32 = -pkin(1) - pkin(2);
t44 = t30 * qJ(2) + t28 * t32;
t17 = -qJ(4) + t44;
t50 = pkin(6) - t17;
t40 = t49 * t27;
t45 = t31 * t29;
t16 = t40 + t45;
t13 = t16 * qJD(5);
t48 = t33 * t13;
t36 = qJD(5) * t49;
t12 = qJD(5) * t46 - t29 * t36;
t47 = t16 * t12;
t41 = t31 * t50;
t37 = -t28 * qJ(2) + t30 * t32;
t35 = pkin(3) - t37;
t34 = t50 * t40;
t11 = t29 * pkin(4) + t35;
t9 = t16 * t28;
t8 = t50 * t29;
t7 = qJD(5) * t10;
t6 = t28 * t13;
t4 = t27 * t41 - t49 * t8;
t3 = t31 * t8 + t34;
t2 = t8 * t36 - t19 * t45 + (-qJD(5) * t41 - t49 * t19) * t27;
t1 = -t19 * t39 - qJD(5) * t34 + (-qJD(5) * t8 + t19 * t27) * t31;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, qJ(2) * t51, 0, 0, 0, 0, 0, 0, t52, 0.2e1 * t42, 0, (-t37 * t28 + t44 * t30) * t51, 0, 0, 0, 0, 0, 0, t29 * t52, -0.2e1 * t27 * t22, -0.2e1 * t38, 0.2e1 * t17 * t38 + 0.2e1 * t35 * t22, -0.2e1 * t47, -0.2e1 * t12 * t33 - 0.2e1 * t16 * t13, 0, -0.2e1 * t48, 0, 0, -0.2e1 * t11 * t13 + 0.2e1 * t22 * t33, 0.2e1 * t11 * t12 - 0.2e1 * t16 * t22, 0.2e1 * t1 * t33 - 0.2e1 * t3 * t12 + 0.2e1 * t4 * t13 + 0.2e1 * t2 * t16, -0.2e1 * t4 * t1 + 0.2e1 * t11 * t22 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t38 - t42) * t28, 0, 0, 0, 0, 0, 0, t30 * t13, -t30 * t12, t10 * t13 + t9 * t12 - t7 * t16 + t33 * t6, -t1 * t10 - t2 * t9 - t28 * t42 - t3 * t7 - t4 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10 * t6 + 0.2e1 * t9 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t16 - t4 * t12 - t3 * t13 + t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t12 + t9 * t13 - t6 * t16 - t33 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t47 - 0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, -t13, t12, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t13, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
