% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:16
% DurationCPUTime: 0.32s
% Computational Cost: add. (134->46), mult. (442->80), div. (0->0), fcn. (241->4), ass. (0->41)
t31 = sin(qJ(3));
t29 = t31 ^ 2;
t33 = cos(qJ(3));
t30 = t33 ^ 2;
t34 = cos(qJ(2));
t46 = pkin(1) * qJD(2);
t42 = t34 * t46;
t55 = (t29 + t30) * t42;
t53 = 0.2e1 * (-t29 + t30) * qJD(3);
t52 = 2 * qJD(4);
t32 = sin(qJ(2));
t43 = t32 * t46;
t27 = t31 * qJD(3);
t28 = t33 * qJD(3);
t8 = pkin(3) * t27 - qJ(4) * t28 - t31 * qJD(4);
t1 = t8 + t43;
t51 = -t1 - t8;
t50 = t34 * pkin(1);
t24 = t32 * pkin(1) + pkin(6);
t49 = t55 * t24;
t48 = t55 * pkin(6);
t25 = -pkin(2) - t50;
t47 = t25 * t28 + t31 * t43;
t45 = pkin(2) * t27;
t44 = pkin(2) * t28;
t41 = pkin(6) * t27;
t40 = pkin(6) * t28;
t39 = t31 * t28;
t36 = -t33 * pkin(3) - t31 * qJ(4);
t35 = t25 * t27 - t33 * t43;
t17 = -pkin(2) + t36;
t7 = qJD(3) * t36 + t33 * qJD(4);
t23 = -0.2e1 * t39;
t22 = 0.2e1 * t39;
t10 = t17 - t50;
t9 = t17 * t27;
t5 = t10 * t27;
t4 = t24 * t28 + t31 * t42;
t3 = t24 * t27 - t33 * t42;
t2 = 0.2e1 * t55;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t43, -0.2e1 * t42, 0, 0, t22, t53, 0, t23, 0, 0, 0.2e1 * t35, 0.2e1 * t47, t2, 0.2e1 * t25 * t43 + 0.2e1 * t49, t22, 0, -t53, 0, 0, t23, -0.2e1 * t1 * t33 + 0.2e1 * t5, t2, -0.2e1 * t1 * t31 - 0.2e1 * t10 * t28, 0.2e1 * t10 * t1 + 0.2e1 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, 0, 0, t22, t53, 0, t23, 0, 0, t35 - t45, -t44 + t47, t55, -pkin(2) * t43 + t48, t22, 0, -t53, 0, 0, t23, t51 * t33 + t5 + t9, t55, t51 * t31 + (-t10 - t17) * t28, t1 * t17 + t10 * t8 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t53, 0, t23, 0, 0, -0.2e1 * t45, -0.2e1 * t44, 0, 0, t22, 0, -t53, 0, 0, t23, -0.2e1 * t8 * t33 + 0.2e1 * t9, 0, -0.2e1 * t17 * t28 - 0.2e1 * t8 * t31, 0.2e1 * t17 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t27, 0, -t4, t3, 0, 0, 0, t28, 0, 0, t27, 0, -t4, t7, -t3, (-pkin(3) * t31 + qJ(4) * t33) * t42 + t7 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t27, 0, -t40, t41, 0, 0, 0, t28, 0, 0, t27, 0, -t40, t7, -t41, t7 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, qJ(4) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
