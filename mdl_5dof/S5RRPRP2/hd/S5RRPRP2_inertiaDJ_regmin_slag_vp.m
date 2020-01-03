% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:44
% DurationCPUTime: 0.23s
% Computational Cost: add. (222->54), mult. (573->92), div. (0->0), fcn. (381->6), ass. (0->48)
t54 = 2 * qJD(5);
t35 = cos(qJ(2));
t25 = t35 * pkin(1) + pkin(2);
t31 = cos(pkin(8));
t30 = sin(pkin(8));
t33 = sin(qJ(2));
t50 = t30 * t33;
t38 = -pkin(1) * t50 + t31 * t25;
t11 = -pkin(3) - t38;
t47 = pkin(1) * qJD(2);
t49 = t31 * t33;
t13 = (t30 * t35 + t49) * t47;
t34 = cos(qJ(4));
t27 = t34 * qJD(4);
t32 = sin(qJ(4));
t53 = t11 * t27 + t13 * t32;
t52 = t31 * pkin(2);
t46 = t32 * qJD(4);
t16 = -pkin(4) * t46 + qJ(5) * t27 + t32 * qJD(5);
t4 = t13 - t16;
t51 = t16 - t4;
t48 = pkin(1) * t49 + t30 * t25;
t45 = t33 * t47;
t44 = t35 * t47;
t24 = -pkin(3) - t52;
t43 = t24 * t46;
t42 = t24 * t27;
t23 = t30 * pkin(2) + pkin(7);
t41 = t23 * t46;
t40 = t23 * t27;
t39 = t11 * t46 - t13 * t34;
t14 = (t31 * t35 - t50) * t47;
t28 = t32 ^ 2;
t29 = t34 ^ 2;
t3 = (t28 + t29) * t14;
t37 = -t34 * pkin(4) - t32 * qJ(5);
t36 = -pkin(3) + t37;
t15 = t37 * qJD(4) + t34 * qJD(5);
t20 = 0.2e1 * t32 * t27;
t18 = 0.2e1 * (-t28 + t29) * qJD(4);
t17 = t36 - t52;
t12 = pkin(7) + t48;
t10 = t17 * t46;
t6 = t36 - t38;
t5 = t6 * t46;
t2 = t12 * t27 + t32 * t14;
t1 = t12 * t46 - t34 * t14;
t7 = [0, 0, 0, 0, -0.2e1 * t45, -0.2e1 * t44, -0.2e1 * t38 * t13 + 0.2e1 * t48 * t14, t20, t18, 0, 0, 0, 0.2e1 * t39, 0.2e1 * t53, -0.2e1 * t4 * t34 + 0.2e1 * t5, 0.2e1 * t3, -0.2e1 * t6 * t27 - 0.2e1 * t4 * t32, 0.2e1 * t12 * t3 + 0.2e1 * t6 * t4; 0, 0, 0, 0, -t45, -t44, (-t13 * t31 + t14 * t30) * pkin(2), t20, t18, 0, 0, 0, t39 + t43, t42 + t53, t51 * t34 + t10 + t5, t3, t51 * t32 + (-t17 - t6) * t27, -t6 * t16 + t4 * t17 + t23 * t3; 0, 0, 0, 0, 0, 0, 0, t20, t18, 0, 0, 0, 0.2e1 * t43, 0.2e1 * t42, 0.2e1 * t16 * t34 + 0.2e1 * t10, 0, 0.2e1 * t16 * t32 - 0.2e1 * t17 * t27, -0.2e1 * t17 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t46, 0, -t2, t1, -t2, t15, -t1, (-pkin(4) * t32 + qJ(5) * t34) * t14 + t15 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t46, 0, -t40, t41, -t40, t15, -t41, t15 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t27, -t46, 0, t27, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, qJ(5) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
