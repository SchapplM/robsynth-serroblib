% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:42
% DurationCPUTime: 0.21s
% Computational Cost: add. (203->47), mult. (545->86), div. (0->0), fcn. (500->8), ass. (0->45)
t53 = qJD(2) + qJD(3);
t31 = sin(pkin(9));
t34 = sin(qJ(3));
t52 = t31 * t34;
t32 = cos(pkin(9));
t51 = t32 * t34;
t35 = sin(qJ(2));
t50 = t34 * t35;
t37 = cos(qJ(3));
t29 = t37 * pkin(2) + pkin(3);
t17 = -pkin(2) * t52 + t32 * t29;
t13 = -pkin(4) - t17;
t48 = pkin(2) * qJD(3);
t15 = (t31 * t37 + t51) * t48;
t36 = cos(qJ(5));
t30 = qJD(5) * t36;
t33 = sin(qJ(5));
t49 = t13 * t30 + t15 * t33;
t18 = pkin(2) * t51 + t31 * t29;
t38 = cos(qJ(2));
t47 = qJD(2) * t38;
t46 = qJD(3) * t37;
t45 = qJD(5) * t33;
t44 = t34 * t48;
t43 = pkin(2) * t46;
t28 = -t32 * pkin(3) - pkin(4);
t42 = t28 * t45;
t41 = t28 * t30;
t40 = t13 * t45 - t15 * t36;
t21 = t34 * t38 + t37 * t35;
t39 = t37 * t38 - t50;
t27 = t31 * pkin(3) + pkin(7);
t23 = 0.2e1 * t33 * t30;
t20 = 0.2e1 * (-t33 ^ 2 + t36 ^ 2) * qJD(5);
t16 = (t32 * t37 - t52) * t48;
t14 = pkin(7) + t18;
t9 = t53 * t21;
t8 = -t37 * t47 - t38 * t46 + t53 * t50;
t7 = t32 * t21 + t31 * t39;
t6 = t31 * t21 - t32 * t39;
t4 = -t31 * t9 - t32 * t8;
t3 = -t31 * t8 + t32 * t9;
t2 = -t3 * t36 + t6 * t45;
t1 = t3 * t33 + t6 * t30;
t5 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * t6 * t3 + 0.2e1 * t7 * t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t35, -t47, 0, -t9, t8, t6 * t15 + t7 * t16 - t3 * t17 + t4 * t18, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -0.2e1 * t44, -0.2e1 * t43, -0.2e1 * t17 * t15 + 0.2e1 * t18 * t16, t23, t20, 0, 0, 0, 0.2e1 * t40, 0.2e1 * t49; 0, 0, 0, 0, 0, -t9, t8, (-t3 * t32 + t31 * t4) * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -t44, -t43, (-t15 * t32 + t16 * t31) * pkin(3), t23, t20, 0, 0, 0, t40 + t42, t41 + t49; 0, 0, 0, 0, 0, 0, 0, 0, t23, t20, 0, 0, 0, 0.2e1 * t42, 0.2e1 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t30 - t33 * t4, -t36 * t4 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t45, 0, -t14 * t30 - t33 * t16, t14 * t45 - t36 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t45, 0, -t27 * t30, t27 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
