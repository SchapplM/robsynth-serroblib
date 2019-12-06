% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:48
% EndTime: 2019-12-05 15:14:50
% DurationCPUTime: 0.24s
% Computational Cost: add. (156->40), mult. (447->80), div. (0->0), fcn. (424->8), ass. (0->37)
t21 = sin(pkin(9));
t24 = sin(qJ(3));
t36 = cos(pkin(9));
t39 = cos(qJ(3));
t9 = t24 * t21 - t39 * t36;
t41 = qJD(4) + qJD(5);
t40 = pkin(6) + pkin(7);
t22 = sin(qJ(5));
t23 = sin(qJ(4));
t38 = t22 * t23;
t25 = cos(qJ(5));
t35 = qJD(5) * t25;
t34 = t23 * qJD(4);
t26 = cos(qJ(4));
t33 = t26 * qJD(4);
t32 = -0.2e1 * pkin(3) * qJD(4);
t31 = pkin(4) * t34;
t30 = qJD(5) * t22 * pkin(4);
t29 = pkin(4) * t35;
t28 = qJD(4) * t40;
t12 = t22 * t26 + t25 * t23;
t11 = -t25 * t26 + t38;
t10 = t21 * t39 + t24 * t36;
t6 = t41 * t12;
t20 = -t26 * pkin(4) - pkin(3);
t16 = t40 * t26;
t15 = t40 * t23;
t14 = t26 * t28;
t13 = t23 * t28;
t8 = t10 * qJD(3);
t7 = t9 * qJD(3);
t5 = -t25 * t33 - t26 * t35 + t41 * t38;
t4 = t22 * t13 - t25 * t14 + (t15 * t22 - t16 * t25) * qJD(5);
t3 = t25 * t13 + t22 * t14 + (t15 * t25 + t16 * t22) * qJD(5);
t2 = t41 * t10 * t11 + t12 * t7;
t1 = t10 * t6 - t11 * t7;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, -t8 * t26 + t34 * t9, t8 * t23 + t33 * t9, 0, 0, 0, 0, 0, t8 * t11 + t9 * t6, t8 * t12 - t9 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0.2e1 * t23 * t33, 0.2e1 * (-t23 ^ 2 + t26 ^ 2) * qJD(4), 0, 0, 0, t23 * t32, t26 * t32, -0.2e1 * t12 * t5, 0.2e1 * t5 * t11 - 0.2e1 * t12 * t6, 0, 0, 0, 0.2e1 * t11 * t31 + 0.2e1 * t20 * t6, 0.2e1 * t12 * t31 - 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t33 + t23 * t7, t10 * t34 + t26 * t7, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t33, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, t33, -t34, 0, -pkin(6) * t33, pkin(6) * t34, 0, 0, -t5, -t6, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t30, -0.2e1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
