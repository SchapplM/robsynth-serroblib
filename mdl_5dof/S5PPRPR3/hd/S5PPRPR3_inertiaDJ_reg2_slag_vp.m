% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:24
% DurationCPUTime: 0.35s
% Computational Cost: add. (189->46), mult. (634->94), div. (0->0), fcn. (654->8), ass. (0->40)
t44 = 2 * qJD(5);
t20 = sin(pkin(9));
t22 = cos(pkin(9));
t21 = sin(pkin(8));
t27 = cos(qJ(3));
t38 = t27 * qJD(3);
t34 = t21 * t38;
t25 = sin(qJ(3));
t40 = t25 * qJD(3);
t35 = t21 * t40;
t6 = t20 * t35 - t22 * t34;
t13 = t20 * t27 + t22 * t25;
t8 = t13 * t21;
t43 = t8 * t6;
t10 = t13 * qJD(3);
t12 = t20 * t25 - t22 * t27;
t42 = t12 * t10;
t24 = sin(qJ(5));
t41 = t24 * qJD(5);
t26 = cos(qJ(5));
t39 = t26 * qJD(5);
t17 = -t22 * pkin(3) - pkin(4);
t37 = t17 * t44;
t36 = t24 * t39;
t11 = t12 * qJD(3);
t18 = t24 ^ 2;
t19 = t26 ^ 2;
t33 = (t18 + t19) * t11;
t32 = t8 * t10 - t6 * t12;
t23 = cos(pkin(8));
t9 = t12 * t21;
t29 = t23 * t24 + t26 * t9;
t30 = t23 * t26 - t24 * t9;
t31 = t24 * t30 - t26 * t29;
t7 = t21 * t10;
t1 = t30 * qJD(5) + t26 * t7;
t2 = t29 * qJD(5) + t24 * t7;
t28 = -t1 * t26 - t2 * t24 + (t24 * t29 + t26 * t30) * qJD(5);
t16 = t20 * pkin(3) + pkin(6);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t9 * t7 - 0.2e1 * t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t1 * t29 - 0.2e1 * t2 * t30 - 0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t11 - t7 * t13 + t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t11 + t28 * t13 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13 * t11 + 0.2e1 * t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13 * t33 + 0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t35, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, (-t20 * t7 + t22 * t6) * pkin(3), 0, 0, 0, 0, 0, 0, t6 * t26 + t8 * t41, -t6 * t24 + t8 * t39, t28, t28 * t16 - t6 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t38, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t11, 0, (-t10 * t22 - t11 * t20) * pkin(3), 0, 0, 0, 0, 0, 0, -t10 * t26 + t12 * t41, t10 * t24 + t12 * t39, -t33, t10 * t17 - t16 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t36, (-t18 + t19) * t44, 0, -0.2e1 * t36, 0, 0, t24 * t37, t26 * t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * qJD(5) - t1 * t24 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t11 - t13 * t39, t26 * t11 + t13 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t41, 0, -t16 * t39, t16 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
