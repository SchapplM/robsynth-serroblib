% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR3
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:09
% EndTime: 2021-01-15 15:42:12
% DurationCPUTime: 0.29s
% Computational Cost: add. (374->64), mult. (936->126), div. (0->0), fcn. (847->6), ass. (0->42)
t40 = sin(pkin(9));
t41 = cos(pkin(9));
t43 = sin(qJ(3));
t44 = cos(qJ(3));
t28 = t40 * t43 - t41 * t44;
t29 = t40 * t44 + t41 * t43;
t42 = sin(qJ(5));
t52 = cos(qJ(5));
t11 = t52 * t28 + t42 * t29;
t53 = pkin(3) * t40;
t50 = -qJ(4) - pkin(6);
t33 = t50 * t43;
t34 = t50 * t44;
t14 = t40 * t33 - t41 * t34;
t49 = t43 * qJD(3);
t48 = t44 * qJD(3);
t47 = -0.2e1 * pkin(2) * qJD(3);
t39 = pkin(3) * t49;
t38 = -t44 * pkin(3) - pkin(2);
t45 = qJD(3) * t50;
t24 = t44 * qJD(4) + t43 * t45;
t25 = -t43 * qJD(4) + t44 * t45;
t7 = -t40 * t24 + t41 * t25;
t13 = t41 * t33 + t40 * t34;
t8 = t41 * t24 + t40 * t25;
t12 = -t42 * t28 + t52 * t29;
t37 = t41 * pkin(3) + pkin(4);
t27 = -t40 * t49 + t41 * t48;
t26 = t29 * qJD(3);
t19 = (-t37 * t42 - t52 * t53) * qJD(5);
t18 = (-t52 * t37 + t42 * t53) * qJD(5);
t16 = t28 * pkin(4) + t38;
t15 = t26 * pkin(4) + t39;
t10 = -t28 * pkin(7) + t14;
t9 = -t29 * pkin(7) + t13;
t6 = -t26 * pkin(7) + t8;
t5 = -t27 * pkin(7) + t7;
t4 = t12 * qJD(5) + t52 * t26 + t42 * t27;
t3 = t11 * qJD(5) + t42 * t26 - t52 * t27;
t2 = t52 * t5 - t42 * t6 + (-t52 * t10 - t42 * t9) * qJD(5);
t1 = -t52 * t6 - t42 * t5 + (t10 * t42 - t52 * t9) * qJD(5);
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t28 * t26 + 0.2e1 * t29 * t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t13 + t27 * t14 - t28 * t7 + t29 * t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0.2e1 * t43 * t48, 0.2e1 * (-t43 ^ 2 + t44 ^ 2) * qJD(3), 0, 0, 0, t43 * t47, t44 * t47, 0.2e1 * t38 * t26 + 0.2e1 * t28 * t39, 0.2e1 * t38 * t27 + 0.2e1 * t29 * t39, -0.2e1 * t13 * t27 - 0.2e1 * t14 * t26 - 0.2e1 * t8 * t28 - 0.2e1 * t7 * t29, 0.2e1 * t13 * t7 + 0.2e1 * t14 * t8 + 0.2e1 * t38 * t39, -0.2e1 * t12 * t3, 0.2e1 * t3 * t11 - 0.2e1 * t12 * t4, 0, 0, 0, 0.2e1 * t15 * t11 + 0.2e1 * t16 * t4, 0.2e1 * t15 * t12 - 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, -t26, -t27, 0, (-t26 * t41 + t27 * t40) * pkin(3), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, t48, -t49, 0, -pkin(6) * t48, pkin(6) * t49, t7, -t8, (-t26 * t40 - t27 * t41) * pkin(3), (t40 * t8 + t41 * t7) * pkin(3), 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, t39, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
