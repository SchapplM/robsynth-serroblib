% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:29
% EndTime: 2018-11-14 13:53:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (175->34), mult. (478->67), div. (0->0), fcn. (342->6), ass. (0->34)
t22 = sin(pkin(7));
t39 = pkin(2) * t22;
t23 = cos(pkin(7));
t25 = sin(qJ(2));
t38 = t23 * t25;
t37 = pkin(1) * qJD(2);
t24 = sin(qJ(4));
t36 = qJD(4) * t24;
t26 = cos(qJ(4));
t35 = qJD(4) * t26;
t27 = cos(qJ(2));
t21 = t27 * pkin(1) + pkin(2);
t28 = -t22 * t25 * pkin(1) + t23 * t21;
t11 = pkin(3) + t28;
t12 = (-t22 * t27 - t38) * t37;
t31 = t27 * t37;
t32 = t25 * t37;
t13 = -t22 * t32 + t23 * t31;
t34 = -t11 * t35 - t24 * t12 - t26 * t13;
t33 = t24 * t39;
t15 = pkin(1) * t38 + t22 * t21;
t30 = t15 + t39;
t29 = t26 * t12 - t24 * t13;
t4 = t24 * t11 + t26 * t15;
t20 = t23 * pkin(2) + pkin(3);
t16 = t24 * t20 + t26 * t39;
t17 = t20 * t35;
t14 = t26 * t20 - t33;
t10 = t16 * qJD(4);
t9 = qJD(4) * t33 - t17;
t3 = t26 * t11 - t24 * t15;
t2 = -t4 * qJD(4) + t29;
t1 = t15 * t36 + t34;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t32, -0.2e1 * t31, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12, -0.2e1 * t13, 0, 0.2e1 * t28 * t12 + 0.2e1 * t15 * t13, 0, 0, 0, 0, 0, 0, 0.2e1 * t2, 0.2e1 * t1, 0, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, 0 (t12 * t23 + t13 * t22) * pkin(2), 0, 0, 0, 0, 0, 0 (-t30 * t26 + (-t11 - t20) * t24) * qJD(4) + t29, t30 * t36 - t17 + t34, 0, -t1 * t16 - t3 * t10 + t2 * t14 - t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10, 0.2e1 * t9, 0, -0.2e1 * t14 * t10 - 0.2e1 * t16 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
