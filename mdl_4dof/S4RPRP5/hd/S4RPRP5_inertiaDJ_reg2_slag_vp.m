% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:02
% DurationCPUTime: 0.20s
% Computational Cost: add. (200->38), mult. (524->76), div. (0->0), fcn. (453->4), ass. (0->33)
t24 = sin(pkin(6));
t25 = cos(pkin(6));
t39 = cos(qJ(3));
t32 = qJD(3) * t39;
t26 = sin(qJ(3));
t36 = qJD(3) * t26;
t11 = t24 * t36 - t25 * t32;
t42 = -0.2e1 * t11;
t38 = t26 * t25;
t16 = t39 * t24 + t38;
t12 = t16 * qJD(3);
t41 = 0.2e1 * t12;
t40 = 2 * qJD(4);
t37 = pkin(5) + qJ(2);
t15 = t26 * t24 - t39 * t25;
t35 = t15 * t41;
t17 = t37 * t25;
t33 = t37 * t24;
t10 = t39 * t17 - t26 * t33;
t28 = t39 * t33;
t31 = t39 * qJD(2);
t4 = qJD(3) * t28 - t25 * t31 + (qJD(2) * t24 + qJD(3) * t17) * t26;
t5 = t17 * t32 + qJD(2) * t38 + (-t37 * t36 + t31) * t24;
t9 = t26 * t17 + t28;
t34 = -t10 * t4 + t9 * t5;
t21 = -t25 * pkin(2) - pkin(1);
t30 = 0.2e1 * (t24 ^ 2 + t25 ^ 2) * qJD(2);
t29 = t11 * t15 - t16 * t12;
t27 = -0.2e1 * t10 * t12 - 0.2e1 * t9 * t11 + 0.2e1 * t4 * t15 + 0.2e1 * t5 * t16;
t8 = t16 * t42;
t7 = t15 * pkin(3) - t16 * qJ(4) + t21;
t3 = t12 * pkin(3) + t11 * qJ(4) - t16 * qJD(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, qJ(2) * t30, t8, 0.2e1 * t29, 0, t35, 0, 0, t21 * t41, t21 * t42, t27, 0.2e1 * t34, t8, 0, -0.2e1 * t29, 0, 0, t35, 0.2e1 * t7 * t12 + 0.2e1 * t3 * t15, t27, 0.2e1 * t7 * t11 - 0.2e1 * t3 * t16, 0.2e1 * t7 * t3 + 0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t11, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, -t5, t4, 0, 0, 0, -t11, 0, 0, t12, 0, -t5, pkin(3) * t11 - t12 * qJ(4) - t15 * qJD(4), -t4, -t5 * pkin(3) - t4 * qJ(4) + t10 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, qJ(4) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
