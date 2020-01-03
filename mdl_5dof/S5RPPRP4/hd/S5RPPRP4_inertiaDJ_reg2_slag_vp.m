% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:21
% DurationCPUTime: 0.35s
% Computational Cost: add. (178->45), mult. (338->83), div. (0->0), fcn. (242->4), ass. (0->37)
t25 = sin(qJ(4));
t20 = t25 * qJD(4);
t26 = cos(qJ(4));
t24 = cos(pkin(7));
t37 = t24 * qJD(2);
t23 = sin(pkin(7));
t27 = -pkin(1) - pkin(2);
t40 = t24 * qJ(2) + t23 * t27;
t8 = -pkin(6) + t40;
t30 = t8 * t20 - t26 * t37;
t1 = -qJ(5) * t20 + t26 * qJD(5) + t30;
t36 = t26 * qJD(4);
t41 = qJ(5) - t8;
t2 = (qJD(5) - t37) * t25 + t41 * t36;
t3 = t41 * t25;
t4 = t41 * t26;
t43 = (-t25 * t4 + t26 * t3) * qJD(4) + t1 * t26 + t2 * t25;
t42 = 0.2e1 * qJD(2);
t21 = t25 ^ 2;
t22 = t26 ^ 2;
t39 = t21 + t22;
t38 = t23 * qJD(2);
t35 = 0.2e1 * t37;
t34 = pkin(4) * t20;
t33 = t25 * t36;
t32 = t23 * t36;
t31 = -t23 * qJ(2) + t24 * t27;
t7 = pkin(3) - t31;
t15 = t24 * t36;
t13 = t24 * t20;
t12 = t23 * t20;
t11 = -0.2e1 * t33;
t10 = 0.2e1 * t33;
t9 = -t34 + t38;
t6 = 0.2e1 * (-t21 + t22) * qJD(4);
t5 = t26 * pkin(4) + t7;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, qJ(2) * t42, 0, 0, 0, 0, 0, 0, 0.2e1 * t38, t35, 0, (-t31 * t23 + t40 * t24) * t42, t10, t6, 0, t11, 0, 0, -0.2e1 * t7 * t20 + 0.2e1 * t26 * t38, -0.2e1 * t25 * t38 - 0.2e1 * t7 * t36, -t39 * t35, (t39 * t8 * t24 + t23 * t7) * t42, t10, t6, 0, t11, 0, 0, -0.2e1 * t5 * t20 + 0.2e1 * t9 * t26, -0.2e1 * t9 * t25 - 0.2e1 * t5 * t36, 0.2e1 * t43, 0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, (-0.1e1 + t39) * t23 * t37, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t43 * t23 - t9 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t25 + t2 * t26 + (-t25 * t3 - t26 * t4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, t20, 0, -t25 * t37 - t8 * t36, t30, 0, 0, 0, 0, -t36, 0, t20, 0, t2, t1, pkin(4) * t36, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t12, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t12, 0, -pkin(4) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t36, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t36, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t36, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
