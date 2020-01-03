% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (141->49), mult. (287->76), div. (0->0), fcn. (163->2), ass. (0->36)
t21 = sin(qJ(3));
t23 = -pkin(3) - pkin(4);
t22 = cos(qJ(3));
t33 = t22 * qJ(4);
t36 = t21 * t23 + t33;
t9 = 0.2e1 * (t21 ^ 2 - t22 ^ 2) * qJD(3);
t35 = 2 * qJD(2);
t25 = 2 * qJD(4);
t34 = qJ(4) * t21;
t24 = -pkin(1) - pkin(6);
t32 = qJ(5) + t24;
t17 = t21 * qJD(3);
t31 = t21 * qJD(4);
t18 = t22 * qJD(3);
t30 = qJ(2) * qJD(3);
t29 = t24 * t17;
t28 = t21 * t18;
t10 = t32 * t21;
t27 = -t22 * qJD(4) + qJD(2);
t26 = -t21 * pkin(3) + t33;
t7 = t26 * qJD(3) + t31;
t20 = qJ(2) * t35;
t19 = qJ(4) * t25;
t16 = t24 * t18;
t15 = -0.2e1 * t28;
t14 = 0.2e1 * t28;
t12 = qJ(2) - t26;
t11 = t32 * t22;
t8 = -qJ(2) + t36;
t6 = (pkin(3) * t22 + t34) * qJD(3) + t27;
t5 = t36 * qJD(3) + t31;
t4 = qJ(5) * t18 + t21 * qJD(5) + t16;
t3 = qJD(3) * t10 - t22 * qJD(5);
t2 = (t23 * t22 - t34) * qJD(3) - t27;
t1 = t4 * t21 - t3 * t22 + (t10 * t22 - t11 * t21) * qJD(3);
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t20, t15, t9, 0, t14, 0, 0, 0.2e1 * qJD(2) * t21 + 0.2e1 * t22 * t30, 0.2e1 * qJD(2) * t22 - 0.2e1 * t21 * t30, 0, t20, t15, 0, -t9, 0, 0, t14, 0.2e1 * t12 * t18 + 0.2e1 * t6 * t21, 0, 0.2e1 * t12 * t17 - 0.2e1 * t6 * t22, 0.2e1 * t12 * t6, t15, -t9, 0, t14, 0, 0, -0.2e1 * t8 * t18 - 0.2e1 * t2 * t21, -0.2e1 * t8 * t17 + 0.2e1 * t2 * t22, 0.2e1 * t1, 0.2e1 * t10 * t4 - 0.2e1 * t11 * t3 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t18, 0, -t29, -t16, 0, 0, 0, -t17, 0, 0, t18, 0, -t29, -t7, t16, t7 * t24, 0, 0, t17, 0, -t18, 0, -t3, t4, t5, t4 * qJ(4) + t10 * qJD(4) + t3 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, t18, t7, 0, 0, 0, 0, 0, 0, -t17, t18, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t19, 0, 0, 0, 0, 0, 0, 0, t25, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, t17, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
