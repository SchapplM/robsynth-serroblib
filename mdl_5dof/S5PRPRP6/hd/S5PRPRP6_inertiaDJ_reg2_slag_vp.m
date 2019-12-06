% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:19
% EndTime: 2019-12-05 15:41:21
% DurationCPUTime: 0.37s
% Computational Cost: add. (85->45), mult. (252->72), div. (0->0), fcn. (162->4), ass. (0->35)
t23 = sin(qJ(2));
t17 = t23 * qJD(2);
t22 = sin(qJ(4));
t20 = t22 ^ 2;
t24 = cos(qJ(4));
t21 = t24 ^ 2;
t35 = t20 + t21;
t4 = t35 * t17;
t39 = 2 * qJD(3);
t38 = 2 * qJD(5);
t26 = -pkin(2) - pkin(6);
t37 = t26 * t4;
t25 = cos(qJ(2));
t33 = t25 * qJD(2);
t36 = qJ(3) * t33 + t23 * qJD(3);
t34 = t22 * qJD(4);
t18 = t24 * qJD(4);
t32 = qJ(3) * qJD(4);
t30 = t26 * t34;
t29 = t22 * t18;
t28 = t26 * t18;
t27 = -t22 * pkin(4) + t24 * qJ(5);
t3 = t27 * qJD(4) + t22 * qJD(5);
t19 = qJ(3) * t39;
t14 = -0.2e1 * t29;
t13 = 0.2e1 * t29;
t10 = (t20 - t21) * qJD(4);
t9 = qJ(3) - t27;
t8 = -t23 * t34 + t24 * t33;
t7 = t23 * t18 + t22 * t33;
t6 = -t22 * t17 + t25 * t18;
t5 = t24 * t17 + t25 * t34;
t2 = -t24 * qJD(5) + qJD(3) + (pkin(4) * t24 + qJ(5) * t22) * qJD(4);
t1 = 0.2e1 * (0.1e1 - t35) * t23 * t33;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t33, -pkin(2) * t17 + t36, 0, 0, 0, 0, 0, 0, t7, t8, -t4, t36 + t37, 0, 0, 0, 0, 0, 0, t7, -t4, -t8, t23 * t2 + t9 * t33 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t19, t14, 0.2e1 * t10, 0, t13, 0, 0, 0.2e1 * qJD(3) * t22 + 0.2e1 * t24 * t32, 0.2e1 * qJD(3) * t24 - 0.2e1 * t22 * t32, 0, t19, t14, 0, -0.2e1 * t10, 0, 0, t13, 0.2e1 * t9 * t18 + 0.2e1 * t2 * t22, 0, -0.2e1 * t2 * t24 + 0.2e1 * t9 * t34, 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, (-qJ(5) * qJD(4) * t25 + pkin(4) * t17) * t24 + (qJ(5) * t17 + (qJD(4) * pkin(4) - qJD(5)) * t25) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t18, 0, -t30, -t28, 0, 0, 0, -t34, 0, 0, t18, 0, -t30, -t3, t28, t3 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t18, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, qJ(5) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
