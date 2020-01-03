% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:13
% DurationCPUTime: 0.19s
% Computational Cost: add. (98->24), mult. (289->47), div. (0->0), fcn. (182->6), ass. (0->26)
t20 = sin(qJ(3));
t24 = cos(pkin(7)) * pkin(1) + pkin(2);
t32 = cos(qJ(3));
t33 = sin(pkin(7)) * pkin(1);
t22 = -t20 * t33 + t32 * t24;
t21 = cos(qJ(4));
t15 = t21 * qJD(4);
t19 = sin(qJ(4));
t31 = t20 * t24 + t32 * t33;
t6 = t31 * qJD(3);
t7 = -pkin(3) - t22;
t34 = t7 * t15 + t6 * t19;
t30 = t19 * qJD(4);
t28 = pkin(3) * t30;
t27 = pkin(3) * t15;
t26 = t19 * t15;
t25 = -t6 * t21 + t7 * t30;
t16 = t19 ^ 2;
t17 = t21 ^ 2;
t5 = t22 * qJD(3);
t1 = (-t16 - t17) * t5;
t13 = -0.2e1 * t26;
t12 = 0.2e1 * t26;
t9 = 0.2e1 * (-t16 + t17) * qJD(4);
t8 = pkin(6) + t31;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t6, -0.2e1 * t5, 0, -0.2e1 * t22 * t6 + 0.2e1 * t31 * t5, t12, t9, 0, t13, 0, 0, 0.2e1 * t25, 0.2e1 * t34, -0.2e1 * t1, -0.2e1 * t8 * t1 + 0.2e1 * t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5, 0, 0, t12, t9, 0, t13, 0, 0, t25 - t28, -t27 + t34, -t1, -t6 * pkin(3) - pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t9, 0, t13, 0, 0, -0.2e1 * t28, -0.2e1 * t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t30, 0, -t8 * t15 - t19 * t5, -t21 * t5 + t8 * t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t30, 0, -pkin(6) * t15, pkin(6) * t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
