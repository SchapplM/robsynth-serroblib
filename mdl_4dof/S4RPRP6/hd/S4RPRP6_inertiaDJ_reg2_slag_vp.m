% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (74->29), mult. (158->52), div. (0->0), fcn. (88->2), ass. (0->22)
t21 = 2 * qJD(2);
t14 = -pkin(1) - pkin(5);
t20 = qJ(4) - t14;
t12 = sin(qJ(3));
t19 = t12 * qJD(3);
t13 = cos(qJ(3));
t18 = t13 * qJD(3);
t17 = qJ(2) * qJD(3);
t16 = pkin(3) * t19;
t15 = t12 * t18;
t6 = t20 * t13;
t11 = qJ(2) * t21;
t10 = t12 * pkin(3) + qJ(2);
t9 = pkin(3) * t18 + qJD(2);
t8 = -0.2e1 * t15;
t7 = 0.2e1 * t15;
t5 = t20 * t12;
t4 = 0.2e1 * (t12 ^ 2 - t13 ^ 2) * qJD(3);
t3 = -qJD(3) * t6 - t12 * qJD(4);
t2 = -t13 * qJD(4) + t20 * t19;
t1 = t3 * t12 + t2 * t13 + (t12 * t6 - t13 * t5) * qJD(3);
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t11, t8, t4, 0, t7, 0, 0, 0.2e1 * qJD(2) * t12 + 0.2e1 * t13 * t17, 0.2e1 * qJD(2) * t13 - 0.2e1 * t12 * t17, 0, t11, t8, t4, 0, t7, 0, 0, 0.2e1 * t10 * t18 + 0.2e1 * t9 * t12, -0.2e1 * t10 * t19 + 0.2e1 * t9 * t13, -0.2e1 * t1, 0.2e1 * t10 * t9 - 0.2e1 * t6 * t2 - 0.2e1 * t5 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, -t18, 0, -t14 * t19, -t14 * t18, 0, 0, 0, 0, -t19, 0, -t18, 0, t2, -t3, t16, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t22;
