% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.18s
% Computational Cost: add. (184->47), mult. (329->70), div. (0->0), fcn. (198->4), ass. (0->29)
t27 = 2 * qJD(3);
t24 = cos(qJ(4));
t21 = qJD(4) * t24;
t26 = -pkin(2) - pkin(3);
t33 = t24 * qJD(3) + t26 * t21;
t32 = pkin(1) * qJD(2);
t23 = sin(qJ(2));
t17 = t23 * pkin(1) + qJ(3);
t31 = qJ(3) + t17;
t22 = sin(qJ(4));
t20 = qJD(4) * t22;
t25 = cos(qJ(2));
t18 = t25 * t32;
t13 = t18 + qJD(3);
t28 = -t25 * pkin(1) - pkin(2);
t15 = -pkin(3) + t28;
t29 = t23 * t32;
t30 = t24 * t13 + t15 * t21 + t22 * t29;
t4 = t22 * t15 + t24 * t17;
t9 = t24 * qJ(3) + t22 * t26;
t16 = -0.2e1 * t29;
t12 = t24 * t29;
t8 = -t22 * qJ(3) + t24 * t26;
t6 = t22 * qJD(3) + qJD(4) * t9;
t5 = qJ(3) * t20 - t33;
t3 = t24 * t15 - t22 * t17;
t2 = qJD(4) * t4 + t22 * t13 - t12;
t1 = t17 * t20 - t30;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -0.2e1 * t18, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0.2e1 * t13, 0.2e1 * t17 * t13 + 0.2e1 * t28 * t29, 0, 0, 0, 0, 0, 0, 0.2e1 * t2, -0.2e1 * t1, 0, -0.2e1 * t4 * t1 - 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t18, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, t27 + t18, -pkin(2) * t29 + t13 * qJ(3) + t17 * qJD(3), 0, 0, 0, 0, 0, 0, -t12 + (qJD(3) + t13) * t22 + (t31 * t24 + (t15 + t26) * t22) * qJD(4), -t31 * t20 + t30 + t33, 0, -t1 * t9 - t2 * t8 - t3 * t6 - t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, qJ(3) * t27, 0, 0, 0, 0, 0, 0, 0.2e1 * t6, -0.2e1 * t5, 0, -0.2e1 * t9 * t5 - 0.2e1 * t8 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, t20, t21, 0, -t1 * t22 - t2 * t24 + (-t22 * t3 + t24 * t4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, -t5 * t22 - t6 * t24 + (-t22 * t8 + t24 * t9) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
