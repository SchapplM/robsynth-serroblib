% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:02
% DurationCPUTime: 0.13s
% Computational Cost: add. (147->37), mult. (380->75), div. (0->0), fcn. (325->4), ass. (0->27)
t19 = sin(pkin(6));
t20 = cos(pkin(6));
t30 = cos(qJ(3));
t25 = qJD(3) * t30;
t21 = sin(qJ(3));
t27 = qJD(3) * t21;
t7 = t19 * t27 - t20 * t25;
t32 = -0.2e1 * t7;
t31 = 2 * qJD(4);
t29 = t21 * t20;
t28 = pkin(5) + qJ(2);
t16 = -t20 * pkin(2) - pkin(1);
t26 = t28 * t19;
t24 = t30 * qJD(2);
t23 = 0.2e1 * (t19 ^ 2 + t20 ^ 2) * qJD(2);
t22 = t30 * t26;
t11 = t30 * t19 + t29;
t12 = t28 * t20;
t10 = t21 * t19 - t30 * t20;
t8 = t11 * qJD(3);
t6 = t30 * t12 - t21 * t26;
t5 = t21 * t12 + t22;
t4 = t10 * pkin(3) - t11 * qJ(4) + t16;
t3 = t12 * t25 + qJD(2) * t29 + (-t28 * t27 + t24) * t19;
t2 = -t20 * t24 + qJD(3) * t22 + (qJD(2) * t19 + qJD(3) * t12) * t21;
t1 = t8 * pkin(3) + t7 * qJ(4) - t11 * qJD(4);
t9 = [0, 0, 0, 0, 0, t23, qJ(2) * t23, t11 * t32, 0.2e1 * t7 * t10 - 0.2e1 * t11 * t8, 0, 0, 0, 0.2e1 * t16 * t8, t16 * t32, 0.2e1 * t1 * t10 + 0.2e1 * t4 * t8, 0.2e1 * t2 * t10 + 0.2e1 * t3 * t11 - 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t1 * t11 + 0.2e1 * t4 * t7, 0.2e1 * t4 * t1 - 0.2e1 * t6 * t2 + 0.2e1 * t5 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, t8, 0, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t3, t2, -t3, pkin(3) * t7 - t8 * qJ(4) - t10 * qJD(4), -t2, -t3 * pkin(3) - t2 * qJ(4) + t6 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, qJ(4) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
