% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:01
% DurationCPUTime: 0.13s
% Computational Cost: add. (56->26), mult. (154->53), div. (0->0), fcn. (121->6), ass. (0->22)
t15 = sin(pkin(8));
t13 = t15 ^ 2;
t16 = cos(pkin(8));
t27 = (t16 ^ 2 + t13) * qJD(3);
t12 = sin(pkin(7)) * pkin(1) + qJ(3);
t26 = -pkin(6) + t12;
t25 = t12 * t27;
t18 = sin(qJ(5));
t23 = qJD(5) * t18;
t19 = cos(qJ(5));
t22 = qJD(5) * t19;
t21 = t15 * qJD(4);
t8 = t15 * t19 - t16 * t18;
t7 = t15 * t18 + t16 * t19;
t20 = cos(pkin(7)) * pkin(1) + t15 * qJ(4) + pkin(2);
t6 = 0.2e1 * t27;
t5 = t26 * t16;
t4 = t26 * t15;
t3 = t15 * t22 - t16 * t23;
t2 = t7 * qJD(5);
t1 = (pkin(3) + pkin(4)) * t16 + t20;
t9 = [0, 0, 0, 0, 0, 0, t6, 0.2e1 * t25, 0.2e1 * t16 * t21, t6, 0.2e1 * t13 * qJD(4), -0.2e1 * (-t16 * pkin(3) - t20) * t21 + 0.2e1 * t25, -0.2e1 * t8 * t2, 0.2e1 * t2 * t7 - 0.2e1 * t8 * t3, 0, 0, 0, 0.2e1 * t1 * t3 + 0.2e1 * t7 * t21, -0.2e1 * t1 * t2 + 0.2e1 * t8 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, -t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, 0, (-t18 * t4 - t19 * t5) * qJD(5) + t8 * qJD(3), (t18 * t5 - t19 * t4) * qJD(5) - t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
