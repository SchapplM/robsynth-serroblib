% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (54->26), mult. (146->49), div. (0->0), fcn. (131->6), ass. (0->20)
t20 = 2 * qJD(4);
t11 = sin(qJ(2));
t19 = qJD(2) * t11;
t13 = cos(qJ(2));
t18 = qJD(2) * t13;
t10 = sin(qJ(5));
t17 = qJD(5) * t10;
t12 = cos(qJ(5));
t16 = qJD(5) * t12;
t9 = cos(pkin(8));
t15 = -t9 * pkin(2) - pkin(3);
t8 = sin(pkin(8));
t4 = t9 * t11 + t8 * t13;
t1 = t4 * qJD(2);
t2 = t9 * t18 - t8 * t19;
t3 = t8 * t11 - t9 * t13;
t14 = 0.2e1 * t3 * t1 + 0.2e1 * t4 * t2;
t7 = t8 * pkin(2) + qJ(4);
t6 = -pkin(6) + t15;
t5 = [0, 0, 0, 0, t14, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t19, -t18, (-t1 * t9 + t2 * t8) * pkin(2), t1, t2, t4 * qJD(4) + t1 * t15 + t2 * t7, 0, 0, 0, 0, 0, t2 * t10 + t4 * t16, t2 * t12 - t4 * t17; 0, 0, 0, 0, 0, 0, t20, t7 * t20, -0.2e1 * t10 * t16, 0.2e1 * (t10 ^ 2 - t12 ^ 2) * qJD(5), 0, 0, 0, 0.2e1 * qJD(4) * t10 + 0.2e1 * t7 * t16, 0.2e1 * qJD(4) * t12 - 0.2e1 * t7 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t1 - t3 * t17, -t10 * t1 - t3 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t16, 0, -t6 * t17, -t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
