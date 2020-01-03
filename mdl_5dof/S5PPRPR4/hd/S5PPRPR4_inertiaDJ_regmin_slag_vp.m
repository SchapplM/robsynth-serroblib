% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:24
% EndTime: 2019-12-31 17:32:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (51->25), mult. (184->53), div. (0->0), fcn. (154->6), ass. (0->21)
t10 = sin(pkin(8));
t11 = cos(pkin(8));
t12 = sin(qJ(5));
t14 = cos(qJ(5));
t3 = t12 * t10 - t14 * t11;
t1 = t3 * qJD(5);
t23 = -0.2e1 * t11 * pkin(4) - (2 * pkin(3));
t22 = t10 ^ 2 + t11 ^ 2;
t21 = pkin(6) + qJ(4);
t13 = sin(qJ(3));
t20 = qJD(3) * t13;
t15 = cos(qJ(3));
t19 = qJD(3) * t15;
t18 = t22 * t15;
t17 = t22 * qJD(4);
t16 = 0.2e1 * t17;
t4 = t14 * t10 + t12 * t11;
t2 = t4 * qJD(5);
t6 = t21 * t11;
t5 = t21 * t10;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t22) * t13 * t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t20, -t19, -t11 * t20, t10 * t20, qJD(3) * t18, t13 * t17 + (-pkin(3) * t13 + qJ(4) * t18) * qJD(3), 0, 0, 0, 0, 0, -t15 * t2 + t3 * t20, t15 * t1 + t4 * t20; 0, 0, 0, 0, 0, 0, 0, t16, qJ(4) * t16, -0.2e1 * t4 * t1, 0.2e1 * t1 * t3 - 0.2e1 * t4 * t2, 0, 0, 0, t2 * t23, -t1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t1 - t4 * t19, t13 * t2 + t3 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (t12 * t5 - t14 * t6) * qJD(5) - t4 * qJD(4), (t12 * t6 + t14 * t5) * qJD(5) + t3 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
