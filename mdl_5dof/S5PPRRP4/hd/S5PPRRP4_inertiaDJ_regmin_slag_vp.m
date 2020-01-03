% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:41
% DurationCPUTime: 0.14s
% Computational Cost: add. (59->27), mult. (186->60), div. (0->0), fcn. (123->4), ass. (0->24)
t10 = sin(qJ(4));
t8 = t10 ^ 2;
t12 = cos(qJ(4));
t9 = t12 ^ 2;
t23 = t8 + t9;
t22 = -qJ(5) - pkin(6);
t21 = t10 * qJD(4);
t11 = sin(qJ(3));
t20 = t11 * qJD(3);
t7 = t12 * qJD(4);
t13 = cos(qJ(3));
t19 = t13 * qJD(3);
t18 = -0.2e1 * pkin(3) * qJD(4);
t6 = pkin(4) * t21;
t17 = qJD(4) * t22;
t3 = t22 * t10;
t4 = t22 * t12;
t16 = t10 * t3 + t12 * t4;
t15 = -t10 * t19 - t11 * t7;
t1 = t12 * qJD(5) + t10 * t17;
t2 = -t10 * qJD(5) + t12 * t17;
t14 = t1 * t12 - t2 * t10 + (t10 * t4 - t12 * t3) * qJD(4);
t5 = -t12 * pkin(4) - pkin(3);
t24 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t23) * t11 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * qJD(4) - t10 * t1 - t12 * t2; 0, 0, 0, -t20, -t19, 0, 0, 0, 0, 0, -t12 * t20 - t13 * t21, t10 * t20 - t13 * t7, t23 * t19, (-t16 * qJD(3) - t6) * t13 + (qJD(3) * t5 + t14) * t11; 0, 0, 0, 0, 0, 0.2e1 * t10 * t7, 0.2e1 * (-t8 + t9) * qJD(4), 0, 0, 0, t10 * t18, t12 * t18, 0.2e1 * t14, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t5 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t7, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t11 * t21 - t12 * t19, 0, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, t7, -t21, 0, -pkin(6) * t7, pkin(6) * t21, -pkin(4) * t7, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t24;
