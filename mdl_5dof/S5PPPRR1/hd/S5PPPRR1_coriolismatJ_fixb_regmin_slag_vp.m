% Calculate minimal parameter regressor of coriolis matrix for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x13]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPPRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:09
% EndTime: 2019-12-05 14:58:10
% DurationCPUTime: 0.15s
% Computational Cost: add. (113->23), mult. (357->49), div. (0->0), fcn. (384->8), ass. (0->27)
t29 = cos(qJ(5));
t39 = qJD(4) * t29;
t23 = sin(pkin(9));
t25 = cos(pkin(9));
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t19 = t28 * t23 - t30 * t25;
t24 = sin(pkin(8));
t14 = t19 * t24;
t38 = t14 * qJD(4);
t20 = t30 * t23 + t28 * t25;
t37 = t20 * qJD(4);
t27 = sin(qJ(5));
t21 = -t27 ^ 2 + t29 ^ 2;
t36 = t21 * qJD(4);
t35 = t27 * qJD(5);
t34 = t29 * qJD(5);
t33 = pkin(4) * t27 * qJD(4);
t32 = pkin(4) * t39;
t31 = t27 * t39;
t26 = cos(pkin(8));
t13 = t20 * t24;
t8 = t19 * t29;
t7 = t19 * t27;
t4 = t13 * t29;
t3 = t13 * t27;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t38, t13 * qJD(4), 0, 0, 0, 0, 0, t3 * qJD(5) + t29 * t38, t4 * qJD(5) - t27 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(4) + (t29 * t14 + t26 * t27) * qJD(5), t4 * qJD(4) + (-t27 * t14 + t26 * t29) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t37, t19 * qJD(4), 0, 0, 0, 0, 0, t7 * qJD(5) - t29 * t37, t8 * qJD(5) + t27 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(4) - t20 * t34, t8 * qJD(4) + t20 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t27 * t34, t21 * qJD(5), 0, 0, 0, -pkin(4) * t35, -pkin(4) * t34; 0, 0, 0, 0, 0, 0, t31, t36, t34, -t35, 0, -pkin(6) * t34 - t33, pkin(6) * t35 - t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t31, -t36, 0, 0, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
