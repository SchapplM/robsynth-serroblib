% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:23
% EndTime: 2019-12-05 15:03:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (59->29), mult. (136->40), div. (0->0), fcn. (152->6), ass. (0->23)
t25 = cos(qJ(3));
t8 = sin(qJ(5));
t24 = qJD(5) * t8;
t20 = cos(pkin(8));
t7 = sin(pkin(8));
t9 = sin(qJ(3));
t3 = -t25 * t20 + t9 * t7;
t23 = t3 * qJD(3);
t4 = t9 * t20 + t25 * t7;
t2 = t4 * qJD(3);
t10 = cos(qJ(5));
t5 = -t10 ^ 2 + t8 ^ 2;
t22 = t5 * qJD(3);
t21 = t8 * qJD(3);
t19 = qJD(5) * t10;
t18 = qJD(5) * (-pkin(3) - pkin(6));
t17 = t10 * qJD(3);
t16 = qJ(4) * qJD(5);
t15 = qJD(3) * qJ(4);
t14 = t8 * t17;
t13 = t8 * t15;
t12 = t10 * t15;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t2, t23, t2, -t23, (-t4 * pkin(3) - t3 * qJ(4)) * qJD(3) + t4 * qJD(4), 0, 0, 0, 0, 0, t4 * t19 - t3 * t21, -t3 * t17 - t4 * t24; 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t17 - t3 * t24, -t3 * t19 - t4 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t8 * t19, t5 * qJD(5), 0, 0, 0, qJD(4) * t8 + t10 * t16, qJD(4) * t10 - t8 * t16; 0, 0, 0, 0, 0, 0, qJD(3), t15, 0, 0, 0, 0, 0, t21, t17; 0, 0, 0, 0, 0, 0, 0, 0, -t14, t22, -t24, -t19, 0, -t8 * t18 + t12, -t10 * t18 - t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -qJD(3), -t15, 0, 0, 0, 0, 0, -t21, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t14, -t22, 0, 0, 0, -t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
