% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP3
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
% cmat_reg [(4*%NQJ)%x13]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:49
% EndTime: 2019-12-31 16:42:49
% DurationCPUTime: 0.12s
% Computational Cost: add. (130->24), mult. (246->36), div. (0->0), fcn. (188->4), ass. (0->27)
t17 = sin(qJ(3));
t14 = t17 ^ 2;
t33 = pkin(3) * t17;
t13 = -cos(pkin(6)) * pkin(1) - pkin(2);
t18 = cos(qJ(3));
t1 = (-pkin(3) * t18 + t13) * t33;
t32 = qJD(1) * t1;
t12 = sin(pkin(6)) * pkin(1) + pkin(5);
t30 = qJ(4) + t12;
t9 = t30 * t18;
t5 = t30 * t14 + t9 * t18;
t31 = t5 * qJD(1);
t29 = qJD(1) * t17;
t28 = qJD(1) * t18;
t15 = t18 ^ 2;
t10 = t14 + t15;
t27 = t10 * qJD(1);
t11 = t15 - t14;
t26 = t11 * qJD(1);
t25 = t17 * qJD(3);
t24 = t18 * qJD(3);
t23 = pkin(3) * t25;
t22 = pkin(3) * t29;
t21 = t13 * t29;
t20 = t13 * t28;
t19 = t17 * t28;
t2 = [0, 0, 0, 0, t17 * t24, t11 * qJD(3), 0, 0, 0, t13 * t25, t13 * t24, qJD(4) * t10, qJD(3) * t1 + qJD(4) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t19, t26, t24, -t25, 0, -t12 * t24 + t21, t12 * t25 + t20, -pkin(3) * t24, -pkin(3) * qJD(3) * t9 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t19, -t26, 0, 0, 0, -t21, -t20, 0, -qJD(4) * t33 - t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t23 - t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
