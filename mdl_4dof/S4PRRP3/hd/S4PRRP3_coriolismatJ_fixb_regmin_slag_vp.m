% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x13]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (82->22), mult. (198->34), div. (0->0), fcn. (140->2), ass. (0->25)
t14 = sin(qJ(3));
t12 = t14 ^ 2;
t30 = pkin(3) * t14;
t29 = -qJ(4) - pkin(5);
t15 = cos(qJ(3));
t1 = (-pkin(3) * t15 - pkin(2)) * t30;
t28 = qJD(2) * t1;
t9 = t29 * t15;
t5 = -t29 * t12 - t9 * t15;
t27 = t5 * qJD(2);
t26 = qJD(2) * t14;
t25 = qJD(2) * t15;
t13 = t15 ^ 2;
t10 = t12 + t13;
t24 = t10 * qJD(2);
t11 = t13 - t12;
t23 = t11 * qJD(2);
t22 = t14 * qJD(3);
t21 = t15 * qJD(3);
t20 = pkin(2) * t26;
t19 = pkin(2) * t25;
t18 = pkin(3) * t22;
t17 = pkin(3) * t26;
t16 = t14 * t25;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t14 * t21, t11 * qJD(3), 0, 0, 0, -pkin(2) * t22, -pkin(2) * t21, qJD(4) * t10, qJD(3) * t1 + qJD(4) * t5; 0, 0, 0, 0, t16, t23, t21, -t22, 0, -pkin(5) * t21 - t20, pkin(5) * t22 - t19, -pkin(3) * t21, pkin(3) * qJD(3) * t9 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t16, -t23, 0, 0, 0, t20, t19, 0, -qJD(4) * t30 - t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t18 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
