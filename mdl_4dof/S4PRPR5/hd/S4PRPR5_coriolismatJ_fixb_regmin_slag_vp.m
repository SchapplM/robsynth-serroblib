% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x12]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:19
% EndTime: 2019-12-31 16:23:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (86->21), mult. (211->38), div. (0->0), fcn. (232->6), ass. (0->22)
t31 = cos(qJ(2));
t30 = sin(qJ(2));
t29 = cos(pkin(7));
t28 = sin(pkin(7));
t18 = sin(qJ(4));
t27 = qJD(2) * t18;
t19 = cos(qJ(4));
t26 = qJD(2) * t19;
t12 = -t18 ^ 2 + t19 ^ 2;
t25 = t12 * qJD(2);
t24 = t18 * qJD(4);
t23 = t19 * qJD(4);
t17 = -pkin(2) * t29 - pkin(3);
t22 = t17 * t27;
t21 = t17 * t26;
t20 = t18 * t26;
t16 = pkin(2) * t28 + pkin(5);
t11 = -t28 * t31 - t29 * t30;
t10 = t28 * t30 - t29 * t31;
t5 = t10 * t19;
t4 = t10 * t18;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t30 * qJD(2), -t31 * qJD(2), (-t10 * t28 + t11 * t29) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, qJD(4) * t4 + t11 * t26, qJD(4) * t5 - t11 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t4 + t11 * t23, qJD(2) * t5 - t11 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t18 * t23, t12 * qJD(4), 0, 0, 0, t17 * t24, t17 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t20, t25, t23, -t24, 0, -t16 * t23 + t22, t16 * t24 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t20, -t25, 0, 0, 0, -t22, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
