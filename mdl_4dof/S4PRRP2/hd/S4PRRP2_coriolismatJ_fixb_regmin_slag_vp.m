% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x8]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (105->18), mult. (235->39), div. (0->0), fcn. (249->4), ass. (0->23)
t16 = sin(qJ(3));
t29 = pkin(2) * t16;
t17 = sin(qJ(2));
t18 = cos(qJ(3));
t19 = cos(qJ(2));
t13 = t16 * t19 + t18 * t17;
t28 = t13 * pkin(3);
t27 = t18 * pkin(2);
t26 = pkin(2) * qJD(3);
t25 = qJD(2) * pkin(2);
t24 = -qJD(2) - qJD(3);
t23 = t16 * t26;
t22 = pkin(2) * t24;
t15 = pkin(3) + t27;
t21 = t27 / 0.2e1 - t15 / 0.2e1;
t10 = (-t15 + t27) * t29;
t4 = (pkin(3) / 0.2e1 + t21) * t13;
t20 = -t4 * qJD(1) - t10 * qJD(2);
t12 = -t16 * t17 + t18 * t19;
t7 = t24 * t12;
t6 = t24 * t13;
t3 = -t28 / 0.2e1 + t21 * t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t17 * qJD(2), -t19 * qJD(2), 0, t6, t7 (t12 * t29 - t13 * t15) * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, 0, t6, t7, t3 * qJD(2) - qJD(3) * t28; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t4 * qJD(3); 0, 0, 0, 0, 0, -t23, -t18 * t26, t10 * qJD(3); 0, 0, 0, 0, 0, t16 * t22, t18 * t22, -pkin(3) * t23 - t20; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t4 * qJD(2); 0, 0, 0, 0, 0, t16 * t25, t18 * t25, t20; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
