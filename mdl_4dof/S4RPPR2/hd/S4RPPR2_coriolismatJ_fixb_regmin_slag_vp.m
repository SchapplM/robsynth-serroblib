% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x12]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:33
% EndTime: 2019-03-08 18:28:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (93->26), mult. (134->34), div. (0->0), fcn. (136->4), ass. (0->25)
t14 = sin(pkin(6));
t15 = cos(pkin(6));
t16 = sin(qJ(4));
t31 = cos(qJ(4));
t10 = t31 * t14 + t16 * t15;
t26 = t10 * qJD(1);
t33 = -t10 * qJD(4) + t26;
t9 = t16 * t14 - t31 * t15;
t29 = t9 * qJD(1);
t32 = -t9 * qJD(4) + t29;
t17 = -pkin(1) - pkin(2);
t11 = t15 * qJ(2) + t14 * t17;
t19 = -t14 * qJ(2) + t15 * t17;
t3 = t11 * t15 - t19 * t14;
t30 = t3 * qJD(1);
t28 = t9 * qJD(2);
t25 = t10 * qJD(2);
t23 = t14 * qJD(1);
t22 = t15 * qJD(1);
t21 = qJ(2) * qJD(1);
t20 = qJD(1) - qJD(4);
t18 = -pkin(3) + t19;
t2 = t31 * t11 + t16 * t18;
t1 = t16 * t11 - t31 * t18;
t4 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t14 * qJD(2), t15 * qJD(2), t3 * qJD(2), 0, t2 * qJD(4) + t25, -t1 * qJD(4) - t28; 0, 0, 0, 0, qJD(1), t21, t23, t22, t30, 0, t26, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t2, -t20 * t1; 0, 0, 0, 0, -qJD(1), -t21, -t23, -t22, -t30, 0, -t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(1) - t25, t1 * qJD(1) + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t4;
