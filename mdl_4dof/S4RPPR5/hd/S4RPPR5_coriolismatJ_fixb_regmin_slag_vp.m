% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR5
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
% cmat_reg [(4*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:50
% EndTime: 2019-12-31 16:39:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (63->29), mult. (113->42), div. (0->0), fcn. (96->4), ass. (0->29)
t11 = -pkin(1) - pkin(2);
t7 = sin(pkin(6));
t8 = cos(pkin(6));
t30 = t8 * qJ(2) + t7 * t11;
t13 = -t7 * qJ(2) + t8 * t11;
t2 = pkin(3) - t13;
t29 = qJD(1) * t2;
t1 = -t13 * t7 + t30 * t8;
t28 = t1 * qJD(1);
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t4 = t10 ^ 2 - t9 ^ 2;
t27 = t4 * qJD(1);
t26 = t7 * qJD(1);
t25 = t7 * qJD(2);
t24 = t8 * qJD(1);
t23 = t8 * qJD(2);
t22 = t9 * qJD(4);
t21 = qJD(1) * t10;
t20 = t10 * qJD(4);
t19 = qJ(2) * qJD(1);
t18 = t9 * t26;
t17 = t9 * t24;
t16 = t7 * t21;
t15 = t8 * t21;
t14 = t9 * t21;
t12 = -t23 + t29;
t3 = -pkin(5) + t30;
t5 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t25, t23, t1 * qJD(2), t9 * t20, t4 * qJD(4), 0, 0, 0, t10 * t25 - t2 * t22, -t2 * t20 - t9 * t25; 0, 0, 0, 0, qJD(1), t19, t26, t24, t28, 0, 0, 0, 0, 0, t16, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t27, -t20, t22, 0, -t3 * t20 - t9 * t29, -t2 * t21 + t3 * t22; 0, 0, 0, 0, -qJD(1), -t19, -t26, -t24, -t28, 0, 0, 0, 0, 0, t8 * t22 - t16, t8 * t20 + t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t20 + t17, t7 * t22 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t27, 0, 0, 0, t12 * t9, t12 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
