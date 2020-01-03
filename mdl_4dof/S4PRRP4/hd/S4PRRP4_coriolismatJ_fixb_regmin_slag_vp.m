% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP4
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
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:59
% EndTime: 2019-12-31 16:27:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (69->29), mult. (147->41), div. (0->0), fcn. (106->2), ass. (0->29)
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t15 = -t12 * pkin(3) - t11 * qJ(4);
t4 = -pkin(2) + t15;
t5 = -t11 * pkin(3) + t12 * qJ(4);
t1 = -t5 * t11 + t4 * t12;
t30 = t1 * qJD(2);
t2 = -t4 * t11 - t5 * t12;
t29 = t2 * qJD(2);
t10 = t11 ^ 2;
t6 = t12 ^ 2 - t10;
t28 = t6 * qJD(2);
t27 = qJD(2) * t11;
t26 = qJD(2) * t12;
t25 = t10 * qJD(2);
t24 = t11 * qJD(3);
t8 = t12 * qJD(3);
t23 = t12 * qJD(4);
t22 = qJD(3) * qJ(4);
t21 = t4 * t5 * qJD(2);
t20 = pkin(2) * t27;
t19 = pkin(2) * t26;
t18 = pkin(5) * t24;
t17 = pkin(5) * t8;
t16 = t4 * t27;
t14 = t5 * qJD(3) + t11 * qJD(4);
t13 = t15 * qJD(3) + t23;
t7 = t11 * t26;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t8, -t24, 0, t8, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t11 * t8, t6 * qJD(3), 0, 0, 0, -pkin(2) * t24, -pkin(2) * t8, -t2 * qJD(3) + t11 * t23, 0, -t1 * qJD(3) + t10 * qJD(4), -t14 * t4; 0, 0, 0, 0, t7, t28, t8, -t24, 0, -t17 - t20, t18 - t19, -t17 - t29, t13, -t18 - t30, t13 * pkin(5) - t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t25, -t16 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t7, -t28, 0, 0, 0, t20, t19, t29, 0, t30, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t25, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
