% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP4
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
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:54
% EndTime: 2019-12-31 16:43:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (95->31), mult. (173->43), div. (0->0), fcn. (132->4), ass. (0->31)
t14 = sin(qJ(3));
t15 = cos(qJ(3));
t18 = -t15 * pkin(3) - t14 * qJ(4);
t9 = -cos(pkin(6)) * pkin(1) - pkin(2);
t4 = t18 + t9;
t5 = -t14 * pkin(3) + t15 * qJ(4);
t1 = -t5 * t14 + t4 * t15;
t33 = t1 * qJD(1);
t2 = -t4 * t14 - t5 * t15;
t32 = t2 * qJD(1);
t12 = t14 ^ 2;
t6 = t15 ^ 2 - t12;
t31 = t6 * qJD(1);
t30 = qJD(1) * t14;
t29 = qJD(1) * t15;
t28 = t12 * qJD(1);
t27 = t14 * qJD(3);
t10 = t15 * qJD(3);
t26 = t15 * qJD(4);
t25 = qJD(3) * qJ(4);
t24 = t4 * t5 * qJD(1);
t8 = sin(pkin(6)) * pkin(1) + pkin(5);
t23 = t8 * t27;
t22 = t8 * t10;
t21 = t4 * t30;
t20 = t9 * t30;
t19 = t9 * t29;
t17 = t5 * qJD(3) + t14 * qJD(4);
t16 = t18 * qJD(3) + t26;
t7 = t14 * t29;
t3 = [0, 0, 0, 0, t14 * t10, t6 * qJD(3), 0, 0, 0, t9 * t27, t9 * t10, -t2 * qJD(3) + t14 * t26, 0, -t1 * qJD(3) + t12 * qJD(4), -t17 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t7, t31, t10, -t27, 0, t20 - t22, t19 + t23, -t22 - t32, t16, -t23 - t33, t16 * t8 - t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t10, t28, -t21 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t10, -t27, 0, t10, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, -t7, -t31, 0, 0, 0, -t20, -t19, t32, 0, t33, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t28, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
