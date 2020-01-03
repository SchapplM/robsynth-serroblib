% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:45
% EndTime: 2019-12-31 16:40:45
% DurationCPUTime: 0.18s
% Computational Cost: add. (100->39), mult. (233->55), div. (0->0), fcn. (232->4), ass. (0->34)
t40 = -pkin(5) + qJ(2);
t21 = sin(pkin(6));
t19 = t21 ^ 2;
t22 = cos(pkin(6));
t17 = t22 ^ 2 + t19;
t23 = sin(qJ(4));
t24 = cos(qJ(4));
t4 = t21 * t23 + t22 * t24;
t6 = t21 * t24 - t22 * t23;
t1 = t4 ^ 2 - t6 ^ 2;
t39 = t1 * qJD(1);
t38 = t4 * qJD(1);
t37 = t4 * qJD(4);
t36 = t6 * qJD(1);
t3 = t6 * qJD(4);
t11 = t17 * qJ(2);
t35 = t11 * qJD(1);
t34 = t17 * qJD(1);
t33 = t19 * qJD(1);
t32 = t21 * qJD(1);
t31 = t21 * qJD(3);
t30 = t4 * t36;
t29 = t4 * t32;
t28 = t6 * t32;
t27 = t22 * t32;
t26 = t21 * qJ(3) + pkin(1);
t2 = (pkin(2) + pkin(3)) * t22 + t26;
t25 = qJD(1) * t2 - qJD(2);
t14 = t17 * qJD(2);
t13 = t40 * t22;
t12 = t40 * t21;
t10 = -t22 * pkin(2) - t26;
t7 = t11 * qJD(2);
t5 = [0, 0, 0, 0, 0, t14, t7, t22 * t31, t14, t19 * qJD(3), -t10 * t31 + t7, -t4 * t3, t1 * qJD(4), 0, 0, 0, t2 * t3 + t4 * t31, -t2 * t37 + t6 * t31; 0, 0, 0, 0, 0, t34, t35, 0, t34, 0, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t27, 0, t33, -t10 * t32, 0, 0, 0, 0, 0, t29, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t39, -t37, -t3, 0, t2 * t36 + (-t23 * t12 - t24 * t13) * qJD(4), -t2 * t38 + (-t24 * t12 + t23 * t13) * qJD(4); 0, 0, 0, 0, 0, -t34, -t35, 0, -t34, 0, -t35 - t31, 0, 0, 0, 0, 0, -t3, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t38; 0, 0, 0, 0, 0, 0, 0, -t27, 0, -t33, (qJD(1) * t10 + qJD(2)) * t21, 0, 0, 0, 0, 0, -t29, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t23, -qJD(4) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t39, 0, 0, 0, -t25 * t6, t25 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
