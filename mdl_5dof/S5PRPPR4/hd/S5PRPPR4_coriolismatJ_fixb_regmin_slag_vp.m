% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:57
% EndTime: 2019-12-31 17:36:58
% DurationCPUTime: 0.21s
% Computational Cost: add. (102->40), mult. (239->55), div. (0->0), fcn. (240->4), ass. (0->34)
t40 = -pkin(6) + qJ(3);
t22 = sin(pkin(8));
t20 = t22 ^ 2;
t23 = cos(pkin(8));
t18 = t23 ^ 2 + t20;
t24 = sin(qJ(5));
t25 = cos(qJ(5));
t5 = t22 * t24 + t23 * t25;
t7 = t22 * t25 - t23 * t24;
t1 = t5 ^ 2 - t7 ^ 2;
t39 = t1 * qJD(2);
t38 = t5 * qJD(2);
t3 = t5 * qJD(5);
t37 = t7 * qJD(2);
t4 = t7 * qJD(5);
t12 = t18 * qJ(3);
t36 = t12 * qJD(2);
t35 = t18 * qJD(2);
t34 = t20 * qJD(2);
t33 = t22 * qJD(2);
t32 = t22 * qJD(4);
t31 = t5 * t37;
t30 = t5 * t33;
t29 = t7 * t33;
t28 = t23 * t33;
t27 = t22 * qJ(4) + pkin(2);
t2 = (pkin(3) + pkin(4)) * t23 + t27;
t26 = qJD(2) * t2 - qJD(3);
t15 = t18 * qJD(3);
t14 = t40 * t23;
t13 = t40 * t22;
t11 = -t23 * pkin(3) - t27;
t8 = t12 * qJD(3);
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t15, t8, t23 * t32, t15, t20 * qJD(4), -t11 * t32 + t8, -t5 * t4, t1 * qJD(5), 0, 0, 0, t2 * t4 + t5 * t32, -t2 * t3 + t7 * t32; 0, 0, 0, 0, 0, 0, t35, t36, 0, t35, 0, t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t34, -t11 * t33, 0, 0, 0, 0, 0, t30, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t39, -t3, -t4, 0, t2 * t37 + (-t24 * t13 - t25 * t14) * qJD(5), -t2 * t38 + (-t25 * t13 + t24 * t14) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t35, 0, -t36 - t32, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t34, (qJD(2) * t11 + qJD(3)) * t22, 0, 0, 0, 0, 0, -t30, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t24, -qJD(5) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t39, 0, 0, 0, -t26 * t7, t26 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
