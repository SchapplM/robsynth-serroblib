% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:17
% EndTime: 2019-12-31 17:45:18
% DurationCPUTime: 0.24s
% Computational Cost: add. (121->44), mult. (237->54), div. (0->0), fcn. (234->6), ass. (0->35)
t21 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t46 = -pkin(6) + t21;
t45 = sin(qJ(5));
t28 = sin(pkin(8));
t26 = t28 ^ 2;
t30 = cos(pkin(8));
t27 = t30 ^ 2;
t44 = t26 + t27;
t31 = cos(qJ(5));
t11 = t31 * t28 + t45 * t30;
t13 = -t45 * t28 + t31 * t30;
t1 = t11 ^ 2 - t13 ^ 2;
t43 = t1 * qJD(1);
t2 = t44 * t21;
t42 = t2 * qJD(1);
t41 = t11 * qJD(1);
t6 = t11 * qJD(5);
t40 = t13 * qJD(1);
t9 = t13 * qJD(5);
t33 = -t26 / 0.2e1 - t27 / 0.2e1;
t16 = -0.1e1 / 0.2e1 + t33;
t39 = t16 * qJD(1);
t38 = t44 * qJD(1);
t22 = sin(pkin(7)) * pkin(1) + qJ(3);
t37 = t22 * qJD(1);
t36 = t28 * qJD(1);
t35 = t30 * qJD(1);
t34 = t11 * t40;
t14 = t28 * pkin(4) + t22;
t32 = qJD(1) * t14 + qJD(4);
t17 = t22 * qJD(3);
t15 = 0.1e1 / 0.2e1 + t33;
t4 = t46 * t30;
t3 = t46 * t28;
t5 = [0, 0, 0, 0, 0, qJD(3), t17, qJD(3) * t28, qJD(3) * t30, t44 * qJD(4), -t2 * qJD(4) + t17, -t11 * t9, t1 * qJD(5), 0, 0, 0, qJD(3) * t11 + t14 * t9, qJD(3) * t13 - t14 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t37, t36, t35, 0, t15 * qJD(4) + t37, 0, 0, 0, 0, 0, t41, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t15 * qJD(3) - t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t43, -t6, -t9, 0, t14 * t40 + (-t31 * t3 - t45 * t4) * qJD(5), -t14 * t41 + (t45 * t3 - t31 * t4) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t6; 0, 0, 0, 0, 0, -qJD(1), -t37, -t36, -t35, 0, t16 * qJD(4) - t37, 0, 0, 0, 0, 0, -t41, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t16 * qJD(3) + t42, 0, 0, 0, 0, 0, t9, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t43, 0, 0, 0, -t32 * t13, t32 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
