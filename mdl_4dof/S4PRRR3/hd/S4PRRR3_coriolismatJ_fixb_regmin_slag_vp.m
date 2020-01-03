% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:40
% EndTime: 2019-12-31 16:31:40
% DurationCPUTime: 0.19s
% Computational Cost: add. (94->38), mult. (254->55), div. (0->0), fcn. (162->4), ass. (0->38)
t24 = cos(qJ(3));
t42 = t24 * pkin(2);
t17 = -pkin(3) - t42;
t47 = pkin(3) / 0.2e1 - t17 / 0.2e1;
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t14 = -t21 ^ 2 + t23 ^ 2;
t36 = qJD(2) + qJD(3);
t46 = t36 * t14;
t33 = -t42 / 0.2e1;
t45 = t33 - t47;
t41 = pkin(2) * qJD(2);
t40 = pkin(2) * qJD(3);
t39 = pkin(3) * qJD(3);
t38 = qJD(2) * t17;
t37 = t21 * qJD(4);
t20 = t23 * qJD(4);
t22 = sin(qJ(3));
t35 = t22 * t40;
t34 = t22 * t41;
t32 = t21 * t38;
t31 = t23 * t38;
t30 = pkin(2) * t36;
t29 = t21 * t34;
t28 = t22 * t30;
t27 = t33 + t47;
t1 = t27 * t21;
t26 = t1 * qJD(2) + t21 * t39;
t2 = t27 * t23;
t25 = t2 * qJD(2) + t23 * t39;
t16 = t22 * pkin(2) + pkin(6);
t15 = t21 * t20;
t13 = t21 * t35;
t10 = t14 * qJD(4);
t7 = t36 * t23 * t21;
t4 = t45 * t23;
t3 = t45 * t21;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t35, -t24 * t40, t15, t10, 0, 0, 0, t17 * t37 - t23 * t35, t17 * t20 + t13; 0, 0, 0, 0, 0, -t28, -t24 * t30, t15, t10, 0, 0, 0, t3 * qJD(4) - t23 * t28, t4 * qJD(4) + t13 + t29; 0, 0, 0, 0, 0, 0, 0, t7, t46, t20, -t37, 0, t3 * qJD(3) - t16 * t20 + t32, t4 * qJD(3) + t16 * t37 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t34, t24 * t41, t15, t10, 0, 0, 0, -t1 * qJD(4) + t23 * t34, -t2 * qJD(4) - t29; 0, 0, 0, 0, 0, 0, 0, t15, t10, 0, 0, 0, -pkin(3) * t37, -pkin(3) * t20; 0, 0, 0, 0, 0, 0, 0, t7, t46, t20, -t37, 0, -pkin(6) * t20 - t26, pkin(6) * t37 - t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t7, -t46, 0, 0, 0, t1 * qJD(3) - t32, t2 * qJD(3) - t31; 0, 0, 0, 0, 0, 0, 0, -t7, -t46, 0, 0, 0, t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
