% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:44
% EndTime: 2019-12-31 16:30:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (112->47), mult. (290->73), div. (0->0), fcn. (221->4), ass. (0->47)
t21 = sin(qJ(3));
t50 = t21 * pkin(3);
t19 = t21 ^ 2;
t23 = cos(qJ(3));
t20 = t23 ^ 2;
t49 = t19 + t20;
t48 = t23 * qJ(4);
t28 = -t23 * pkin(3) - t21 * qJ(4);
t12 = -pkin(2) + t28;
t13 = -t48 + t50;
t3 = t12 * t23 + t13 * t21;
t47 = t3 * qJD(2);
t4 = -t12 * t21 + t13 * t23;
t46 = t4 * qJD(2);
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t5 = (-0.1e1 + t49) * t24 * t22;
t45 = t5 * qJD(1);
t44 = qJD(2) * t21;
t43 = qJD(2) * t23;
t16 = t20 - t19;
t42 = t16 * qJD(2);
t41 = t19 * qJD(2);
t40 = t21 * qJD(3);
t39 = t21 * qJD(4);
t38 = t22 * qJD(2);
t18 = t23 * qJD(3);
t37 = t23 * qJD(4);
t36 = t24 * qJD(2);
t35 = qJD(3) * qJ(4);
t34 = pkin(2) * t44;
t33 = pkin(2) * t43;
t32 = pkin(5) * t40;
t31 = pkin(5) * t18;
t30 = t12 * t44;
t29 = t49 * qJD(2);
t27 = -t50 / 0.2e1 + t48 / 0.2e1;
t1 = (t13 / 0.2e1 + t27) * t24;
t26 = -t12 * t13 * qJD(2) + t1 * qJD(1);
t25 = t28 * qJD(3) + t37;
t17 = t21 * t43;
t9 = -t23 * t38 - t24 * t40;
t8 = -t24 * t18 + t21 * t38;
t7 = t22 * t18 + t21 * t36;
t6 = t22 * t40 - t23 * t36;
t2 = (-t13 / 0.2e1 + t27) * t24;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2); 0, 0, -t38, -t36, 0, 0, 0, 0, 0, t9, t8, t9, t24 * t29, -t8, t12 * t38 + t45 + t2 * qJD(3) + (pkin(5) * t29 + t39) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t7, 0, -t6, t2 * qJD(2) + t25 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * qJD(3) - t45; 0, 0, 0, 0, t21 * t18, t16 * qJD(3), 0, 0, 0, -pkin(2) * t40, -pkin(2) * t18, -t4 * qJD(3) + t21 * t37, 0, -t3 * qJD(3) + t19 * qJD(4), (qJD(3) * t13 - t39) * t12; 0, 0, 0, 0, t17, t42, t18, -t40, 0, -t31 - t34, t32 - t33, -t31 - t46, t25, -t32 - t47, t25 * pkin(5) - t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, t41, -t30 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2); 0, 0, 0, 0, -t17, -t42, 0, 0, 0, t34, t33, t46, 0, t47, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t41, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
