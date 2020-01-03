% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:03
% EndTime: 2019-12-31 17:44:04
% DurationCPUTime: 0.23s
% Computational Cost: add. (127->42), mult. (264->57), div. (0->0), fcn. (265->6), ass. (0->35)
t19 = sin(pkin(7)) * pkin(1) + qJ(3);
t42 = -pkin(6) + t19;
t23 = sin(pkin(8));
t21 = t23 ^ 2;
t24 = cos(pkin(8));
t18 = t24 ^ 2 + t21;
t26 = sin(qJ(5));
t27 = cos(qJ(5));
t12 = t23 * t26 + t24 * t27;
t14 = t23 * t27 - t24 * t26;
t1 = t12 ^ 2 - t14 ^ 2;
t41 = t1 * qJD(1);
t5 = t18 * t19;
t40 = t5 * qJD(1);
t39 = t12 * qJD(1);
t10 = t12 * qJD(5);
t38 = t14 * qJD(1);
t11 = t14 * qJD(5);
t37 = t18 * qJD(1);
t36 = t21 * qJD(1);
t35 = t23 * qJD(1);
t34 = t23 * qJD(4);
t33 = t12 * t38;
t32 = t12 * t35;
t31 = t14 * t35;
t30 = t24 * t35;
t28 = cos(pkin(7)) * pkin(1) + t23 * qJ(4) + pkin(2);
t6 = (pkin(3) + pkin(4)) * t24 + t28;
t29 = qJD(1) * t6 - qJD(3);
t15 = t18 * qJD(3);
t9 = t42 * t24;
t8 = t42 * t23;
t7 = -t24 * pkin(3) - t28;
t2 = t5 * qJD(3);
t3 = [0, 0, 0, 0, 0, 0, t15, t2, t24 * t34, t15, t21 * qJD(4), -t7 * t34 + t2, -t12 * t11, t1 * qJD(5), 0, 0, 0, t6 * t11 + t12 * t34, -t6 * t10 + t14 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t37, t40, 0, t37, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t36, -t7 * t35, 0, 0, 0, 0, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t41, -t10, -t11, 0, t6 * t38 + (-t26 * t8 - t27 * t9) * qJD(5), -t6 * t39 + (t26 * t9 - t27 * t8) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, -t37, -t40, 0, -t37, 0, -t40 - t34, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t39; 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, -t36, (qJD(1) * t7 + qJD(3)) * t23, 0, 0, 0, 0, 0, -t32, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t26, -qJD(5) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t41, 0, 0, 0, -t29 * t14, t29 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
