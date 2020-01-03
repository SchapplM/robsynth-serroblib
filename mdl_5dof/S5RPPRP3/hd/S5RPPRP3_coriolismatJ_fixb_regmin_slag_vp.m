% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:04
% EndTime: 2019-12-31 17:51:05
% DurationCPUTime: 0.22s
% Computational Cost: add. (249->43), mult. (367->49), div. (0->0), fcn. (266->4), ass. (0->35)
t19 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t48 = -qJ(5) + t19;
t28 = cos(qJ(4));
t25 = t28 ^ 2;
t47 = pkin(4) * t28;
t20 = sin(pkin(7)) * pkin(1) + qJ(3);
t27 = sin(qJ(4));
t14 = t27 * pkin(4) + t20;
t1 = t14 * t47;
t46 = t1 * qJD(1);
t12 = t48 * t27;
t6 = t12 * t27 + t48 * t25;
t45 = t6 * qJD(1);
t44 = qJD(4) * t27;
t43 = qJD(4) * t28;
t42 = t14 * qJD(1);
t24 = t27 ^ 2;
t29 = -t24 / 0.2e1 - t25 / 0.2e1;
t16 = -0.1e1 / 0.2e1 + t29;
t41 = t16 * qJD(1);
t17 = t24 - t25;
t40 = t17 * qJD(1);
t18 = -t24 - t25;
t39 = t18 * qJD(1);
t38 = t20 * qJD(1);
t37 = t27 * qJD(1);
t36 = t28 * qJD(1);
t35 = pkin(4) * t44;
t34 = pkin(4) * t43;
t33 = pkin(4) * t36;
t32 = t20 * t37;
t31 = t20 * t36;
t30 = t27 * t36;
t15 = 0.1e1 / 0.2e1 + t29;
t2 = [0, 0, 0, 0, 0, qJD(3), t20 * qJD(3), -t27 * t43, t17 * qJD(4), 0, 0, 0, qJD(3) * t27 + t20 * t43, qJD(3) * t28 - t20 * t44, -t18 * qJD(5), t14 * qJD(3) + t1 * qJD(4) - t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t38, 0, 0, 0, 0, 0, t37, t36, 0, t15 * qJD(5) + t42; 0, 0, 0, 0, 0, 0, 0, -t30, t40, -t44, -t43, 0, -t19 * t44 + t31, -t19 * t43 - t32, t35, -t12 * pkin(4) * qJD(4) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t15 * qJD(3) - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t44, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(1), -t38, 0, 0, 0, 0, 0, -t37, -t36, 0, t16 * qJD(5) - t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43, 0, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, t30, -t40, 0, 0, 0, -t31, t32, 0, -qJD(5) * t47 - t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t16 * qJD(3) + t34 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
