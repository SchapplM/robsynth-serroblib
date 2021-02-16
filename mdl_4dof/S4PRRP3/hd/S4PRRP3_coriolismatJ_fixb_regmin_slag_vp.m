% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP3
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
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:06
% EndTime: 2021-01-14 22:27:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (104->33), mult. (244->45), div. (0->0), fcn. (175->2), ass. (0->34)
t20 = cos(qJ(3));
t39 = t20 * pkin(3);
t15 = -pkin(2) - t39;
t19 = sin(qJ(3));
t38 = t15 * t19;
t37 = pkin(5) + qJ(4);
t1 = pkin(3) * t38;
t36 = t1 * qJD(2);
t10 = t37 * t19;
t11 = t37 * t20;
t5 = t10 * t19 + t11 * t20;
t35 = t5 * qJD(2);
t6 = t19 * t39 - t38;
t34 = t6 * qJD(2);
t17 = t19 ^ 2;
t9 = t17 * pkin(3) + t15 * t20;
t33 = t9 * qJD(2);
t32 = qJD(3) * t19;
t16 = qJD(3) * t20;
t31 = t11 * qJD(3);
t18 = t20 ^ 2;
t13 = t17 + t18;
t30 = t13 * qJD(2);
t14 = t18 - t17;
t29 = t14 * qJD(2);
t28 = t19 * qJD(2);
t27 = t19 * qJD(4);
t26 = t20 * qJD(2);
t25 = pkin(2) * t28;
t24 = pkin(2) * t26;
t23 = pkin(3) * t32;
t22 = pkin(3) * t28;
t21 = t19 * t26;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t16, -t32, -t16, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t19 * t16, t14 * qJD(3), 0, 0, 0, -pkin(2) * t32, -pkin(2) * t16, -t6 * qJD(3), t9 * qJD(3), t13 * qJD(4), t1 * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, t21, t29, t16, -t32, 0, -pkin(5) * t16 - t25, pkin(5) * t32 - t24, -t31 - t34, t10 * qJD(3) + t33, -pkin(3) * t16, -pkin(3) * t31 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t21, -t29, 0, 0, 0, t25, t24, -t27 + t34, -t20 * qJD(4) - t33, 0, -pkin(3) * t27 - t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t26, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t16, -t30, t23 - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t26, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
