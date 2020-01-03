% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:15
% EndTime: 2019-12-31 16:46:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (150->38), mult. (243->48), div. (0->0), fcn. (161->2), ass. (0->34)
t22 = -pkin(1) - pkin(5);
t42 = -qJ(4) + t22;
t21 = cos(qJ(3));
t19 = t21 ^ 2;
t41 = pkin(3) * t21;
t20 = sin(qJ(3));
t14 = t20 * pkin(3) + qJ(2);
t1 = t14 * t41;
t40 = t1 * qJD(1);
t8 = t42 * t20;
t4 = t42 * t19 + t8 * t20;
t39 = t4 * qJD(1);
t18 = t20 ^ 2;
t25 = -t18 / 0.2e1 - t19 / 0.2e1;
t11 = -0.1e1 / 0.2e1 + t25;
t38 = t11 * qJD(1);
t12 = t18 + t19;
t37 = t12 * qJD(1);
t13 = t18 - t19;
t36 = t13 * qJD(1);
t35 = t14 * qJD(1);
t34 = t20 * qJD(1);
t33 = t20 * qJD(3);
t32 = t21 * qJD(1);
t31 = t21 * qJD(3);
t30 = qJ(2) * qJD(3);
t29 = qJD(1) * qJ(2);
t28 = pkin(3) * t33;
t27 = pkin(3) * t32;
t26 = t20 * t32;
t24 = t20 * t29;
t23 = t21 * t29;
t10 = 0.1e1 / 0.2e1 + t25;
t2 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t20 * t31, t13 * qJD(3), 0, 0, 0, qJD(2) * t20 + t21 * t30, qJD(2) * t21 - t20 * t30, t12 * qJD(4), t14 * qJD(2) + t1 * qJD(3) - t4 * qJD(4); 0, 0, 0, 0, qJD(1), t29, 0, 0, 0, 0, 0, t34, t32, 0, t10 * qJD(4) + t35; 0, 0, 0, 0, 0, 0, -t26, t36, -t33, -t31, 0, -t22 * t33 + t23, -t22 * t31 - t24, t28, -t8 * pkin(3) * qJD(3) + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t10 * qJD(2) - t39; 0, 0, 0, 0, -qJD(1), -t29, 0, 0, 0, 0, 0, -t34, -t32, 0, t11 * qJD(4) - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t31, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, t26, -t36, 0, 0, 0, -t23, t24, 0, -qJD(4) * t41 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, pkin(3) * t31 - t11 * qJD(2) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
