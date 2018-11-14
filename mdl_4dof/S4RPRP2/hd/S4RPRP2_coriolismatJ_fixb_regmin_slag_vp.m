% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP2
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
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:35
% EndTime: 2018-11-14 13:49:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (165->29), mult. (217->41), div. (0->0), fcn. (155->2), ass. (0->25)
t30 = -pkin(1) - pkin(2);
t16 = sin(qJ(3));
t29 = t16 * pkin(3);
t17 = cos(qJ(3));
t10 = t16 * qJ(2) - t17 * t30;
t9 = -pkin(3) - t10;
t28 = t9 * t16;
t1 = (-pkin(3) / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t16;
t27 = t1 * qJD(1);
t11 = t17 * qJ(2) + t16 * t30;
t3 = (-t10 - t9) * t11;
t26 = t3 * qJD(1);
t4 = t11 * t17 - t28;
t25 = t4 * qJD(1);
t24 = t11 * qJD(3);
t23 = t16 * qJD(1);
t22 = t16 * qJD(2);
t21 = t17 * qJD(1);
t20 = t17 * qJD(2);
t19 = qJ(2) * qJD(1);
t18 = qJD(1) - qJD(3);
t13 = t18 * t17;
t12 = t18 * t16;
t2 = -t10 * t16 / 0.2e1 - t28 / 0.2e1 - t29 / 0.2e1;
t5 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, t22 + t24, -t10 * qJD(3) + t20, t4 * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, qJD(1), t19, 0, t23, t21, t2 * qJD(3) + t25; 0, 0, 0, 0, 0, 0, 0, t18 * t11, -t18 * t10, -pkin(3) * t24 + t2 * qJD(2) + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -qJD(1), -t19, 0, -t12, -t13, -t1 * qJD(3) - t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t12, t13, -qJD(3) * t29 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t11 * qJD(1) - t22, t10 * qJD(1) - t20, t1 * qJD(2) - t26; 0, 0, 0, 0, 0, 0, 0, -t23, -t21, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
