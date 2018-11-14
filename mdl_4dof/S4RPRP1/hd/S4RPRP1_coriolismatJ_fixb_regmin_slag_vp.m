% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:39
% EndTime: 2018-11-14 13:48:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (104->24), mult. (180->23), div. (0->0), fcn. (140->4), ass. (0->20)
t25 = pkin(1) * sin(pkin(6));
t24 = cos(qJ(3));
t18 = sin(qJ(3));
t19 = cos(pkin(6)) * pkin(1) + pkin(2);
t12 = t18 * t19 + t24 * t25;
t10 = qJ(4) + t12;
t11 = t18 * t25 - t24 * t19;
t1 = -t10 * t11 + (-pkin(3) + t11) * t12;
t23 = t1 * qJD(1);
t20 = t11 * qJD(3);
t22 = -t20 + qJD(4);
t21 = t11 * qJD(1);
t8 = t12 * qJD(1);
t15 = qJD(1) + qJD(3);
t16 = qJ(4) * qJD(4);
t14 = t15 * qJ(4);
t9 = t12 * qJD(3);
t3 = t10 * qJD(4);
t2 = -t9 - t8;
t4 = [0, 0, 0, 0, 0, -t9, t20, -t9, t22, t1 * qJD(3) + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t2, t21 + t20, t2, -t21 + t22, t23 + (-t12 * pkin(3) - t11 * qJ(4)) * qJD(3) + t3; 0, 0, 0, 0, 0, 0, 0, 0, t15, t15 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t8, -t21, t8, t21 + qJD(4), t16 - t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t16; 0, 0, 0, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, -t15, -qJ(4) * qJD(3) - t10 * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t4;
