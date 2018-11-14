% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x8]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:22
% EndTime: 2018-11-14 14:02:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (103->15), mult. (232->26), div. (0->0), fcn. (272->6), ass. (0->20)
t34 = qJD(2) + qJD(4);
t19 = sin(pkin(6));
t20 = cos(pkin(6));
t29 = sin(qJ(2));
t30 = cos(qJ(2));
t14 = -t19 * t29 + t20 * t30;
t15 = -t19 * t30 - t20 * t29;
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t33 = t34 * (-t21 * t14 + t22 * t15);
t32 = t34 * (-t22 * t14 - t21 * t15);
t31 = pkin(2) * t19;
t18 = t20 * pkin(2) + pkin(3);
t12 = -t22 * t18 + t21 * t31;
t24 = t12 * qJD(2);
t13 = t21 * t18 + t22 * t31;
t23 = t13 * qJD(2);
t7 = t13 * qJD(4);
t6 = t12 * qJD(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t29 * qJD(2), -t30 * qJD(2) (t14 * t19 + t15 * t20) * qJD(2) * pkin(2), 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t7 - t23, t6 + t24; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
