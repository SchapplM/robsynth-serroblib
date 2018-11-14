% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:41
% EndTime: 2018-11-14 13:51:41
% DurationCPUTime: 0.14s
% Computational Cost: add. (124->42), mult. (275->56), div. (0->0), fcn. (205->4), ass. (0->34)
t17 = sin(pkin(6));
t37 = t17 * pkin(2);
t20 = cos(qJ(2));
t36 = t20 * pkin(1);
t15 = pkin(2) + t36;
t35 = t17 * t15;
t34 = t17 * t20;
t18 = cos(pkin(6));
t19 = sin(qJ(2));
t33 = t18 * t19;
t32 = pkin(1) * qJD(1);
t31 = pkin(1) * qJD(2);
t12 = t17 * t19 * pkin(1);
t24 = t18 * t15 - t12;
t21 = pkin(1) * t33 + t35;
t5 = qJ(4) + t21;
t8 = (t33 + t34) * pkin(1);
t9 = t18 * t36 - t12;
t1 = t5 * t9 + (-pkin(3) - t24) * t8;
t30 = t1 * qJD(1);
t2 = t21 * t9 - t24 * t8;
t29 = t2 * qJD(1);
t28 = t5 * qJD(1);
t27 = t8 * qJD(1);
t26 = t9 * qJD(1);
t25 = t9 * qJD(2) + qJD(4);
t16 = qJD(1) + qJD(2);
t23 = pkin(1) * t16;
t13 = qJ(4) + t37;
t4 = -qJ(4) + (t36 / 0.2e1 - pkin(2) / 0.2e1 - t15 / 0.2e1) * t17;
t22 = t4 * qJD(1) - t13 * qJD(2);
t6 = t8 * qJD(2);
t3 = t37 / 0.2e1 + qJ(4) + t35 / 0.2e1 + (t33 + t34 / 0.2e1) * pkin(1);
t7 = [0, 0, 0, 0, -t19 * t31, -t20 * t31, t2 * qJD(2), -t6, t25, t1 * qJD(2) + t5 * qJD(4); 0, 0, 0, 0, -t19 * t23, -t20 * t23, t29 + (t17 * t9 - t18 * t8) * qJD(2) * pkin(2), -t6 - t27, t25 + t26, t30 + (t9 * t13 + t8 * (-t18 * pkin(2) - pkin(3))) * qJD(2) + t3 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t16, t3 * qJD(2) + t28; 0, 0, 0, 0, t19 * t32, t20 * t32, -t29, t27, qJD(4) - t26, -t4 * qJD(4) - t30; 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t13 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t16, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t16, t4 * qJD(2) - t28; 0, 0, 0, 0, 0, 0, 0, 0, -t16, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
