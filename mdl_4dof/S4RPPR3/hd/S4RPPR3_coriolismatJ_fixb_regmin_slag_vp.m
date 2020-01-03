% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:58
% EndTime: 2019-12-31 16:37:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (75->26), mult. (165->40), div. (0->0), fcn. (178->6), ass. (0->26)
t12 = sin(pkin(6)) * pkin(1) + qJ(3);
t29 = pkin(5) + t12;
t28 = cos(qJ(4));
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t11 = t16 ^ 2 + t17 ^ 2;
t18 = sin(qJ(4));
t7 = t18 * t16 - t28 * t17;
t9 = t28 * t16 + t18 * t17;
t1 = t7 ^ 2 - t9 ^ 2;
t27 = t1 * qJD(1);
t2 = t11 * t12;
t26 = t2 * qJD(1);
t25 = t7 * qJD(1);
t5 = t7 * qJD(4);
t24 = t9 * qJD(1);
t6 = t9 * qJD(4);
t10 = -cos(pkin(6)) * pkin(1) - pkin(2) - t17 * pkin(3);
t23 = qJD(1) * t10;
t22 = qJD(4) * t10;
t21 = t11 * qJD(1);
t20 = t7 * t24;
t19 = qJD(3) + t23;
t4 = t29 * t17;
t3 = t29 * t16;
t8 = [0, 0, 0, 0, 0, 0, t11 * qJD(3), t2 * qJD(3), -t7 * t6, t1 * qJD(4), 0, 0, 0, t9 * t22, -t7 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t21, t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t20, t27, -t5, -t6, 0, t9 * t23 + (t18 * t3 - t28 * t4) * qJD(4), -t7 * t23 + (t18 * t4 + t28 * t3) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, -t21, -t26, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25; 0, 0, 0, 0, 0, 0, 0, 0, t20, -t27, 0, 0, 0, -t19 * t9, t19 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
