% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:56
% EndTime: 2019-12-31 16:25:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (28->26), mult. (69->35), div. (0->0), fcn. (54->4), ass. (0->19)
t3 = sin(qJ(4));
t5 = cos(qJ(4));
t1 = t3 ^ 2 - t5 ^ 2;
t18 = t1 * qJD(2);
t17 = t3 * qJD(2);
t16 = t3 * qJD(4);
t4 = sin(qJ(2));
t2 = t4 * qJD(2);
t15 = t5 * qJD(2);
t14 = t5 * qJD(4);
t6 = cos(qJ(2));
t13 = t6 * qJD(2);
t12 = qJ(3) * qJD(4);
t11 = qJD(2) * qJ(3);
t10 = t3 * t15;
t9 = t3 * t11;
t8 = t5 * t11;
t7 = -pkin(2) - pkin(5);
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t2, -t13, t2, t13, (-t4 * pkin(2) + t6 * qJ(3)) * qJD(2) + t4 * qJD(3), 0, 0, 0, 0, 0, t3 * t13 + t4 * t14, t5 * t13 - t4 * t16; 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t15 + t6 * t16, t6 * t14 - t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t3 * t14, t1 * qJD(4), 0, 0, 0, qJD(3) * t3 + t5 * t12, qJD(3) * t5 - t3 * t12; 0, 0, 0, 0, 0, qJD(2), t11, 0, 0, 0, 0, 0, t17, t15; 0, 0, 0, 0, 0, 0, 0, -t10, t18, -t16, -t14, 0, -t7 * t16 + t8, -t7 * t14 - t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t11, 0, 0, 0, 0, 0, -t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t10, -t18, 0, 0, 0, -t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t19;
