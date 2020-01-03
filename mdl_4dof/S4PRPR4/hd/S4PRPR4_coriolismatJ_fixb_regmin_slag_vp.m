% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPR4
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:01
% EndTime: 2019-12-31 16:22:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (21->19), mult. (46->21), div. (0->0), fcn. (32->2), ass. (0->15)
t2 = sin(qJ(4));
t3 = cos(qJ(4));
t1 = t2 ^ 2 - t3 ^ 2;
t14 = t1 * qJD(2);
t13 = t2 * qJD(2);
t12 = t2 * qJD(4);
t11 = t3 * qJD(2);
t10 = t3 * qJD(4);
t9 = qJ(3) * qJD(2);
t8 = qJ(3) * qJD(4);
t7 = t2 * t11;
t6 = t2 * t9;
t5 = t3 * t9;
t4 = -pkin(2) - pkin(5);
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t2 * t10, t1 * qJD(4), 0, 0, 0, qJD(3) * t2 + t3 * t8, qJD(3) * t3 - t2 * t8; 0, 0, 0, 0, 0, qJD(2), t9, 0, 0, 0, 0, 0, t13, t11; 0, 0, 0, 0, 0, 0, 0, -t7, t14, -t12, -t10, 0, -t4 * t12 + t5, -t4 * t10 - t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t9, 0, 0, 0, 0, 0, -t13, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t7, -t14, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
