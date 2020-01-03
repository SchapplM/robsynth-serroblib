% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:55
% EndTime: 2019-12-31 16:38:55
% DurationCPUTime: 0.10s
% Computational Cost: add. (31->21), mult. (57->22), div. (0->0), fcn. (43->4), ass. (0->15)
t4 = sin(qJ(4));
t5 = cos(qJ(4));
t1 = t4 ^ 2 - t5 ^ 2;
t14 = t1 * qJD(1);
t3 = sin(pkin(6)) * pkin(1) + qJ(3);
t13 = t3 * qJD(1);
t12 = t4 * qJD(1);
t11 = t4 * qJD(4);
t10 = t5 * qJD(1);
t9 = t5 * qJD(4);
t8 = t3 * t12;
t7 = t3 * t10;
t6 = t4 * t10;
t2 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5);
t15 = [0, 0, 0, 0, 0, qJD(3), t3 * qJD(3), -t4 * t9, t1 * qJD(4), 0, 0, 0, qJD(3) * t4 + t3 * t9, qJD(3) * t5 - t3 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t13, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, -t6, t14, -t11, -t9, 0, -t2 * t11 + t7, -t2 * t9 - t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t11; 0, 0, 0, 0, 0, -qJD(1), -t13, 0, 0, 0, 0, 0, -t12, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t9; 0, 0, 0, 0, 0, 0, 0, t6, -t14, 0, 0, 0, -t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
