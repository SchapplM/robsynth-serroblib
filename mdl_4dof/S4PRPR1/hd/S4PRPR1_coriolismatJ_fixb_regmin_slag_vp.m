% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPR1
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
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:03
% EndTime: 2019-03-08 18:21:04
% DurationCPUTime: 0.09s
% Computational Cost: add. (30->15), mult. (43->18), div. (0->0), fcn. (32->2), ass. (0->14)
t14 = -pkin(2) - pkin(3);
t6 = sin(qJ(4));
t13 = t6 * qJD(2);
t12 = t6 * qJD(3);
t7 = cos(qJ(4));
t11 = t7 * qJD(2);
t10 = t7 * qJD(3);
t9 = qJ(3) * qJD(2);
t8 = qJD(2) - qJD(4);
t4 = t8 * t7;
t3 = t8 * t6;
t2 = t7 * qJ(3) + t6 * t14;
t1 = t6 * qJ(3) - t7 * t14;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, t2 * qJD(4) + t12, -t1 * qJD(4) + t10; 0, 0, 0, 0, 0, qJD(2), t9, 0, t13, t11; 0, 0, 0, 0, 0, 0, 0, 0, t8 * t2, -t8 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t9, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(2) - t12, t1 * qJD(2) - t10; 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
