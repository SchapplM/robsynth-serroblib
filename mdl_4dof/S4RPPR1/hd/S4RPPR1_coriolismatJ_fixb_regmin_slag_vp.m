% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR1
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
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:32
% EndTime: 2019-03-08 18:27:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (48->17), mult. (62->20), div. (0->0), fcn. (51->4), ass. (0->15)
t6 = sin(pkin(6)) * pkin(1) + qJ(3);
t15 = t6 * qJD(1);
t7 = sin(qJ(4));
t14 = t7 * qJD(1);
t13 = t7 * qJD(3);
t8 = cos(qJ(4));
t12 = t8 * qJD(1);
t11 = t8 * qJD(3);
t10 = qJD(1) - qJD(4);
t9 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t5 = t10 * t8;
t4 = t10 * t7;
t2 = t8 * t6 + t7 * t9;
t1 = t7 * t6 - t8 * t9;
t3 = [0, 0, 0, 0, 0, qJD(3), t6 * qJD(3), 0, t2 * qJD(4) + t13, -t1 * qJD(4) + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t15, 0, t14, t12; 0, 0, 0, 0, 0, 0, 0, 0, t10 * t2, -t10 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(1), -t15, 0, -t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(1) - t13, t1 * qJD(1) - t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t3;
