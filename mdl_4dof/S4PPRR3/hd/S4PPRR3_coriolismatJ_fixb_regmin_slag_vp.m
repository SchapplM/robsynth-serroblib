% Calculate minimal parameter regressor of coriolis matrix for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x12]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PPRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:28
% EndTime: 2019-12-31 16:17:28
% DurationCPUTime: 0.10s
% Computational Cost: add. (15->14), mult. (53->25), div. (0->0), fcn. (42->4), ass. (0->15)
t5 = cos(qJ(4));
t14 = qJD(3) * t5;
t3 = sin(qJ(4));
t1 = -t3 ^ 2 + t5 ^ 2;
t13 = t1 * qJD(3);
t12 = t3 * qJD(4);
t4 = sin(qJ(3));
t11 = t4 * qJD(3);
t2 = t5 * qJD(4);
t6 = cos(qJ(3));
t10 = t6 * qJD(3);
t9 = pkin(3) * t3 * qJD(3);
t8 = pkin(3) * t14;
t7 = t3 * t14;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t11, -t10, 0, 0, 0, 0, 0, -t11 * t5 - t12 * t6, t11 * t3 - t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t3 - t2 * t4, -t10 * t5 + t12 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t3 * t2, t1 * qJD(4), 0, 0, 0, -pkin(3) * t12, -pkin(3) * t2; 0, 0, 0, 0, 0, t7, t13, t2, -t12, 0, -pkin(5) * t2 - t9, pkin(5) * t12 - t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t7, -t13, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
