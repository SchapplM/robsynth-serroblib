% Calculate minimal parameter regressor of coriolis matrix for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PPRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:53
% EndTime: 2019-12-31 16:19:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (20->17), mult. (69->25), div. (0->0), fcn. (58->4), ass. (0->19)
t8 = cos(qJ(4));
t18 = qJD(3) * t8;
t6 = sin(qJ(4));
t5 = -t6 ^ 2 + t8 ^ 2;
t17 = t5 * qJD(3);
t16 = t6 * qJD(4);
t7 = sin(qJ(3));
t15 = t7 * qJD(3);
t14 = t8 * qJD(4);
t9 = cos(qJ(3));
t13 = t9 * qJD(3);
t12 = pkin(3) * t6 * qJD(3);
t11 = pkin(3) * t18;
t10 = t6 * t18;
t4 = -t9 * t14 + t6 * t15;
t3 = -t6 * t13 - t7 * t14;
t2 = t8 * t15 + t9 * t16;
t1 = -t8 * t13 + t7 * t16;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t13, t15, 0, 0, 0, 0, 0, t1, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t15, -t13, 0, 0, 0, 0, 0, -t2, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t6 * t14, t5 * qJD(4), 0, 0, 0, -pkin(3) * t16, -pkin(3) * t14; 0, 0, 0, 0, 0, t10, t17, t14, -t16, 0, -pkin(5) * t14 - t12, pkin(5) * t16 - t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t10, -t17, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t19;
