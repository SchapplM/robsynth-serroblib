% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:31
% EndTime: 2019-12-31 17:33:32
% DurationCPUTime: 0.16s
% Computational Cost: add. (29->27), mult. (71->35), div. (0->0), fcn. (56->4), ass. (0->19)
t3 = sin(qJ(5));
t5 = cos(qJ(5));
t1 = t3 ^ 2 - t5 ^ 2;
t18 = t1 * qJD(3);
t17 = t3 * qJD(3);
t16 = t3 * qJD(5);
t4 = sin(qJ(3));
t2 = t4 * qJD(3);
t15 = t5 * qJD(3);
t14 = t5 * qJD(5);
t6 = cos(qJ(3));
t13 = t6 * qJD(3);
t12 = qJ(4) * qJD(5);
t11 = qJD(3) * qJ(4);
t10 = t3 * t15;
t9 = t3 * t11;
t8 = t5 * t11;
t7 = -pkin(3) - pkin(6);
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t2, -t13, t2, t13, (-t4 * pkin(3) + t6 * qJ(4)) * qJD(3) + t4 * qJD(4), 0, 0, 0, 0, 0, t3 * t13 + t4 * t14, t5 * t13 - t4 * t16; 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t15 + t6 * t16, t6 * t14 - t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t3 * t14, t1 * qJD(5), 0, 0, 0, qJD(4) * t3 + t5 * t12, qJD(4) * t5 - t3 * t12; 0, 0, 0, 0, 0, 0, qJD(3), t11, 0, 0, 0, 0, 0, t17, t15; 0, 0, 0, 0, 0, 0, 0, 0, -t10, t18, -t16, -t14, 0, -t7 * t16 + t8, -t7 * t14 - t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -qJD(3), -t11, 0, 0, 0, 0, 0, -t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t10, -t18, 0, 0, 0, -t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t19;
