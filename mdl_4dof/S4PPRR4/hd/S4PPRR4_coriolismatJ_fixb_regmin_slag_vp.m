% Calculate minimal parameter regressor of coriolis matrix for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x12]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:41
% EndTime: 2019-12-31 16:18:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (61->17), mult. (165->31), div. (0->0), fcn. (172->6), ass. (0->20)
t26 = cos(qJ(3));
t25 = cos(pkin(7));
t16 = cos(qJ(4));
t24 = qJD(3) * t16;
t13 = sin(pkin(7));
t15 = sin(qJ(3));
t10 = t26 * t13 + t15 * t25;
t23 = t10 * qJD(3);
t14 = sin(qJ(4));
t11 = -t14 ^ 2 + t16 ^ 2;
t22 = t11 * qJD(3);
t21 = t14 * qJD(4);
t20 = t16 * qJD(4);
t19 = pkin(3) * t14 * qJD(3);
t18 = pkin(3) * t24;
t17 = t14 * t24;
t9 = t15 * t13 - t26 * t25;
t4 = t9 * t16;
t3 = t9 * t14;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t23, t9 * qJD(3), 0, 0, 0, 0, 0, t3 * qJD(4) - t16 * t23, t4 * qJD(4) + t14 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(3) - t10 * t20, t4 * qJD(3) + t10 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t14 * t20, t11 * qJD(4), 0, 0, 0, -pkin(3) * t21, -pkin(3) * t20; 0, 0, 0, 0, 0, t17, t22, t20, -t21, 0, -pkin(5) * t20 - t19, pkin(5) * t21 - t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t17, -t22, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
