% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:22
% EndTime: 2019-12-05 15:03:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (39->22), mult. (117->43), div. (0->0), fcn. (104->6), ass. (0->16)
t11 = cos(qJ(3));
t6 = sin(pkin(8));
t7 = cos(pkin(8));
t9 = sin(qJ(3));
t3 = -t11 * t7 + t9 * t6;
t19 = 2 * qJD(4);
t8 = sin(qJ(5));
t16 = qJD(5) * t8;
t10 = cos(qJ(5));
t15 = qJD(5) * t10;
t14 = qJD(5) * (-pkin(3) - pkin(6));
t13 = qJ(4) * qJD(5);
t4 = t11 * t6 + t9 * t7;
t2 = t4 * qJD(3);
t1 = t3 * qJD(3);
t5 = [0, 0, 0, 0, 0, 0, 0, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t2, t1, t2, -t1, -t2 * pkin(3) - t1 * qJ(4) + t4 * qJD(4), 0, 0, 0, 0, 0, -t1 * t8 + t4 * t15, -t1 * t10 - t4 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t19, qJ(4) * t19, -0.2e1 * t8 * t15, 0.2e1 * (-t10 ^ 2 + t8 ^ 2) * qJD(5), 0, 0, 0, 0.2e1 * qJD(4) * t8 + 0.2e1 * t10 * t13, 0.2e1 * qJD(4) * t10 - 0.2e1 * t8 * t13; 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t2 - t3 * t16, -t3 * t15 - t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15, 0, -t8 * t14, -t10 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
