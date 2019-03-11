% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:43
% EndTime: 2019-03-08 18:29:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (42->14), mult. (110->20), div. (0->0), fcn. (64->4), ass. (0->13)
t10 = cos(qJ(3));
t15 = pkin(1) * sin(pkin(6));
t7 = cos(pkin(6)) * pkin(1) + pkin(2);
t9 = sin(qJ(3));
t16 = -t10 * t7 + t9 * t15;
t11 = 2 * qJD(4);
t3 = t16 * qJD(3);
t12 = t10 * t15 + t9 * t7;
t5 = qJ(4) + t12;
t4 = t12 * qJD(3);
t2 = qJD(4) - t3;
t1 = 0.2e1 * t4;
t6 = [0, 0, 0, 0, 0, -t1, 0.2e1 * t3, -t1, 0.2e1 * t2, 0.2e1 * t5 * t2 + 0.2e1 * (-pkin(3) + t16) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t4, t3, -t4, t11 - t3, -t4 * pkin(3) + t2 * qJ(4) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t11, qJ(4) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
