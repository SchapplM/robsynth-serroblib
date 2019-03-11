% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (43->19), mult. (120->30), div. (0->0), fcn. (103->4), ass. (0->17)
t18 = qJD(2) + qJD(3);
t8 = sin(qJ(3));
t9 = sin(qJ(2));
t17 = t8 * t9;
t10 = cos(qJ(3));
t16 = t10 * pkin(2);
t15 = qJD(3) * t10;
t11 = cos(qJ(2));
t14 = t11 * qJD(2);
t13 = qJD(3) * t8 * pkin(2);
t12 = pkin(2) * t15;
t4 = t10 * t9 + t8 * t11;
t7 = pkin(3) + t16;
t3 = t10 * t11 - t17;
t2 = t18 * t4;
t1 = -t10 * t14 - t11 * t15 + t18 * t17;
t5 = [0, 0, 0, 0, 0, 0, 0, -0.2e1 * t4 * t1 - 0.2e1 * t3 * t2; 0, 0, -t9 * qJD(2), -t14, 0, -t2, t1, -t2 * t7 + (-t1 * t8 + (t10 * t4 - t3 * t8) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, -0.2e1 * t13, -0.2e1 * t12, 0.2e1 * (-t7 + t16) * t13; 0, 0, 0, 0, 0, -t2, t1, -t2 * pkin(3); 0, 0, 0, 0, 0, -t13, -t12, -pkin(3) * t13; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
