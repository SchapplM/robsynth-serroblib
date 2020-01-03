% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (19->12), mult. (67->27), div. (0->0), fcn. (60->6), ass. (0->14)
t11 = cos(qJ(3));
t6 = sin(pkin(7));
t7 = cos(pkin(7));
t9 = sin(qJ(3));
t3 = -t11 * t7 + t9 * t6;
t8 = sin(qJ(4));
t14 = t8 * qJD(4);
t10 = cos(qJ(4));
t13 = t10 * qJD(4);
t12 = -0.2e1 * pkin(3) * qJD(4);
t4 = t11 * t6 + t9 * t7;
t2 = t4 * qJD(3);
t1 = t3 * qJD(3);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, -t2 * t10 + t3 * t14, t3 * t13 + t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0.2e1 * t8 * t13, 0.2e1 * (t10 ^ 2 - t8 ^ 2) * qJD(4), 0, 0, 0, t8 * t12, t10 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t1 - t4 * t13, t10 * t1 + t4 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t13; 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, -pkin(5) * t13, pkin(5) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
