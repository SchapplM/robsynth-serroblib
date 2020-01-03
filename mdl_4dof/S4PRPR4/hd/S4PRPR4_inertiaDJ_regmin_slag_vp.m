% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRPR4
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
% MMD_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (12->11), mult. (32->21), div. (0->0), fcn. (16->2), ass. (0->8)
t7 = 2 * qJD(3);
t1 = sin(qJ(4));
t6 = t1 * qJD(4);
t2 = cos(qJ(4));
t5 = t2 * qJD(4);
t4 = qJ(3) * qJD(4);
t3 = -pkin(2) - pkin(5);
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t7, qJ(3) * t7, -0.2e1 * t1 * t5, 0.2e1 * (t1 ^ 2 - t2 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t1 + 0.2e1 * t2 * t4, 0.2e1 * qJD(3) * t2 - 0.2e1 * t1 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5, 0, -t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
