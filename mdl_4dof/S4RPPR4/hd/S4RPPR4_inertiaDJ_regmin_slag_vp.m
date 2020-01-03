% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPPR4
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
% MMD_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (17->13), mult. (37->22), div. (0->0), fcn. (21->4), ass. (0->8)
t7 = 2 * qJD(3);
t3 = sin(qJ(4));
t6 = t3 * qJD(4);
t4 = cos(qJ(4));
t5 = t4 * qJD(4);
t2 = sin(pkin(6)) * pkin(1) + qJ(3);
t1 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5);
t8 = [0, 0, 0, 0, 0, t7, t2 * t7, -0.2e1 * t3 * t5, 0.2e1 * (t3 ^ 2 - t4 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t3 + 0.2e1 * t2 * t5, 0.2e1 * qJD(3) * t4 - 0.2e1 * t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5, 0, -t1 * t6, -t1 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
