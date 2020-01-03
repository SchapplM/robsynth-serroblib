% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:52
% DurationCPUTime: 0.08s
% Computational Cost: add. (14->12), mult. (57->23), div. (0->0), fcn. (46->4), ass. (0->14)
t5 = sin(qJ(4));
t13 = t5 * qJD(4);
t6 = sin(qJ(3));
t12 = t6 * qJD(3);
t7 = cos(qJ(4));
t11 = t7 * qJD(4);
t8 = cos(qJ(3));
t10 = t8 * qJD(3);
t9 = -0.2e1 * pkin(3) * qJD(4);
t4 = -t8 * t11 + t5 * t12;
t3 = -t5 * t10 - t6 * t11;
t2 = t7 * t12 + t8 * t13;
t1 = -t7 * t10 + t6 * t13;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t10, t12, 0, 0, 0, 0, 0, t1, -t3; 0, 0, 0, -t12, -t10, 0, 0, 0, 0, 0, -t2, t4; 0, 0, 0, 0, 0, 0.2e1 * t5 * t11, 0.2e1 * (-t5 ^ 2 + t7 ^ 2) * qJD(4), 0, 0, 0, t5 * t9, t7 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t1; 0, 0, 0, 0, 0, 0, 0, t11, -t13, 0, -pkin(5) * t11, pkin(5) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
