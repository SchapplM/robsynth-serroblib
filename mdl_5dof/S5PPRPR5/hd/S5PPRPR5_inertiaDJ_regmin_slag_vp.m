% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:30
% DurationCPUTime: 0.11s
% Computational Cost: add. (20->19), mult. (58->35), div. (0->0), fcn. (40->4), ass. (0->12)
t11 = 2 * qJD(4);
t2 = sin(qJ(5));
t10 = t2 * qJD(5);
t3 = sin(qJ(3));
t1 = t3 * qJD(3);
t4 = cos(qJ(5));
t9 = t4 * qJD(5);
t5 = cos(qJ(3));
t8 = t5 * qJD(3);
t7 = qJ(4) * qJD(5);
t6 = -pkin(3) - pkin(6);
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t1, -t8, t1, t8, t3 * qJD(4) + (-pkin(3) * t3 + qJ(4) * t5) * qJD(3), 0, 0, 0, 0, 0, t2 * t8 + t3 * t9, -t3 * t10 + t4 * t8; 0, 0, 0, 0, 0, 0, t11, qJ(4) * t11, -0.2e1 * t2 * t9, 0.2e1 * (t2 ^ 2 - t4 ^ 2) * qJD(5), 0, 0, 0, 0.2e1 * qJD(4) * t2 + 0.2e1 * t4 * t7, 0.2e1 * qJD(4) * t4 - 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t1 + t5 * t10, -t2 * t1 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t9, 0, -t6 * t10, -t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
