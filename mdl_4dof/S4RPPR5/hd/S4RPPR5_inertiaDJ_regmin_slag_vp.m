% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:48
% EndTime: 2019-12-31 16:39:49
% DurationCPUTime: 0.11s
% Computational Cost: add. (26->15), mult. (66->37), div. (0->0), fcn. (44->4), ass. (0->15)
t16 = 2 * qJD(2);
t5 = sin(pkin(6));
t6 = cos(pkin(6));
t9 = -pkin(1) - pkin(2);
t15 = t6 * qJ(2) + t5 * t9;
t14 = t5 * qJD(2);
t13 = t6 * qJD(2);
t7 = sin(qJ(4));
t12 = t7 * qJD(4);
t8 = cos(qJ(4));
t11 = t8 * qJD(4);
t10 = -t5 * qJ(2) + t6 * t9;
t2 = -pkin(5) + t15;
t1 = pkin(3) - t10;
t3 = [0, 0, 0, 0, t16, qJ(2) * t16, 0.2e1 * t14, 0.2e1 * t13, (-t10 * t5 + t15 * t6) * t16, 0.2e1 * t7 * t11, 0.2e1 * (-t7 ^ 2 + t8 ^ 2) * qJD(4), 0, 0, 0, -0.2e1 * t1 * t12 + 0.2e1 * t8 * t14, -0.2e1 * t1 * t11 - 0.2e1 * t7 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t12, t6 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t12, 0, -t2 * t11 - t7 * t13, t2 * t12 - t8 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t11, t5 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
