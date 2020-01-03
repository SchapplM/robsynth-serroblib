% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (24->16), mult. (97->30), div. (0->0), fcn. (52->4), ass. (0->18)
t10 = sin(qJ(3));
t19 = pkin(2) * qJD(3);
t15 = t10 * t19;
t12 = cos(qJ(3));
t7 = -t12 * pkin(2) - pkin(3);
t11 = cos(qJ(4));
t8 = t11 * qJD(4);
t9 = sin(qJ(4));
t20 = t9 * t15 + t7 * t8;
t18 = t9 * qJD(4);
t17 = pkin(3) * t18;
t16 = pkin(3) * t8;
t14 = t12 * t19;
t13 = -t11 * t15 + t7 * t18;
t6 = t10 * pkin(2) + pkin(6);
t5 = 0.2e1 * t9 * t8;
t1 = 0.2e1 * (t11 ^ 2 - t9 ^ 2) * qJD(4);
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t15, -0.2e1 * t14, t5, t1, 0, 0, 0, 0.2e1 * t13, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t15, -t14, t5, t1, 0, 0, 0, t13 - t17, -t16 + t20; 0, 0, 0, 0, 0, 0, 0, t5, t1, 0, 0, 0, -0.2e1 * t17, -0.2e1 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t18, 0, -t9 * t14 - t6 * t8, -t11 * t14 + t6 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t18, 0, -pkin(6) * t8, pkin(6) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
