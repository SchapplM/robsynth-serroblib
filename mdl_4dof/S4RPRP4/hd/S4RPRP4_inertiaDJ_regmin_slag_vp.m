% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRP4
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
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (41->17), mult. (100->35), div. (0->0), fcn. (60->4), ass. (0->15)
t16 = 2 * qJD(4);
t8 = sin(qJ(3));
t15 = t8 * qJD(3);
t9 = cos(qJ(3));
t6 = t9 * qJD(3);
t14 = 0.2e1 * t6;
t3 = sin(pkin(6)) * pkin(1) + pkin(5);
t13 = t3 * t15;
t12 = t3 * t6;
t4 = -cos(pkin(6)) * pkin(1) - pkin(2);
t11 = -t9 * pkin(3) - t8 * qJ(4);
t10 = t11 * qJD(3) + t9 * qJD(4);
t2 = t11 + t4;
t1 = -pkin(3) * t15 + qJ(4) * t6 + t8 * qJD(4);
t5 = [0, 0, 0, 0, t8 * t14, 0.2e1 * (-t8 ^ 2 + t9 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t4 * t15, t4 * t14, 0.2e1 * t1 * t9 + 0.2e1 * t2 * t15, 0, 0.2e1 * t1 * t8 - 0.2e1 * t2 * t6, -0.2e1 * t2 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t6, -t15, 0, -t12, t13, -t12, t10, -t13, t10 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t6, -t15, 0, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, qJ(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
