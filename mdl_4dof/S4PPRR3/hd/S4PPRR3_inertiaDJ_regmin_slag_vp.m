% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (9->9), mult. (41->23), div. (0->0), fcn. (30->4), ass. (0->10)
t2 = sin(qJ(4));
t9 = t2 * qJD(4);
t3 = sin(qJ(3));
t8 = t3 * qJD(3);
t4 = cos(qJ(4));
t1 = t4 * qJD(4);
t5 = cos(qJ(3));
t7 = t5 * qJD(3);
t6 = -0.2e1 * pkin(3) * qJD(4);
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t8, -t7, 0, 0, 0, 0, 0, -t4 * t8 - t5 * t9, -t5 * t1 + t2 * t8; 0, 0, 0, 0, 0, 0.2e1 * t2 * t1, 0.2e1 * (-t2 ^ 2 + t4 ^ 2) * qJD(4), 0, 0, 0, t2 * t6, t4 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t1 - t2 * t7, t3 * t9 - t4 * t7; 0, 0, 0, 0, 0, 0, 0, t1, -t9, 0, -pkin(5) * t1, pkin(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
