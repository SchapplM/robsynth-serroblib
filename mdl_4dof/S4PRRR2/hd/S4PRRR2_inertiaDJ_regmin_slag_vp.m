% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:24
% EndTime: 2019-07-18 13:27:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (30->15), mult. (104->28), div. (0->0), fcn. (58->4), ass. (0->18)
t18 = pkin(1) * qJD(3);
t7 = sin(qJ(3));
t15 = t7 * t18;
t6 = sin(qJ(4));
t17 = qJD(4) * t6;
t19 = t7 * pkin(1) * t17 + t6 * t15;
t8 = cos(qJ(4));
t16 = qJD(4) * t8;
t9 = cos(qJ(3));
t14 = t9 * t18;
t13 = pkin(2) * t17;
t12 = pkin(2) * t16;
t5 = t9 * pkin(1) + pkin(2);
t11 = (-pkin(2) - t5) * qJD(4);
t10 = (-t7 * t16 + (-t6 * t9 - t7 * t8) * qJD(3)) * pkin(1);
t2 = -t5 * t17 + t10;
t1 = (-qJD(4) * t5 - t14) * t8 + t19;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t15, -0.2e1 * t14, 0, 0.2e1 * t2, 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t15, -t14, 0, t6 * t11 + t10, (t11 - t14) * t8 + t19; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13, -0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
