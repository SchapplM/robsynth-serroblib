% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:35
% EndTime: 2018-11-14 13:50:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (84->18), mult. (206->30), div. (0->0), fcn. (142->6), ass. (0->21)
t12 = sin(qJ(3));
t15 = cos(pkin(7)) * pkin(1) + pkin(2);
t22 = cos(qJ(3));
t24 = pkin(1) * sin(pkin(7));
t26 = t12 * t24 - t22 * t15;
t7 = pkin(3) - t26;
t25 = -pkin(3) - t7;
t13 = cos(qJ(4));
t8 = t12 * t15 + t22 * t24;
t23 = t13 * t8;
t11 = sin(qJ(4));
t21 = qJD(4) * t11;
t19 = pkin(3) * t21;
t18 = qJD(4) * t13 * pkin(3);
t5 = t26 * qJD(3);
t6 = t8 * qJD(3);
t17 = t11 * t5 - t13 * t6;
t16 = t11 * t6 + t8 * t21;
t2 = (-t11 * t7 - t23) * qJD(4) + t17;
t1 = (-qJD(4) * t7 + t5) * t13 + t16;
t3 = [0, 0, 0, 0, 0, -0.2e1 * t6, 0.2e1 * t5, 0, 0.2e1 * t2, 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t6, t5, 0 (t25 * t11 - t23) * qJD(4) + t17 (t25 * qJD(4) + t5) * t13 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t19, -0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
