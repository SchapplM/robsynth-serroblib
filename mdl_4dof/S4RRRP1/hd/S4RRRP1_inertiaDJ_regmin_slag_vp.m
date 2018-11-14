% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:31
% EndTime: 2018-11-14 13:54:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (68->24), mult. (208->43), div. (0->0), fcn. (123->4), ass. (0->22)
t11 = cos(qJ(3));
t23 = t11 * pkin(2);
t10 = sin(qJ(2));
t22 = t10 * t11;
t21 = pkin(1) * qJD(2);
t9 = sin(qJ(3));
t20 = qJD(3) * t9;
t19 = qJD(3) * t11;
t18 = t9 * t10 * pkin(1);
t17 = pkin(2) * t20;
t16 = t10 * t21;
t12 = cos(qJ(2));
t15 = t12 * t21;
t14 = pkin(2) * t19;
t8 = t12 * pkin(1) + pkin(2);
t1 = -t8 * t19 - t11 * t15 + (qJD(2) + qJD(3)) * t18;
t13 = (-t10 * t19 + (-t12 * t9 - t22) * qJD(2)) * pkin(1);
t7 = pkin(3) + t23;
t4 = pkin(1) * t22 + t9 * t8;
t3 = t11 * t8 + pkin(3) - t18;
t2 = -t8 * t20 + t13;
t5 = [0, 0, 0, 0, -0.2e1 * t16, -0.2e1 * t15, 0, 0.2e1 * t2, 0.2e1 * t1, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2; 0, 0, 0, 0, -t16, -t15, 0 (-pkin(2) - t8) * t20 + t13, t1 - t14, t2 * t7 + (-t1 * t9 + (t11 * t4 - t3 * t9) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17, -0.2e1 * t14, 0.2e1 * (-t7 + t23) * t17; 0, 0, 0, 0, 0, 0, 0, t2, t1, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, -t17, -t14, -pkin(3) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
