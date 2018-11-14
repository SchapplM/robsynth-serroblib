% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:25
% EndTime: 2018-11-14 14:11:26
% DurationCPUTime: 0.09s
% Computational Cost: add. (46->15), mult. (130->35), div. (0->0), fcn. (124->6), ass. (0->17)
t10 = sin(pkin(6));
t16 = pkin(2) * t10;
t11 = cos(pkin(6));
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t8 = t10 * t15 + t11 * t13;
t7 = -t10 * t13 + t11 * t15;
t14 = cos(qJ(4));
t12 = sin(qJ(4));
t9 = t11 * pkin(2) + pkin(3);
t6 = t7 * qJD(2);
t5 = t8 * qJD(2);
t4 = (-t12 * t9 - t14 * t16) * qJD(4);
t3 = (t12 * t16 - t14 * t9) * qJD(4);
t2 = -t12 * t6 - t14 * t5 + (-t12 * t7 - t14 * t8) * qJD(4);
t1 = t12 * t5 - t14 * t6 + (t12 * t8 - t14 * t7) * qJD(4);
t17 = [0, 0, 0, 0, -0.2e1 * t7 * t5 + 0.2e1 * t8 * t6, 0, 0, 0; 0, 0, -t13 * qJD(2), -t15 * qJD(2) (t10 * t6 - t11 * t5) * pkin(2), 0, t2, t1; 0, 0, 0, 0, 0, 0, 0.2e1 * t4, 0.2e1 * t3; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t17;
