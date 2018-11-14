% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPPR2
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
% MMD_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:42
% EndTime: 2018-11-14 13:47:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (41->14), mult. (89->27), div. (0->0), fcn. (64->4), ass. (0->16)
t16 = 2 * qJD(2);
t12 = -pkin(1) - pkin(2);
t8 = sin(pkin(6));
t9 = cos(pkin(6));
t15 = -t8 * qJ(2) + t9 * t12;
t10 = sin(qJ(4));
t11 = cos(qJ(4));
t14 = t10 * t9 + t11 * t8;
t13 = t10 * t8 - t11 * t9;
t6 = t9 * qJ(2) + t8 * t12;
t5 = -pkin(3) + t15;
t4 = t14 * qJD(4);
t3 = t13 * qJD(4);
t2 = (t10 * t5 + t11 * t6) * qJD(4) + t14 * qJD(2);
t1 = (t10 * t6 - t11 * t5) * qJD(4) + t13 * qJD(2);
t7 = [0, 0, 0, 0, t16, qJ(2) * t16, t8 * t16, t9 * t16 (-t15 * t8 + t6 * t9) * t16, 0, 0.2e1 * t2, -0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
