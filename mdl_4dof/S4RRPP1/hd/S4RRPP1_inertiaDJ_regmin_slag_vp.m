% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:40
% EndTime: 2018-11-14 13:51:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (40->19), mult. (135->36), div. (0->0), fcn. (81->4), ass. (0->18)
t14 = 2 * qJD(4);
t10 = sin(pkin(6));
t11 = cos(pkin(6));
t12 = sin(qJ(2));
t19 = t11 * t12;
t13 = cos(qJ(2));
t9 = t13 * pkin(1) + pkin(2);
t20 = pkin(1) * t19 + t10 * t9;
t18 = pkin(1) * qJD(2);
t17 = t12 * t18;
t16 = t13 * t18;
t4 = -t10 * t17 + t11 * t16;
t15 = -t10 * t12 * pkin(1) + t11 * t9;
t8 = t10 * pkin(2) + qJ(4);
t3 = (t10 * t13 + t19) * t18;
t2 = qJ(4) + t20;
t1 = qJD(4) + t4;
t5 = [0, 0, 0, 0, -0.2e1 * t17, -0.2e1 * t16, -0.2e1 * t15 * t3 + 0.2e1 * t20 * t4, -0.2e1 * t3, 0.2e1 * t1, 0.2e1 * t2 * t1 + 0.2e1 * (-pkin(3) - t15) * t3; 0, 0, 0, 0, -t17, -t16 (t10 * t4 - t11 * t3) * pkin(2), -t3, t14 + t4, t1 * t8 + t2 * qJD(4) + t3 * (-t11 * pkin(2) - pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, t14, t8 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
