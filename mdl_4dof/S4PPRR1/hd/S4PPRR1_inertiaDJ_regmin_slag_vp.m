% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PPRR1
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
% MMD_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4PPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:19
% EndTime: 2018-11-14 13:40:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (16->7), mult. (44->13), div. (0->0), fcn. (38->4), ass. (0->11)
t12 = qJD(3) + qJD(4);
t11 = pkin(3) * qJD(4);
t3 = sin(qJ(4));
t10 = t3 * t11;
t5 = cos(qJ(4));
t9 = t5 * t11;
t6 = cos(qJ(3));
t4 = sin(qJ(3));
t2 = t12 * (-t3 * t6 - t4 * t5);
t1 = t12 * (t3 * t4 - t5 * t6);
t7 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t4 * qJD(3), -t6 * qJD(3), 0, t2, t1; 0, 0, 0, 0, 0, 0, -0.2e1 * t10, -0.2e1 * t9; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
