% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:23
% EndTime: 2018-11-14 14:10:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (9->8), mult. (24->8), div. (0->0), fcn. (14->2), ass. (0->7)
t8 = 2 * qJD(3);
t4 = cos(qJ(2)) * qJD(2);
t6 = sin(qJ(2));
t9 = qJ(3) * t4 + t6 * qJD(3);
t3 = t6 * qJD(2);
t5 = qJ(3) * t8;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t3, -t4, -t3, t4, -pkin(2) * t3 + t9, -t3, t4 (-pkin(2) - pkin(3)) * t3 + t9; 0, 0, 0, 0, 0, t8, t5, 0, t8, t5; 0, 0, 0, 0, 0, 0, t3, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
