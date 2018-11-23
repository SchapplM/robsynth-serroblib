% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% MMD_reg [((2+1)*2/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S2RR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_inertiaDJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_inertiaDJ_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_inertiaDJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:26
% EndTime: 2018-11-16 16:44:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (3->3), mult. (15->10), div. (0->0), fcn. (8->2), ass. (0->5)
t1 = sin(qJ(2));
t4 = t1 * qJD(2);
t2 = cos(qJ(2));
t3 = t2 * qJD(2);
t5 = [0, 0, 0, 0.2e1 * t1 * t3, 0.2e1 * (-t1 ^ 2 + t2 ^ 2) * qJD(2), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t3, t4, 0, pkin(1) * t3, -pkin(1) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
