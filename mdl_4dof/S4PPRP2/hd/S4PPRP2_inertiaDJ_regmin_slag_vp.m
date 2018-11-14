% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S4PPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:32
% EndTime: 2018-11-14 13:57:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->8), mult. (52->15), div. (0->0), fcn. (48->4), ass. (0->10)
t5 = sin(pkin(5));
t6 = cos(pkin(5));
t7 = sin(qJ(3));
t8 = cos(qJ(3));
t12 = t7 * t5 - t8 * t6;
t11 = 2 * qJD(4);
t3 = t8 * t5 + t7 * t6;
t2 = t3 * qJD(3);
t1 = t12 * qJD(3);
t4 = [0, 0, 0, 0, 0, 0, 0, -0.2e1 * t3 * t1 + 0.2e1 * t12 * t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t2, t1, -t2, -t1, -t2 * pkin(3) - t1 * qJ(4) + t3 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t11, qJ(4) * t11; 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
