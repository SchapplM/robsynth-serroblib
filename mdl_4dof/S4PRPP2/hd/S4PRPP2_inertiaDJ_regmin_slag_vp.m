% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:00
% EndTime: 2019-03-08 18:19:00
% DurationCPUTime: 0.07s
% Computational Cost: add. (26->11), mult. (77->22), div. (0->0), fcn. (71->4), ass. (0->12)
t13 = 2 * qJD(4);
t10 = cos(qJ(2));
t7 = sin(pkin(5));
t8 = cos(pkin(5));
t9 = sin(qJ(2));
t12 = t10 * t8 - t7 * t9;
t5 = t10 * t7 + t8 * t9;
t2 = t5 * qJD(2);
t3 = t12 * qJD(2);
t11 = -0.2e1 * t2 * t12 + 0.2e1 * t5 * t3;
t6 = pkin(2) * t7 + qJ(4);
t1 = [0, 0, 0, 0, t11, 0, 0, t11; 0, 0, -t9 * qJD(2), -t10 * qJD(2) (-t2 * t8 + t3 * t7) * pkin(2), -t2, t3, t3 * t6 + t5 * qJD(4) + t2 * (-pkin(2) * t8 - pkin(3)); 0, 0, 0, 0, 0, 0, t13, t6 * t13; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
