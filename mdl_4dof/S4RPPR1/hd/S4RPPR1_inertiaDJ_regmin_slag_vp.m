% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:31
% EndTime: 2019-03-08 18:27:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (23->10), mult. (42->16), div. (0->0), fcn. (25->4), ass. (0->10)
t9 = 2 * qJD(3);
t5 = sin(qJ(4));
t8 = qJD(4) * t5;
t6 = cos(qJ(4));
t7 = qJD(4) * t6;
t4 = sin(pkin(6)) * pkin(1) + qJ(3);
t3 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t2 = t5 * qJD(3) + (t3 * t5 + t4 * t6) * qJD(4);
t1 = -t6 * qJD(3) + (-t3 * t6 + t4 * t5) * qJD(4);
t10 = [0, 0, 0, 0, 0, t9, t4 * t9, 0, 0.2e1 * t2, -0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
