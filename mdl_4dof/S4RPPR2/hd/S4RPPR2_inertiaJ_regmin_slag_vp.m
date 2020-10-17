% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:11:04
% EndTime: 2019-05-04 19:11:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (46->18), mult. (62->26), div. (0->0), fcn. (58->4), ass. (0->13)
t14 = -pkin(1) - pkin(2);
t10 = sin(pkin(6));
t11 = cos(pkin(6));
t6 = t10 * qJ(2) - t11 * t14;
t13 = cos(qJ(4));
t12 = sin(qJ(4));
t8 = t11 * qJ(2) + t10 * t14;
t5 = -pkin(3) - t6;
t4 = t13 * t10 + t12 * t11;
t3 = t12 * t10 - t13 * t11;
t2 = t12 * t5 + t13 * t8;
t1 = t12 * t8 - t13 * t5;
t7 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t6, 0.2e1 * t8, t6 ^ 2 + t8 ^ 2, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, -1, 0, -pkin(1), -t11, t10, t8 * t10 - t6 * t11, 0, t3, t4; 0, 0, 0, 0, 0, 1, 0, 0, t10 ^ 2 + t11 ^ 2, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t7;
