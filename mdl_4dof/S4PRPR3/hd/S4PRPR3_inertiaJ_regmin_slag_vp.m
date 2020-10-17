% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:54
% EndTime: 2019-12-31 16:20:55
% DurationCPUTime: 0.09s
% Computational Cost: add. (25->13), mult. (62->29), div. (0->0), fcn. (68->4), ass. (0->12)
t10 = cos(pkin(7));
t16 = -0.2e1 * t10 * pkin(3) - (2 * pkin(2));
t9 = sin(pkin(7));
t15 = t10 ^ 2 + t9 ^ 2;
t14 = pkin(5) + qJ(3);
t12 = cos(qJ(4));
t11 = sin(qJ(4));
t4 = t14 * t10;
t3 = t14 * t9;
t2 = t11 * t10 + t12 * t9;
t1 = -t12 * t10 + t11 * t9;
t5 = [1, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0.2e1 * pkin(2) * t10, -0.2e1 * pkin(2) * t9, 0.2e1 * t15 * qJ(3), t15 * qJ(3) ^ 2 + (pkin(2) ^ 2), t2 ^ 2, -0.2e1 * t2 * t1, 0, 0, 0, t1 * t16, t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t10, t9, 0, -pkin(2), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, 0, -t11 * t4 - t12 * t3, t11 * t3 - t12 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
