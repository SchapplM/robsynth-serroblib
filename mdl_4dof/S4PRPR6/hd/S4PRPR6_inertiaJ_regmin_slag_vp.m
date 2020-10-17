% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRPR6
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (36->19), mult. (90->41), div. (0->0), fcn. (101->6), ass. (0->15)
t11 = cos(pkin(7));
t20 = -0.2e1 * t11 * pkin(3) - (2 * pkin(2));
t10 = sin(pkin(7));
t19 = t10 ^ 2 + t11 ^ 2;
t18 = pkin(5) + qJ(3);
t17 = t19 * qJ(3);
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t2 = t14 * t10 + t12 * t11;
t1 = t12 * t10 - t14 * t11;
t15 = cos(qJ(2));
t13 = sin(qJ(2));
t4 = t18 * t11;
t3 = t18 * t10;
t5 = [1, 0, 0, 0, 0, 0, 0, t19 * t13 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t15, -t13, t15 * t11, -t15 * t10, t19 * t13, t15 * pkin(2) + t13 * t17, 0, 0, 0, 0, 0, -t15 * t1, -t15 * t2; 0, 1, 0, 0, 0.2e1 * pkin(2) * t11, -0.2e1 * pkin(2) * t10, 0.2e1 * t17, t19 * qJ(3) ^ 2 + (pkin(2) ^ 2), t2 ^ 2, -0.2e1 * t2 * t1, 0, 0, 0, t1 * t20, t2 * t20; 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t11, t10, 0, -pkin(2), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t13, t1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, 0, -t12 * t4 - t14 * t3, t12 * t3 - t14 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
