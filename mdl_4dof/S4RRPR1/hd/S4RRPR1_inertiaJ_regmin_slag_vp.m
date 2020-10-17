% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:21:41
% EndTime: 2019-05-04 19:21:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (62->18), mult. (126->34), div. (0->0), fcn. (108->6), ass. (0->20)
t14 = sin(pkin(7));
t23 = pkin(2) * t14;
t22 = sin(qJ(2)) * pkin(1);
t13 = cos(qJ(2)) * pkin(1);
t12 = t13 + pkin(2);
t15 = cos(pkin(7));
t7 = t14 * t12 + t15 * t22;
t21 = -t7 - t23;
t5 = t15 * t12 - t14 * t22;
t18 = cos(qJ(4));
t16 = sin(qJ(4));
t11 = t15 * pkin(2) + pkin(3);
t10 = t18 * t11;
t8 = -t16 * t11 - t18 * t23;
t6 = -t16 * t23 + t10;
t4 = pkin(3) + t5;
t3 = t18 * t4;
t2 = -t16 * t4 - t18 * t7;
t1 = -t16 * t7 + t3;
t9 = [1, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t22, t5 ^ 2 + t7 ^ 2, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 1, t13, -t22 (t14 * t7 + t15 * t5) * pkin(2), 1, t21 * t16 + t10 + t3, t21 * t18 + (-t11 - t4) * t16; 0, 0, 0, 1, 0, 0 (t14 ^ 2 + t15 ^ 2) * pkin(2) ^ 2, 1, 0.2e1 * t6, 0.2e1 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, t6, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
