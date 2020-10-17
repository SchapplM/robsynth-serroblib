% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:45
% DurationCPUTime: 0.12s
% Computational Cost: add. (29->21), mult. (63->28), div. (0->0), fcn. (85->6), ass. (0->19)
t12 = cos(qJ(5));
t19 = 0.2e1 * t12;
t13 = cos(qJ(4));
t16 = t13 * pkin(3);
t7 = -pkin(4) - t16;
t18 = pkin(4) - t7;
t10 = sin(qJ(4));
t17 = t10 * pkin(3);
t11 = sin(qJ(3));
t14 = cos(qJ(3));
t2 = t10 * t11 - t13 * t14;
t15 = t2 * t12;
t9 = sin(qJ(5));
t8 = t9 ^ 2;
t6 = pkin(7) + t17;
t4 = t9 * t19;
t3 = t10 * t14 + t13 * t11;
t1 = t2 * t9;
t5 = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t14, -t11, 0, -t2, -t3, 0, 0, 0, 0, 0, -t15, t1; 0, 0, 1, 0, 0, 1, 0.2e1 * t16, -0.2e1 * t17, t8, t4, 0, 0, 0, -0.2e1 * t7 * t12, 0.2e1 * t7 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t2, -t3, 0, 0, 0, 0, 0, -t15, t1; 0, 0, 0, 0, 0, 1, t16, -t17, t8, t4, 0, 0, 0, t18 * t12, -t18 * t9; 0, 0, 0, 0, 0, 1, 0, 0, t8, t4, 0, 0, 0, pkin(4) * t19, -0.2e1 * pkin(4) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t3, -t12 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t12, 0, -t9 * t6, -t12 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t12, 0, -t9 * pkin(7), -t12 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
