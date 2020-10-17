% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (31->17), mult. (60->26), div. (0->0), fcn. (46->4), ass. (0->15)
t16 = sin(qJ(2)) * pkin(1);
t3 = qJ(3) + t16;
t17 = 0.2e1 * t3;
t12 = 0.2e1 * qJ(3);
t15 = cos(qJ(2)) * pkin(1);
t14 = qJ(3) + t3;
t5 = -pkin(2) - t15;
t13 = -0.2e1 * pkin(2);
t11 = -pkin(2) - pkin(6);
t9 = cos(qJ(4));
t7 = sin(qJ(4));
t6 = t9 ^ 2;
t2 = -0.2e1 * t9 * t7;
t1 = -pkin(6) + t5;
t4 = [1, 0, 0, 1, 0.2e1 * t15, -0.2e1 * t16, 0.2e1 * t5, t17, t3 ^ 2 + t5 ^ 2, t6, t2, 0, 0, 0, t7 * t17, t9 * t17; 0, 0, 0, 1, t15, -t16, t13 - t15, t12 + t16, -t5 * pkin(2) + t3 * qJ(3), t6, t2, 0, 0, 0, t14 * t7, t14 * t9; 0, 0, 0, 1, 0, 0, t13, t12, pkin(2) ^ 2 + qJ(3) ^ 2, t6, t2, 0, 0, 0, t7 * t12, t9 * t12; 0, 0, 0, 0, 0, 0, 1, 0, t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t7, 0, t9 * t1, -t7 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t7, 0, t9 * t11, -t7 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
