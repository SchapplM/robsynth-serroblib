% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (31->16), mult. (74->30), div. (0->0), fcn. (93->4), ass. (0->16)
t12 = cos(qJ(3));
t17 = -0.2e1 * t12 * pkin(3) - (2 * pkin(2));
t16 = 0.2e1 * t12;
t9 = sin(qJ(4));
t15 = t9 * pkin(3);
t14 = pkin(5) + pkin(6);
t11 = cos(qJ(4));
t13 = t11 * pkin(3);
t10 = sin(qJ(3));
t6 = t14 * t12;
t5 = t14 * t10;
t4 = t11 * t10 + t9 * t12;
t3 = t9 * t10 - t11 * t12;
t2 = -t11 * t6 + t9 * t5;
t1 = -t11 * t5 - t9 * t6;
t7 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, t10 ^ 2, t10 * t16, 0, 0, 0, pkin(2) * t16, -0.2e1 * pkin(2) * t10, t4 ^ 2, -0.2e1 * t4 * t3, 0, 0, 0, t3 * t17, t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, t10, t12, 0, -t10 * pkin(5), -t12 * pkin(5), 0, 0, t4, -t3, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t13, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
