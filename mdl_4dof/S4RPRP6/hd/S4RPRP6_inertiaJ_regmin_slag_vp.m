% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:40
% EndTime: 2021-01-15 10:27:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (37->17), mult. (57->28), div. (0->0), fcn. (49->2), ass. (0->14)
t8 = sin(qJ(3));
t5 = t8 * pkin(3) + qJ(2);
t13 = 0.2e1 * t5;
t12 = 0.2e1 * qJ(2);
t9 = cos(qJ(3));
t11 = t9 * pkin(3);
t10 = -pkin(1) - pkin(5);
t7 = t9 ^ 2;
t6 = t9 * t10;
t4 = t8 ^ 2 + t7;
t3 = -t9 * qJ(4) + t6;
t2 = (-qJ(4) + t10) * t8;
t1 = t2 * t8 + t3 * t9;
t14 = [1, 0, 0, -2 * pkin(1), t12, (pkin(1) ^ 2) + qJ(2) ^ 2, t7, -0.2e1 * t9 * t8, 0, 0, 0, t8 * t12, t9 * t12, t8 * t13, t9 * t13, -0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t1; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t6, -t8 * t10, t3, -t2, -t11, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, t9, -t8, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t14;
