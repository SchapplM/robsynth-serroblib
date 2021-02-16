% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:36
% EndTime: 2021-01-15 10:20:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (42->19), mult. (74->33), div. (0->0), fcn. (66->4), ass. (0->15)
t9 = sin(qJ(3));
t15 = 0.2e1 * t9;
t10 = cos(qJ(3));
t14 = -0.2e1 * t10;
t13 = t10 * pkin(3);
t7 = sin(pkin(6));
t4 = t7 * pkin(1) + pkin(5);
t12 = qJ(4) + t4;
t8 = cos(pkin(6));
t5 = -t8 * pkin(1) - pkin(2);
t6 = t9 ^ 2;
t3 = t5 - t13;
t2 = t12 * t10;
t1 = t12 * t9;
t11 = [1, 0, 0, (t7 ^ 2 + t8 ^ 2) * pkin(1) ^ 2, t6, t10 * t15, 0, 0, 0, t5 * t14, t5 * t15, t3 * t14, t3 * t15, 0.2e1 * t1 * t9 + 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t10 + t2 * t9; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t6; 0, 0, 0, 0, 0, 0, t9, t10, 0, -t9 * t4, -t10 * t4, -t1, -t2, -t9 * pkin(3), -t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t10, -t9, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t11;
