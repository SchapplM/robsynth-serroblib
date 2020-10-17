% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP7
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (37->17), mult. (59->29), div. (0->0), fcn. (42->2), ass. (0->13)
t10 = cos(qJ(3));
t14 = -0.2e1 * t10;
t13 = 2 * qJ(2);
t7 = t10 ^ 2;
t9 = sin(qJ(3));
t4 = t9 ^ 2 + t7;
t11 = -pkin(1) - pkin(5);
t12 = t9 * t11;
t3 = t10 * pkin(3) + t9 * qJ(4);
t5 = t10 * t11;
t2 = t9 * pkin(3) - t10 * qJ(4) + qJ(2);
t1 = t4 * t11;
t6 = [1, 0, 0, -2 * pkin(1), t13, pkin(1) ^ 2 + qJ(2) ^ 2, t7, t9 * t14, 0, 0, 0, t9 * t13, t10 * t13, 0.2e1 * t2 * t9, -0.2e1 * t1, t2 * t14, t4 * t11 ^ 2 + t2 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, t1; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t5, -t12, t5, -t3, t12, t3 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t10, 0, t9, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
