% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:08
% EndTime: 2021-01-14 22:36:08
% DurationCPUTime: 0.10s
% Computational Cost: add. (35->24), mult. (73->37), div. (0->0), fcn. (74->4), ass. (0->17)
t10 = cos(qJ(3));
t18 = 0.2e1 * t10;
t8 = sin(qJ(3));
t9 = sin(qJ(2));
t17 = t8 * t9;
t5 = t8 ^ 2;
t16 = t10 ^ 2 + t5;
t15 = t10 * t9;
t11 = cos(qJ(2));
t14 = t11 * t8;
t13 = qJ(4) + pkin(5);
t1 = t13 * t8;
t2 = t13 * t10;
t12 = t1 * t8 + t2 * t10;
t4 = -t10 * pkin(3) - pkin(2);
t3 = t11 * t10;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t9 ^ 2 + t11 ^ 2; 0, 0, t11, -t9, 0, 0, 0, 0, 0, t3, -t14, t3, -t14, t16 * t9, -t11 * t4 + t12 * t9; 0, 1, 0, 0, t5, t8 * t18, 0, 0, 0, pkin(2) * t18, -0.2e1 * pkin(2) * t8, -0.2e1 * t4 * t10, 0.2e1 * t4 * t8, 0.2e1 * t12, t1 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t15, -t17, -t15, 0, -pkin(3) * t17; 0, 0, 0, 0, 0, 0, t8, t10, 0, -t8 * pkin(5), -t10 * pkin(5), -t1, -t2, -t8 * pkin(3), -t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t8, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
