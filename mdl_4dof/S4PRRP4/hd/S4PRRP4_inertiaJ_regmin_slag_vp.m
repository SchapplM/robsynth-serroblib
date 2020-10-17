% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRP4
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (20->12), mult. (49->26), div. (0->0), fcn. (38->2), ass. (0->12)
t4 = sin(qJ(3));
t13 = -0.2e1 * t4;
t5 = cos(qJ(3));
t12 = 0.2e1 * t5;
t11 = t4 * pkin(5);
t10 = t5 * pkin(5);
t2 = t4 ^ 2;
t9 = t5 ^ 2 + t2;
t8 = t5 * pkin(3) + t4 * qJ(4);
t7 = -t4 * pkin(3) + t5 * qJ(4);
t1 = -pkin(2) - t8;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, t2, t4 * t12, 0, 0, 0, pkin(2) * t12, pkin(2) * t13, -0.2e1 * t1 * t5, 0.2e1 * t9 * pkin(5), t1 * t13, t9 * pkin(5) ^ 2 + t1 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, t5, 0, t4, t8; 0, 0, 0, 0, 0, 0, t4, t5, 0, -t11, -t10, -t11, t7, t10, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
