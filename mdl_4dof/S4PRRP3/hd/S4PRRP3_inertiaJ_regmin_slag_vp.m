% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRP3
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
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:05
% EndTime: 2021-01-14 22:27:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (24->13), mult. (51->29), div. (0->0), fcn. (47->2), ass. (0->10)
t6 = cos(qJ(3));
t9 = 0.2e1 * t6;
t8 = t6 * pkin(3);
t7 = -qJ(4) - pkin(5);
t5 = sin(qJ(3));
t4 = t5 ^ 2;
t3 = -pkin(2) - t8;
t2 = t7 * t6;
t1 = t7 * t5;
t10 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 ^ 2 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t1 - t5 * t2; 0, 1, 0, 0, t4, t5 * t9, 0, 0, 0, pkin(2) * t9, -0.2e1 * pkin(2) * t5, -0.2e1 * t3 * t6, 0.2e1 * t3 * t5, -0.2e1 * t1 * t5 - 0.2e1 * t2 * t6, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, 0, t8; 0, 0, 0, 0, 0, 0, t5, t6, 0, -t5 * pkin(5), -t6 * pkin(5), t1, t2, -t5 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t10;
