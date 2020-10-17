% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (26->16), mult. (49->25), div. (0->0), fcn. (71->6), ass. (0->11)
t11 = 2 * qJ(4);
t10 = -pkin(3) - pkin(6);
t9 = cos(qJ(3));
t8 = cos(qJ(5));
t7 = sin(qJ(3));
t6 = sin(qJ(5));
t5 = cos(pkin(8));
t4 = sin(pkin(8));
t2 = t9 * t4 + t7 * t5;
t1 = t7 * t4 - t9 * t5;
t3 = [1, t4 ^ 2 + t5 ^ 2, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t1, -t2, t1, t2, -t1 * pkin(3) + t2 * qJ(4), 0, 0, 0, 0, 0, t2 * t6, t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, -0.2e1 * pkin(3), t11, pkin(3) ^ 2 + (qJ(4) ^ 2), t8 ^ 2, -0.2e1 * t8 * t6, 0, 0, 0, t6 * t11, t8 * t11; 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t1, -t6 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t6, 0, t8 * t10, -t6 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
