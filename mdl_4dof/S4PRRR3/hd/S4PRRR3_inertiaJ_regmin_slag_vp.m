% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR3
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
% MM_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t7 = cos(qJ(4));
t12 = 0.2e1 * t7;
t11 = sin(qJ(3)) * pkin(2);
t10 = cos(qJ(3)) * pkin(2);
t3 = -pkin(3) - t10;
t9 = pkin(3) - t3;
t5 = sin(qJ(4));
t4 = t5 ^ 2;
t2 = pkin(6) + t11;
t1 = t5 * t12;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t10, -0.2e1 * t11, t4, t1, 0, 0, 0, -0.2e1 * t3 * t7, 0.2e1 * t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t10, -t11, t4, t1, 0, 0, 0, t9 * t7, -t9 * t5; 0, 0, 0, 0, 1, 0, 0, t4, t1, 0, 0, 0, pkin(3) * t12, -0.2e1 * pkin(3) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t7, 0, -t5 * t2, -t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t7, 0, -t5 * pkin(6), -t7 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
