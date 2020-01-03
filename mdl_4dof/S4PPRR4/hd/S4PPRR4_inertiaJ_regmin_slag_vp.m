% Calculate minimal parameter regressor of joint inertia matrix for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t7 = cos(qJ(4));
t9 = 0.2e1 * t7;
t8 = cos(qJ(3));
t6 = sin(qJ(3));
t5 = sin(qJ(4));
t4 = cos(pkin(7));
t3 = sin(pkin(7));
t2 = t8 * t3 + t6 * t4;
t1 = t6 * t3 - t8 * t4;
t10 = [1, t3 ^ 2 + t4 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t1, -t2, 0, 0, 0, 0, 0, -t1 * t7, t1 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, t5 ^ 2, t5 * t9, 0, 0, 0, pkin(3) * t9, -0.2e1 * pkin(3) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t2, -t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t5; 0, 0, 0, 0, 0, 0, 0, t5, t7, 0, -t5 * pkin(5), -t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
