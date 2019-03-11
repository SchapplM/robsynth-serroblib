% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t8 = cos(pkin(6));
t5 = pkin(1) * t8 + pkin(2);
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t7 = sin(pkin(6));
t4 = pkin(1) * t7 + qJ(3);
t3 = -pkin(3) - t5;
t2 = t10 * t4 + t9 * t3;
t1 = -t10 * t3 + t9 * t4;
t6 = [1, 0, 0 (t7 ^ 2 + t8 ^ 2) * pkin(1) ^ 2, 0.2e1 * t5, 0.2e1 * t4, t4 ^ 2 + t5 ^ 2, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, -1, 0, -t5, 0, -t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -1, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
