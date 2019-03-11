% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t5 = sin(qJ(3));
t9 = t5 * pkin(2);
t8 = cos(qJ(2));
t7 = cos(qJ(3));
t6 = sin(qJ(2));
t4 = t7 * pkin(2);
t3 = t4 + pkin(3);
t2 = t5 * t8 + t7 * t6;
t1 = -t5 * t6 + t7 * t8;
t10 = [1, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2; 0, 0, t8, -t6, 0, t1, -t2, t1 * t3 + t2 * t9; 0, 1, 0, 0, 1, 0.2e1 * t4, -0.2e1 * t9, t5 ^ 2 * pkin(2) ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t1, -t2, t1 * pkin(3); 0, 0, 0, 0, 1, t4, -t9, t3 * pkin(3); 0, 0, 0, 0, 1, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
