% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t9 = sin(qJ(3));
t15 = t9 * pkin(2);
t14 = sin(qJ(2)) * pkin(1);
t11 = cos(qJ(3));
t13 = t11 * t14;
t8 = cos(qJ(2)) * pkin(1);
t6 = t8 + pkin(2);
t2 = t11 * t6 - t9 * t14;
t7 = t11 * pkin(2);
t5 = t7 + pkin(3);
t3 = t9 * t6 + t13;
t1 = pkin(3) + t2;
t4 = [1, 0, 0, 1, 0.2e1 * t8, -0.2e1 * t14, 1, 0.2e1 * t2, -0.2e1 * t3, t1 ^ 2 + t3 ^ 2; 0, 0, 0, 1, t8, -t14, 1, t2 + t7, -t13 + (-pkin(2) - t6) * t9, t1 * t5 + t3 * t15; 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t7, -0.2e1 * t15, t9 ^ 2 * pkin(2) ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 1, t2, -t3, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 1, t7, -t15, t5 * pkin(3); 0, 0, 0, 0, 0, 0, 1, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
