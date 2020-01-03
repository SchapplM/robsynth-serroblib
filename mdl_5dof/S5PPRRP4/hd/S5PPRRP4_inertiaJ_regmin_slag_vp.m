% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t9 = cos(qJ(4));
t16 = 0.2e1 * t9;
t7 = sin(qJ(4));
t8 = sin(qJ(3));
t15 = t7 * t8;
t14 = t9 * pkin(4);
t4 = t7 ^ 2;
t13 = t9 ^ 2 + t4;
t12 = -qJ(5) - pkin(6);
t1 = t12 * t7;
t2 = t12 * t9;
t11 = -t1 * t7 - t2 * t9;
t10 = cos(qJ(3));
t3 = -pkin(3) - t14;
t5 = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t8 ^ 2 + t10 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t1 + t7 * t2; 0, 0, 0, t10, -t8, 0, 0, 0, 0, 0, t10 * t9, -t10 * t7, t13 * t8, -t10 * t3 + t11 * t8; 0, 0, 1, 0, 0, t4, t7 * t16, 0, 0, 0, pkin(3) * t16, -0.2e1 * pkin(3) * t7, 0.2e1 * t11, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t7, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t9 * t8, 0, -pkin(4) * t15; 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, -t7 * pkin(6), -t9 * pkin(6), -t7 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
