% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_inertiaJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t7 = sin(qJ(4));
t13 = t7 * pkin(2);
t12 = sin(qJ(3)) * pkin(1);
t9 = cos(qJ(4));
t11 = t9 * t12;
t6 = cos(qJ(3)) * pkin(1);
t4 = t6 + pkin(2);
t1 = -t7 * t12 + t9 * t4;
t5 = t9 * pkin(2);
t2 = -t7 * t4 - t11;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t6, -0.2e1 * t12, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t6, -t12, 1, t1 + t5, -t11 + (-pkin(2) - t4) * t7; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t5, -0.2e1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, t5, -t13; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
