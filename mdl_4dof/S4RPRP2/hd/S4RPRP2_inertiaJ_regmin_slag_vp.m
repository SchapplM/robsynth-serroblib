% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t5 = sin(qJ(3));
t6 = cos(qJ(3));
t7 = -pkin(1) - pkin(2);
t2 = t5 * qJ(2) - t6 * t7;
t3 = t6 * qJ(2) + t5 * t7;
t1 = -pkin(3) - t2;
t4 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 1, 0.2e1 * t2, 0.2e1 * t3, t1 ^ 2 + t3 ^ 2; 0, 0, 0, -1, 0, -pkin(1), 0, -t6, t5, t1 * t6 + t3 * t5; 0, 0, 0, 0, 0, 1, 0, 0, 0, t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, -1, -t2, -t3, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, t6, -t5, t6 * pkin(3); 0, 0, 0, 0, 0, 0, 1, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
