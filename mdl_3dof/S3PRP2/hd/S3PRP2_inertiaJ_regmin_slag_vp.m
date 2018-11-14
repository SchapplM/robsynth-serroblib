% Calculate minimal parameter regressor of joint inertia matrix for
% S3PRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% 
% Output:
% MM_reg [((3+1)*3/2)x7]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MM_reg = S3PRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_inertiaJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t2 = cos(qJ(2));
t1 = sin(qJ(2));
t3 = [1, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2; 0, 0, t2, -t1, t2, t1, pkin(2) * t2 + qJ(3) * t1; 0, 1, 0, 0, 0.2e1 * pkin(2), 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2; 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, -1, 0, -pkin(2); 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
