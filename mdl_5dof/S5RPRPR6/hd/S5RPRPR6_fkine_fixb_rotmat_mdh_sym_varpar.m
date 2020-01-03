% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:37
% EndTime: 2019-12-31 18:17:37
% DurationCPUTime: 0.07s
% Computational Cost: add. (104->25), mult. (28->14), div. (0->0), fcn. (52->8), ass. (0->20)
t26 = pkin(5) + 0;
t14 = qJ(1) + pkin(8);
t16 = sin(qJ(1));
t25 = t16 * pkin(1) + 0;
t18 = cos(qJ(1));
t24 = t18 * pkin(1) + 0;
t9 = sin(t14);
t23 = pkin(2) * t9 + t25;
t10 = cos(t14);
t22 = pkin(2) * t10 + t24;
t21 = qJ(2) + t26;
t8 = pkin(6) + t21;
t11 = qJ(3) + t14;
t6 = sin(t11);
t7 = cos(t11);
t20 = t7 * pkin(3) + t6 * qJ(4) + t22;
t19 = t6 * pkin(3) - t7 * qJ(4) + t23;
t17 = cos(qJ(5));
t15 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t16, 0, 0; t16, t18, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t10, -t9, 0, t24; t9, t10, 0, t25; 0, 0, 1, t21; 0, 0, 0, 1; t7, -t6, 0, t22; t6, t7, 0, t23; 0, 0, 1, t8; 0, 0, 0, 1; 0, -t7, t6, t20; 0, -t6, -t7, t19; 1, 0, 0, t8; 0, 0, 0, 1; t6 * t15, t6 * t17, t7, t7 * pkin(7) + t20; -t7 * t15, -t7 * t17, t6, t6 * pkin(7) + t19; t17, -t15, 0, pkin(4) + t8; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
