% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 13:06:46
% EndTime: 2019-12-29 13:06:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (70->25), mult. (44->28), div. (0->0), fcn. (73->8), ass. (0->21)
t11 = sin(qJ(3));
t9 = qJ(1) + pkin(7);
t4 = sin(t9);
t25 = t4 * t11;
t5 = cos(t9);
t24 = t5 * t11;
t10 = sin(qJ(4));
t14 = cos(qJ(3));
t23 = t10 * t14;
t13 = cos(qJ(4));
t22 = t13 * t14;
t21 = pkin(4) + 0;
t12 = sin(qJ(1));
t20 = t12 * pkin(1) + 0;
t15 = cos(qJ(1));
t19 = t15 * pkin(1) + 0;
t6 = qJ(2) + t21;
t18 = t5 * pkin(2) + t4 * pkin(5) + t19;
t17 = pkin(3) * t14 + pkin(6) * t11;
t16 = t4 * pkin(2) - t5 * pkin(5) + t20;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t12, 0, 0; t12, t15, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t5, -t4, 0, t19; t4, t5, 0, t20; 0, 0, 1, t6; 0, 0, 0, 1; t5 * t14, -t24, t4, t18; t4 * t14, -t25, -t5, t16; t11, t14, 0, t6; 0, 0, 0, 1; t4 * t10 + t5 * t22, t4 * t13 - t5 * t23, t24, t17 * t5 + t18; -t5 * t10 + t4 * t22, -t5 * t13 - t4 * t23, t25, t17 * t4 + t16; t11 * t13, -t11 * t10, -t14, t11 * pkin(3) - t14 * pkin(6) + t6; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
