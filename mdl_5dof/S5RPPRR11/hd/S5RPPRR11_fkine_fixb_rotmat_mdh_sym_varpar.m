% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-29 16:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 16:25:43
% EndTime: 2019-12-29 16:25:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (56->30), mult. (56->28), div. (0->0), fcn. (89->6), ass. (0->23)
t10 = sin(qJ(1));
t9 = sin(qJ(4));
t28 = t10 * t9;
t13 = cos(qJ(1));
t27 = t13 * t9;
t11 = cos(qJ(5));
t26 = t10 * t11;
t12 = cos(qJ(4));
t25 = t10 * t12;
t24 = t13 * t11;
t23 = t13 * t12;
t7 = pkin(5) + 0;
t22 = pkin(2) + t7;
t21 = t13 * pkin(1) + t10 * qJ(2) + 0;
t20 = pkin(3) + t22;
t19 = t13 * qJ(3) + t21;
t18 = pkin(4) * t9 - pkin(7) * t12;
t17 = t10 * pkin(1) - t13 * qJ(2) + 0;
t16 = t10 * qJ(3) + t17;
t15 = -t10 * pkin(6) + t19;
t14 = t13 * pkin(6) + t16;
t8 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t13, -t10, 0, 0; t10, t13, 0, 0; 0, 0, 1, t7; 0, 0, 0, 1; 0, -t13, t10, t21; 0, -t10, -t13, t17; 1, 0, 0, t7; 0, 0, 0, 1; 0, t10, t13, t19; 0, -t13, t10, t16; 1, 0, 0, t22; 0, 0, 0, 1; t27, t23, -t10, t15; t28, t25, t13, t14; t12, -t9, 0, t20; 0, 0, 0, 1; -t10 * t8 + t9 * t24, -t8 * t27 - t26, -t23, t18 * t13 + t15; t13 * t8 + t9 * t26, -t8 * t28 + t24, -t25, t18 * t10 + t14; t12 * t11, -t12 * t8, t9, t12 * pkin(4) + t9 * pkin(7) + t20; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
