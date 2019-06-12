% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a4]';
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
% Datum: 2019-06-03 15:11
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5PRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-03 15:11:16
% EndTime: 2019-06-03 15:11:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (42->21), mult. (34->21), div. (0->0), fcn. (67->8), ass. (0->22)
t10 = sin(qJ(2));
t7 = qJ(3) + qJ(4);
t4 = sin(t7);
t21 = t10 * t4;
t8 = sin(qJ(5));
t20 = t10 * t8;
t13 = cos(qJ(2));
t19 = t13 * t4;
t18 = t13 * t8;
t11 = cos(qJ(5));
t17 = t10 * t11;
t12 = cos(qJ(3));
t16 = t10 * t12;
t15 = t13 * t11;
t14 = t13 * t12;
t6 = qJ(1) + 0;
t9 = sin(qJ(3));
t5 = cos(t7);
t3 = -t9 * pkin(1) + 0;
t2 = pkin(1) * t14 + 0;
t1 = pkin(1) * t16 + t6;
t22 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, t6; 0, 0, 0, 1; t13, -t10, 0, 0; 0, 0, -1, 0; t10, t13, 0, t6; 0, 0, 0, 1; t14, -t13 * t9, t10, 0; -t9, -t12, 0, 0; t16, -t10 * t9, -t13, t6; 0, 0, 0, 1; t13 * t5, -t19, t10, t2; -t4, -t5, 0, t3; t10 * t5, -t21, -t13, t1; 0, 0, 0, 1; t5 * t15 + t20, -t5 * t18 + t17, t19, t2; -t4 * t11, t4 * t8, t5, t3; t5 * t17 - t18, -t5 * t20 - t15, t21, t1; 0, 0, 0, 1;];
T_ges = t22;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
