% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:07
% EndTime: 2019-07-18 13:30:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (64->18), mult. (16->10), div. (0->0), fcn. (36->8), ass. (0->19)
t12 = qJ(2) + qJ(3);
t22 = pkin(1) + 0;
t14 = sin(qJ(2));
t21 = t14 * pkin(2) + 0;
t11 = qJ(1) + 0;
t6 = sin(t12);
t20 = pkin(3) * t6 + t21;
t16 = cos(qJ(2));
t19 = t16 * pkin(2) + t22;
t18 = pkin(4) + t11;
t7 = cos(t12);
t17 = pkin(3) * t7 + t19;
t15 = cos(qJ(5));
t13 = sin(qJ(5));
t8 = qJ(4) + t12;
t5 = pkin(5) + t18;
t4 = cos(t8);
t3 = sin(t8);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t16, -t14, 0, t22; t14, t16, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t7, -t6, 0, t19; t6, t7, 0, t21; 0, 0, 1, t18; 0, 0, 0, 1; t4, -t3, 0, t17; t3, t4, 0, t20; 0, 0, 1, t5; 0, 0, 0, 1; t4 * t15, -t4 * t13, t3, t3 * pkin(6) + t17; t3 * t15, -t3 * t13, -t4, -t4 * pkin(6) + t20; t13, t15, 0, t5; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
