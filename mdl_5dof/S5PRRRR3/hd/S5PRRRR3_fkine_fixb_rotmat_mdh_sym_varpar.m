% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-10-24 10:36
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:36:10
% EndTime: 2019-10-24 10:36:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (111->24), mult. (26->14), div. (0->0), fcn. (50->10), ass. (0->24)
t16 = pkin(9) + qJ(2);
t17 = sin(pkin(9));
t29 = t17 * pkin(1) + 0;
t18 = cos(pkin(9));
t28 = t18 * pkin(1) + 0;
t27 = qJ(1) + 0;
t11 = sin(t16);
t26 = pkin(2) * t11 + t29;
t12 = cos(t16);
t25 = pkin(2) * t12 + t28;
t24 = pkin(5) + t27;
t13 = qJ(3) + t16;
t8 = sin(t13);
t23 = pkin(3) * t8 + t26;
t9 = cos(t13);
t22 = pkin(3) * t9 + t25;
t21 = pkin(6) + t24;
t20 = cos(qJ(5));
t19 = sin(qJ(5));
t10 = qJ(4) + t13;
t5 = pkin(7) + t21;
t4 = cos(t10);
t3 = sin(t10);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t17, 0, 0; t17, t18, 0, 0; 0, 0, 1, t27; 0, 0, 0, 1; t12, -t11, 0, t28; t11, t12, 0, t29; 0, 0, 1, t24; 0, 0, 0, 1; t9, -t8, 0, t25; t8, t9, 0, t26; 0, 0, 1, t21; 0, 0, 0, 1; t4, -t3, 0, t22; t3, t4, 0, t23; 0, 0, 1, t5; 0, 0, 0, 1; t4 * t20, -t4 * t19, t3, t4 * pkin(4) + t3 * pkin(8) + t22; t3 * t20, -t3 * t19, -t4, t3 * pkin(4) - t4 * pkin(8) + t23; t19, t20, 0, t5; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
