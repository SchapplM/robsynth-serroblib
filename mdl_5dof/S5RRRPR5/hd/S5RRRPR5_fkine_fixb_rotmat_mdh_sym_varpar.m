% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:12:54
% EndTime: 2019-12-31 21:12:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (122->42), mult. (69->44), div. (0->0), fcn. (110->10), ass. (0->29)
t23 = -pkin(7) - pkin(6);
t19 = sin(qJ(1));
t16 = qJ(2) + qJ(3);
t9 = pkin(9) + t16;
t5 = sin(t9);
t34 = t19 * t5;
t22 = cos(qJ(1));
t33 = t22 * t5;
t21 = cos(qJ(2));
t8 = t21 * pkin(2) + pkin(1);
t17 = sin(qJ(5));
t32 = t19 * t17;
t20 = cos(qJ(5));
t31 = t19 * t20;
t30 = t22 * t17;
t29 = t22 * t20;
t15 = pkin(5) + 0;
t14 = -qJ(4) + t23;
t11 = cos(t16);
t3 = pkin(3) * t11 + t8;
t28 = t22 * t14 + t19 * t3 + 0;
t18 = sin(qJ(2));
t27 = t18 * pkin(2) + t15;
t10 = sin(t16);
t26 = pkin(3) * t10 + t27;
t6 = cos(t9);
t25 = pkin(4) * t6 + pkin(8) * t5;
t24 = -t19 * t14 + t22 * t3 + 0;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t22 * t21, -t22 * t18, t19, t22 * pkin(1) + t19 * pkin(6) + 0; t19 * t21, -t19 * t18, -t22, t19 * pkin(1) - t22 * pkin(6) + 0; t18, t21, 0, t15; 0, 0, 0, 1; t22 * t11, -t22 * t10, t19, -t19 * t23 + t22 * t8 + 0; t19 * t11, -t19 * t10, -t22, t19 * t8 + t22 * t23 + 0; t10, t11, 0, t27; 0, 0, 0, 1; t22 * t6, -t33, t19, t24; t19 * t6, -t34, -t22, t28; t5, t6, 0, t26; 0, 0, 0, 1; t6 * t29 + t32, -t6 * t30 + t31, t33, t25 * t22 + t24; t6 * t31 - t30, -t6 * t32 - t29, t34, t25 * t19 + t28; t5 * t20, -t5 * t17, -t6, t5 * pkin(4) - t6 * pkin(8) + t26; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
