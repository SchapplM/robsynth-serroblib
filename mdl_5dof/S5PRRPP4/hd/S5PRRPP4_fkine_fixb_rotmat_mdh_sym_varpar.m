% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRPP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:41
% EndTime: 2019-12-31 17:40:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->31), mult. (59->22), div. (0->0), fcn. (91->6), ass. (0->23)
t18 = pkin(7) + qJ(2);
t12 = sin(t18);
t21 = sin(qJ(3));
t31 = qJ(4) * t21;
t22 = cos(qJ(3));
t6 = t12 * t22;
t32 = pkin(3) * t6 + t12 * t31;
t13 = cos(t18);
t8 = t13 * t22;
t19 = sin(pkin(7));
t30 = t19 * pkin(1) + 0;
t20 = cos(pkin(7));
t29 = t20 * pkin(1) + 0;
t28 = qJ(1) + 0;
t14 = pkin(5) + t28;
t27 = t12 * pkin(2) + t30;
t26 = t13 * pkin(2) + t12 * pkin(6) + t29;
t25 = pkin(3) * t8 + t13 * t31 + t26;
t24 = -t13 * pkin(6) + t27;
t23 = t21 * pkin(3) - t22 * qJ(4) + t14;
t7 = t13 * t21;
t5 = t12 * t21;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t19, 0, 0; t19, t20, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t13, -t12, 0, t29; t12, t13, 0, t30; 0, 0, 1, t14; 0, 0, 0, 1; t8, -t7, t12, t26; t6, -t5, -t13, t24; t21, t22, 0, t14; 0, 0, 0, 1; t8, t12, t7, t25; t6, -t13, t5, t24 + t32; t21, 0, -t22, t23; 0, 0, 0, 1; t8, t7, -t12, pkin(4) * t8 - t12 * qJ(5) + t25; t6, t5, t13, pkin(4) * t6 + (-pkin(6) + qJ(5)) * t13 + t27 + t32; t21, -t22, 0, t21 * pkin(4) + t23; 0, 0, 0, 1;];
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
