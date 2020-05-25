% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:11
% EndTime: 2018-11-14 13:44:11
% DurationCPUTime: 0.06s
% Computational Cost: add. (58->17), mult. (12->6), div. (0->0), fcn. (28->8), ass. (0->19)
t13 = pkin(7) + qJ(2);
t14 = sin(pkin(7));
t22 = t14 * pkin(1) + 0;
t15 = cos(pkin(7));
t21 = t15 * pkin(1) + 0;
t20 = qJ(1) + 0;
t8 = sin(t13);
t19 = pkin(2) * t8 + t22;
t9 = cos(t13);
t18 = pkin(2) * t9 + t21;
t17 = pkin(4) + t20;
t10 = qJ(3) + t13;
t16 = pkin(5) + t17;
t7 = qJ(4) + t10;
t6 = cos(t10);
t5 = sin(t10);
t2 = cos(t7);
t1 = sin(t7);
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t14, 0, 0; t14, t15, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t9, -t8, 0, t21; t8, t9, 0, t22; 0, 0, 1, t17; 0, 0, 0, 1; t6, -t5, 0, t18; t5, t6, 0, t19; 0, 0, 1, t16; 0, 0, 0, 1; t2, -t1, 0, pkin(3) * t6 + t18; t1, t2, 0, pkin(3) * t5 + t19; 0, 0, 1, pkin(6) + t16; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
