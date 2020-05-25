% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:01
% EndTime: 2018-11-23 15:37:01
% DurationCPUTime: 0.09s
% Computational Cost: add. (131->35), mult. (66->32), div. (0->0), fcn. (103->8), ass. (0->27)
t17 = cos(qJ(5));
t12 = qJ(1) + pkin(9);
t7 = sin(t12);
t34 = t7 * t17;
t8 = cos(t12);
t33 = t8 * t17;
t13 = sin(qJ(6));
t14 = sin(qJ(5));
t32 = t13 * t14;
t16 = cos(qJ(6));
t31 = t14 * t16;
t30 = pkin(6) + 0;
t15 = sin(qJ(1));
t29 = t15 * pkin(1) + 0;
t18 = cos(qJ(1));
t28 = t18 * pkin(1) + 0;
t9 = qJ(2) + t30;
t27 = t8 * pkin(2) + t7 * qJ(3) + t28;
t26 = pkin(3) + t9;
t25 = pkin(5) * t14 - pkin(8) * t17;
t24 = t8 * qJ(4) + t27;
t23 = pkin(4) + t26;
t22 = t7 * pkin(2) - t8 * qJ(3) + t29;
t21 = t7 * qJ(4) + t22;
t20 = -t7 * pkin(7) + t24;
t19 = t8 * pkin(7) + t21;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t15, 0, 0; t15, t18, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t8, -t7, 0, t28; t7, t8, 0, t29; 0, 0, 1, t9; 0, 0, 0, 1; 0, -t8, t7, t27; 0, -t7, -t8, t22; 1, 0, 0, t9; 0, 0, 0, 1; 0, t7, t8, t24; 0, -t8, t7, t21; 1, 0, 0, t26; 0, 0, 0, 1; t8 * t14, t33, -t7, t20; t7 * t14, t34, t8, t19; t17, -t14, 0, t23; 0, 0, 0, 1; -t7 * t13 + t8 * t31, -t7 * t16 - t8 * t32, -t33, t25 * t8 + t20; t8 * t13 + t7 * t31, t8 * t16 - t7 * t32, -t34, t25 * t7 + t19; t17 * t16, -t17 * t13, t14, t17 * pkin(5) + t14 * pkin(8) + t23; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
