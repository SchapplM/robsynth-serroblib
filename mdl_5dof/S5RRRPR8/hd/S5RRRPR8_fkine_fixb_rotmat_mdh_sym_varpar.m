% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:18:59
% EndTime: 2019-12-31 21:18:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->44), mult. (80->40), div. (0->0), fcn. (121->8), ass. (0->29)
t22 = cos(qJ(1));
t16 = qJ(2) + qJ(3);
t12 = sin(t16);
t30 = qJ(4) * t12;
t13 = cos(t16);
t8 = t22 * t13;
t37 = pkin(3) * t8 + t22 * t30;
t19 = sin(qJ(1));
t36 = t19 * t12;
t7 = t19 * t13;
t17 = sin(qJ(5));
t35 = t19 * t17;
t20 = cos(qJ(5));
t34 = t19 * t20;
t33 = t22 * t12;
t32 = t22 * t17;
t31 = t22 * t20;
t15 = pkin(5) + 0;
t21 = cos(qJ(2));
t10 = t21 * pkin(2) + pkin(1);
t29 = t22 * t10 + 0;
t18 = sin(qJ(2));
t28 = t18 * pkin(2) + t15;
t23 = -pkin(7) - pkin(6);
t27 = t19 * t10 + t22 * t23 + 0;
t26 = pkin(3) * t7 + t19 * t30 + t27;
t25 = -t19 * t23 + t29;
t24 = t12 * pkin(3) - t13 * qJ(4) + t28;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t22 * t21, -t22 * t18, t19, t22 * pkin(1) + t19 * pkin(6) + 0; t19 * t21, -t19 * t18, -t22, t19 * pkin(1) - t22 * pkin(6) + 0; t18, t21, 0, t15; 0, 0, 0, 1; t8, -t33, t19, t25; t7, -t36, -t22, t27; t12, t13, 0, t28; 0, 0, 0, 1; t19, -t8, t33, t25 + t37; -t22, -t7, t36, t26; 0, -t12, -t13, t24; 0, 0, 0, 1; t12 * t32 + t34, t12 * t31 - t35, t8, pkin(8) * t8 + (pkin(4) - t23) * t19 + t29 + t37; t12 * t35 - t31, t12 * t34 + t32, t7, -t22 * pkin(4) + pkin(8) * t7 + t26; -t13 * t17, -t13 * t20, t12, t12 * pkin(8) + t24; 0, 0, 0, 1;];
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
