% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:35
% EndTime: 2019-12-31 18:23:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (112->38), mult. (71->34), div. (0->0), fcn. (108->8), ass. (0->27)
t16 = qJ(1) + pkin(8);
t10 = sin(t16);
t18 = sin(qJ(3));
t31 = qJ(4) * t18;
t21 = cos(qJ(3));
t5 = t10 * t21;
t36 = pkin(3) * t5 + t10 * t31;
t35 = t10 * t18;
t11 = cos(t16);
t34 = t11 * t18;
t6 = t11 * t21;
t17 = sin(qJ(5));
t33 = t17 * t18;
t20 = cos(qJ(5));
t32 = t18 * t20;
t30 = pkin(5) + 0;
t19 = sin(qJ(1));
t29 = t19 * pkin(1) + 0;
t22 = cos(qJ(1));
t28 = t22 * pkin(1) + 0;
t27 = t10 * pkin(2) + t29;
t12 = qJ(2) + t30;
t26 = t11 * pkin(2) + t10 * pkin(6) + t28;
t25 = pkin(3) * t6 + t11 * t31 + t26;
t24 = -t11 * pkin(6) + t27;
t23 = t18 * pkin(3) - t21 * qJ(4) + t12;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t19, 0, 0; t19, t22, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t11, -t10, 0, t28; t10, t11, 0, t29; 0, 0, 1, t12; 0, 0, 0, 1; t6, -t34, t10, t26; t5, -t35, -t11, t24; t18, t21, 0, t12; 0, 0, 0, 1; t10, -t6, t34, t25; -t11, -t5, t35, t24 + t36; 0, -t18, -t21, t23; 0, 0, 0, 1; t10 * t20 + t11 * t33, -t10 * t17 + t11 * t32, t6, t10 * pkin(4) + pkin(7) * t6 + t25; t10 * t33 - t11 * t20, t10 * t32 + t11 * t17, t5, pkin(7) * t5 + (-pkin(4) - pkin(6)) * t11 + t27 + t36; -t21 * t17, -t21 * t20, t18, t18 * pkin(7) + t23; 0, 0, 0, 1;];
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
