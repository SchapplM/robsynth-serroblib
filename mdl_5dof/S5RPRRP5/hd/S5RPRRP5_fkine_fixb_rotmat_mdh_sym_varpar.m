% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:35
% EndTime: 2019-12-31 18:40:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (114->26), mult. (40->18), div. (0->0), fcn. (68->8), ass. (0->25)
t17 = sin(qJ(4));
t16 = qJ(1) + pkin(8);
t13 = qJ(3) + t16;
t8 = sin(t13);
t31 = t8 * t17;
t9 = cos(t13);
t30 = t9 * t17;
t29 = pkin(5) + 0;
t18 = sin(qJ(1));
t28 = t18 * pkin(1) + 0;
t20 = cos(qJ(1));
t27 = t20 * pkin(1) + 0;
t11 = sin(t16);
t26 = pkin(2) * t11 + t28;
t12 = cos(t16);
t25 = pkin(2) * t12 + t27;
t24 = qJ(2) + t29;
t10 = pkin(6) + t24;
t23 = t9 * pkin(3) + t8 * pkin(7) + t25;
t19 = cos(qJ(4));
t22 = pkin(4) * t19 + qJ(5) * t17;
t21 = t8 * pkin(3) - t9 * pkin(7) + t26;
t2 = t9 * t19;
t1 = t8 * t19;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t18, 0, 0; t18, t20, 0, 0; 0, 0, 1, t29; 0, 0, 0, 1; t12, -t11, 0, t27; t11, t12, 0, t28; 0, 0, 1, t24; 0, 0, 0, 1; t9, -t8, 0, t25; t8, t9, 0, t26; 0, 0, 1, t10; 0, 0, 0, 1; t2, -t30, t8, t23; t1, -t31, -t9, t21; t17, t19, 0, t10; 0, 0, 0, 1; t2, t8, t30, t22 * t9 + t23; t1, -t9, t31, t22 * t8 + t21; t17, 0, -t19, t17 * pkin(4) - t19 * qJ(5) + t10; 0, 0, 0, 1;];
T_ges = t3;
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
