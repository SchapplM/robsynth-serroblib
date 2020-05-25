% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:31
% EndTime: 2019-12-31 21:59:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (110->50), mult. (117->54), div. (0->0), fcn. (168->8), ass. (0->34)
t26 = -pkin(8) - pkin(7);
t20 = sin(qJ(3));
t37 = t20 * pkin(3);
t23 = cos(qJ(3));
t10 = t23 * pkin(3) + pkin(2);
t19 = qJ(3) + qJ(4);
t11 = sin(t19);
t21 = sin(qJ(2));
t36 = t21 * t11;
t22 = sin(qJ(1));
t35 = t22 * t20;
t24 = cos(qJ(2));
t34 = t22 * t24;
t25 = cos(qJ(1));
t33 = t25 * t24;
t18 = pkin(5) + 0;
t32 = t22 * pkin(1) + 0;
t31 = t25 * pkin(1) + t22 * pkin(6) + 0;
t30 = pkin(2) * t24 + pkin(7) * t21;
t17 = -qJ(5) + t26;
t12 = cos(t19);
t5 = pkin(4) * t12 + t10;
t29 = -t17 * t21 + t24 * t5;
t28 = t10 * t24 - t21 * t26;
t27 = -t25 * pkin(6) + t32;
t9 = t25 * t21;
t8 = t22 * t21;
t7 = t21 * t12;
t6 = pkin(4) * t11 + t37;
t4 = t22 * t11 + t12 * t33;
t3 = -t11 * t33 + t22 * t12;
t2 = -t25 * t11 + t12 * t34;
t1 = -t11 * t34 - t25 * t12;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t22, 0, 0; t22, t25, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t33, -t9, t22, t31; t34, -t8, -t25, t27; t21, t24, 0, t18; 0, 0, 0, 1; t23 * t33 + t35, -t20 * t33 + t22 * t23, t9, t30 * t25 + t31; -t25 * t20 + t23 * t34, -t20 * t34 - t25 * t23, t8, t30 * t22 + t27; t21 * t23, -t21 * t20, -t24, t21 * pkin(2) - t24 * pkin(7) + t18; 0, 0, 0, 1; t4, t3, t9, pkin(3) * t35 + t28 * t25 + t31; t2, t1, t8, (-pkin(6) - t37) * t25 + t28 * t22 + t32; t7, -t36, -t24, t21 * t10 + t24 * t26 + t18; 0, 0, 0, 1; t4, t3, t9, t22 * t6 + t29 * t25 + t31; t2, t1, t8, (-pkin(6) - t6) * t25 + t29 * t22 + t32; t7, -t36, -t24, t24 * t17 + t21 * t5 + t18; 0, 0, 0, 1;];
T_ges = t13;
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
