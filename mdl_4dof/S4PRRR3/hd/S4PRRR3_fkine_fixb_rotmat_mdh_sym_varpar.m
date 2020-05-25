% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR3
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:31
% EndTime: 2019-12-31 16:31:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (65->19), mult. (18->12), div. (0->0), fcn. (38->8), ass. (0->18)
t12 = sin(pkin(7));
t21 = t12 * pkin(1) + 0;
t11 = pkin(7) + qJ(2);
t13 = cos(pkin(7));
t20 = t13 * pkin(1) + 0;
t19 = qJ(1) + 0;
t6 = sin(t11);
t18 = pkin(2) * t6 + t21;
t7 = cos(t11);
t17 = pkin(2) * t7 + t20;
t16 = pkin(4) + t19;
t15 = cos(qJ(4));
t14 = sin(qJ(4));
t8 = qJ(3) + t11;
t5 = pkin(5) + t16;
t4 = cos(t8);
t3 = sin(t8);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t13, -t12, 0, 0; t12, t13, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t7, -t6, 0, t20; t6, t7, 0, t21; 0, 0, 1, t16; 0, 0, 0, 1; t4, -t3, 0, t17; t3, t4, 0, t18; 0, 0, 1, t5; 0, 0, 0, 1; t4 * t15, -t4 * t14, t3, t4 * pkin(3) + t3 * pkin(6) + t17; t3 * t15, -t3 * t14, -t4, t3 * pkin(3) - t4 * pkin(6) + t18; t14, t15, 0, t5; 0, 0, 0, 1;];
T_ges = t1;
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
