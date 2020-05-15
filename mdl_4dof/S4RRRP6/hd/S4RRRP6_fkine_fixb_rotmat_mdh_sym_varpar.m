% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:04
% EndTime: 2019-12-31 17:18:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (52->32), mult. (75->34), div. (0->0), fcn. (113->6), ass. (0->26)
t14 = sin(qJ(3));
t15 = sin(qJ(2));
t28 = t15 * t14;
t16 = sin(qJ(1));
t27 = t16 * t14;
t18 = cos(qJ(2));
t26 = t16 * t18;
t19 = cos(qJ(1));
t25 = t19 * t18;
t12 = pkin(4) + 0;
t24 = t16 * pkin(1) + 0;
t23 = t19 * pkin(1) + t16 * pkin(5) + 0;
t22 = pkin(2) * t18 + pkin(6) * t15;
t13 = -qJ(4) - pkin(6);
t17 = cos(qJ(3));
t8 = t17 * pkin(3) + pkin(2);
t21 = -t13 * t15 + t18 * t8;
t20 = -t19 * pkin(5) + t24;
t7 = t19 * t15;
t6 = t16 * t15;
t5 = t15 * t17;
t4 = t17 * t25 + t27;
t3 = -t14 * t25 + t16 * t17;
t2 = -t19 * t14 + t17 * t26;
t1 = -t14 * t26 - t19 * t17;
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t16, 0, 0; t16, t19, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; t25, -t7, t16, t23; t26, -t6, -t19, t20; t15, t18, 0, t12; 0, 0, 0, 1; t4, t3, t7, t22 * t19 + t23; t2, t1, t6, t22 * t16 + t20; t5, -t28, -t18, t15 * pkin(2) - t18 * pkin(6) + t12; 0, 0, 0, 1; t4, t3, t7, pkin(3) * t27 + t21 * t19 + t23; t2, t1, t6, (-pkin(3) * t14 - pkin(5)) * t19 + t21 * t16 + t24; t5, -t28, -t18, t18 * t13 + t15 * t8 + t12; 0, 0, 0, 1;];
T_ges = t9;
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
