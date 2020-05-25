% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:02
% EndTime: 2019-12-31 18:47:02
% DurationCPUTime: 0.10s
% Computational Cost: add. (79->35), mult. (104->24), div. (0->0), fcn. (160->6), ass. (0->22)
t17 = sin(qJ(4));
t19 = cos(qJ(1));
t27 = sin(qJ(3));
t28 = sin(qJ(1));
t29 = cos(qJ(3));
t5 = -t19 * t29 - t27 * t28;
t31 = t5 * t17;
t6 = t19 * t27 - t28 * t29;
t30 = t6 * t17;
t16 = pkin(5) + 0;
t11 = -pkin(6) + t16;
t26 = t19 * pkin(1) + t28 * qJ(2) + 0;
t25 = t19 * pkin(2) + t26;
t18 = cos(qJ(4));
t24 = -pkin(4) * t18 - qJ(5) * t17;
t23 = t28 * pkin(1) - t19 * qJ(2) + 0;
t22 = t28 * pkin(2) + t23;
t21 = -t5 * pkin(3) + t6 * pkin(7) + t25;
t20 = -t6 * pkin(3) - t5 * pkin(7) + t22;
t2 = t5 * t18;
t1 = t6 * t18;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t28, 0, 0; t28, t19, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; t19, 0, t28, t26; t28, 0, -t19, t23; 0, 1, 0, t16; 0, 0, 0, 1; -t5, -t6, 0, t25; -t6, t5, 0, t22; 0, 0, -1, t11; 0, 0, 0, 1; -t2, t31, t6, t21; -t1, t30, -t5, t20; -t17, -t18, 0, t11; 0, 0, 0, 1; -t2, t6, -t31, t24 * t5 + t21; -t1, -t5, -t30, t24 * t6 + t20; -t17, 0, t18, -t17 * pkin(4) + t18 * qJ(5) + t11; 0, 0, 0, 1;];
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
