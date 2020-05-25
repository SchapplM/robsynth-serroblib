% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RRPP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:48
% EndTime: 2019-12-31 16:58:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (41->26), mult. (51->20), div. (0->0), fcn. (79->4), ass. (0->17)
t15 = sin(qJ(1));
t14 = sin(qJ(2));
t23 = qJ(3) * t14;
t16 = cos(qJ(2));
t6 = t15 * t16;
t24 = pkin(2) * t6 + t15 * t23;
t17 = cos(qJ(1));
t8 = t17 * t16;
t13 = pkin(4) + 0;
t22 = t15 * pkin(1) + 0;
t21 = t17 * pkin(1) + t15 * pkin(5) + 0;
t20 = -t17 * pkin(5) + t22;
t19 = pkin(2) * t8 + t17 * t23 + t21;
t18 = t14 * pkin(2) - t16 * qJ(3) + t13;
t7 = t17 * t14;
t5 = t15 * t14;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t15, 0, 0; t15, t17, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; t8, -t7, t15, t21; t6, -t5, -t17, t20; t14, t16, 0, t13; 0, 0, 0, 1; t8, t15, t7, t19; t6, -t17, t5, t20 + t24; t14, 0, -t16, t18; 0, 0, 0, 1; t8, t7, -t15, pkin(3) * t8 - t15 * qJ(4) + t19; t6, t5, t17, pkin(3) * t6 + (-pkin(5) + qJ(4)) * t17 + t22 + t24; t14, -t16, 0, t14 * pkin(3) + t18; 0, 0, 0, 1;];
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
