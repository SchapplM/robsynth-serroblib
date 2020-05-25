% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RPRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:58
% EndTime: 2019-12-31 16:46:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (37->24), mult. (34->16), div. (0->0), fcn. (58->4), ass. (0->17)
t11 = cos(qJ(1));
t8 = sin(qJ(3));
t20 = t11 * t8;
t10 = cos(qJ(3));
t9 = sin(qJ(1));
t19 = t9 * t10;
t18 = t11 * t10;
t7 = pkin(4) + 0;
t17 = t9 * pkin(1) + 0;
t16 = pkin(2) + t7;
t15 = t11 * pkin(1) + t9 * qJ(2) + 0;
t14 = t11 * pkin(5) + t15;
t13 = pkin(3) * t8 - qJ(4) * t10;
t12 = -t11 * qJ(2) + t17;
t3 = t9 * pkin(5);
t1 = t9 * t8;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t11, -t9, 0, 0; t9, t11, 0, 0; 0, 0, 1, t7; 0, 0, 0, 1; 0, -t11, t9, t15; 0, -t9, -t11, t12; 1, 0, 0, t7; 0, 0, 0, 1; t1, t19, t11, t14; -t20, -t18, t9, t12 + t3; t10, -t8, 0, t16; 0, 0, 0, 1; t1, t11, -t19, t13 * t9 + t14; -t20, t9, t18, t3 + (-qJ(2) - t13) * t11 + t17; t10, 0, t8, t10 * pkin(3) + t8 * qJ(4) + t16; 0, 0, 0, 1;];
T_ges = t2;
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
