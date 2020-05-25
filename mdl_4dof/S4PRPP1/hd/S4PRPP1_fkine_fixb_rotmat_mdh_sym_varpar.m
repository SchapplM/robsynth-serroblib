% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:41
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PRPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:57
% EndTime: 2018-11-14 13:40:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->18), mult. (16->8), div. (0->0), fcn. (32->4), ass. (0->12)
t10 = sin(pkin(5));
t16 = t10 * pkin(1) + 0;
t11 = cos(pkin(5));
t15 = t11 * pkin(1) + 0;
t14 = qJ(1) + 0;
t6 = pkin(4) + t14;
t9 = pkin(5) + qJ(2);
t4 = sin(t9);
t5 = cos(t9);
t13 = t5 * pkin(2) + t4 * qJ(3) + t15;
t12 = t4 * pkin(2) - t5 * qJ(3) + t16;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t11, -t10, 0, 0; t10, t11, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t5, -t4, 0, t15; t4, t5, 0, t16; 0, 0, 1, t6; 0, 0, 0, 1; 0, -t5, t4, t13; 0, -t4, -t5, t12; 1, 0, 0, t6; 0, 0, 0, 1; 0, t4, t5, t5 * qJ(4) + t13; 0, -t5, t4, t4 * qJ(4) + t12; 1, 0, 0, pkin(3) + t6; 0, 0, 0, 1;];
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
