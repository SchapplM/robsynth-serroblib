% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2018-11-14 13:56
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PPPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:56:21
% EndTime: 2018-11-14 13:56:21
% DurationCPUTime: 0.05s
% Computational Cost: add. (23->11), mult. (2->2), div. (0->0), fcn. (10->4), ass. (0->11)
t10 = -pkin(1) + 0;
t4 = qJ(1) + 0;
t3 = qJ(2) + 0;
t9 = pkin(2) + t4;
t8 = -qJ(3) + t10;
t7 = cos(pkin(5));
t6 = sin(pkin(5));
t5 = pkin(5) + qJ(4);
t2 = cos(t5);
t1 = sin(t5);
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, -1, 0, 0; 0, 0, 1, t4; -1, 0, 0, 0; 0, 0, 0, 1; 0, 0, 1, t3; 1, 0, 0, t4; 0, 1, 0, t10; 0, 0, 0, 1; t6, t7, 0, t3; t7, -t6, 0, t9; 0, 0, -1, t8; 0, 0, 0, 1; t1, t2, 0, t6 * pkin(3) + t3; t2, -t1, 0, t7 * pkin(3) + t9; 0, 0, -1, -pkin(4) + t8; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
