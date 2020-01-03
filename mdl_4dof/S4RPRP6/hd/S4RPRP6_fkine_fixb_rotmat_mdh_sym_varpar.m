% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRP6
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:55
% EndTime: 2019-12-31 16:45:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (37->25), mult. (29->16), div. (0->0), fcn. (53->4), ass. (0->15)
t8 = sin(qJ(3));
t9 = sin(qJ(1));
t1 = t9 * t8;
t11 = cos(qJ(1));
t17 = t11 * t8;
t10 = cos(qJ(3));
t16 = t11 * t10;
t6 = pkin(4) + 0;
t15 = t9 * pkin(1) + 0;
t14 = pkin(2) + t6;
t13 = t11 * pkin(1) + t9 * qJ(2) + 0;
t12 = -t11 * qJ(2) + t15;
t7 = -qJ(4) - pkin(5);
t2 = t9 * t10;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t11, -t9, 0, 0; t9, t11, 0, 0; 0, 0, 1, t6; 0, 0, 0, 1; 0, -t11, t9, t13; 0, -t9, -t11, t12; 1, 0, 0, t6; 0, 0, 0, 1; t1, t2, t11, t11 * pkin(5) + t13; -t17, -t16, t9, t9 * pkin(5) + t12; t10, -t8, 0, t14; 0, 0, 0, 1; t1, t2, t11, pkin(3) * t1 - t11 * t7 + t13; -t17, -t16, t9, -t9 * t7 + (-pkin(3) * t8 - qJ(2)) * t11 + t15; t10, -t8, 0, t10 * pkin(3) + t14; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
