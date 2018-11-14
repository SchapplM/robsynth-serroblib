% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:01
% EndTime: 2018-11-14 13:39:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->22), mult. (44->14), div. (0->0), fcn. (72->4), ass. (0->13)
t20 = cos(qJ(3));
t19 = sin(qJ(3));
t12 = qJ(1) + 0;
t13 = sin(pkin(5));
t14 = cos(pkin(5));
t18 = t14 * pkin(1) + t13 * qJ(2) + 0;
t17 = t14 * pkin(2) + t18;
t16 = t13 * pkin(1) - t14 * qJ(2) + 0;
t15 = t13 * pkin(2) + t16;
t7 = -pkin(4) + t12;
t2 = -t13 * t20 + t14 * t19;
t1 = -t13 * t19 - t14 * t20;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t14, -t13, 0, 0; t13, t14, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; t14, 0, t13, t18; t13, 0, -t14, t16; 0, 1, 0, t12; 0, 0, 0, 1; -t1, -t2, 0, t17; -t2, t1, 0, t15; 0, 0, -1, t7; 0, 0, 0, 1; -t1, 0, t2, -t1 * pkin(3) + t2 * qJ(4) + t17; -t2, 0, -t1, -t2 * pkin(3) - t1 * qJ(4) + t15; 0, -1, 0, t7; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
