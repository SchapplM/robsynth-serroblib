% Calculate homogenous joint transformation matrices for
% S5RRRRP2
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
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RRRRP2_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:36
% EndTime: 2022-01-20 11:48:36
% DurationCPUTime: 0.02s
% Computational Cost: add. (5->5), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t35 = cos(qJ(1));
t34 = cos(qJ(2));
t33 = cos(qJ(3));
t32 = cos(qJ(4));
t31 = sin(qJ(1));
t30 = sin(qJ(2));
t29 = sin(qJ(3));
t28 = sin(qJ(4));
t1 = [t35, -t31, 0, 0; t31, t35, 0, 0; 0, 0, 1, pkin(5); t34, -t30, 0, pkin(1); t30, t34, 0, 0; 0, 0, 1, pkin(6); t33, -t29, 0, pkin(2); 0, 0, -1, -pkin(7); t29, t33, 0, 0; t32, -t28, 0, pkin(3); t28, t32, 0, 0; 0, 0, 1, pkin(8); 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, qJ(5);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
