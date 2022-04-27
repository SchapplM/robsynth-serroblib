% Calculate homogenous joint transformation matrices for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-31 23:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PPRRR4_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-31 23:19:33
% EndTime: 2022-01-31 23:19:33
% DurationCPUTime: 0.03s
% Computational Cost: add. (11->11), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t106 = cos(qJ(3));
t105 = cos(qJ(4));
t104 = cos(qJ(5));
t103 = sin(qJ(3));
t102 = sin(qJ(4));
t101 = sin(qJ(5));
t100 = cos(pkin(5));
t99 = cos(pkin(6));
t98 = cos(pkin(10));
t97 = cos(pkin(11));
t96 = sin(pkin(5));
t95 = sin(pkin(6));
t94 = sin(pkin(10));
t93 = sin(pkin(11));
t1 = [t98, -t94, 0, 0; t94, t98, 0, 0; 0, 0, 1, qJ(1); t97, -t93, 0, pkin(1); t100 * t93, t100 * t97, -t96, -t96 * qJ(2); t96 * t93, t96 * t97, t100, t100 * qJ(2); t106, -t103, 0, pkin(2); t99 * t103, t99 * t106, -t95, -t95 * pkin(7); t95 * t103, t95 * t106, t99, t99 * pkin(7); t105, -t102, 0, pkin(3); 0, 0, -1, -pkin(8); t102, t105, 0, 0; t104, -t101, 0, pkin(4); 0, 0, -1, -pkin(9); t101, t104, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
