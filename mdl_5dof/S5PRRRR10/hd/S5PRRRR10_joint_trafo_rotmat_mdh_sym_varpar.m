% Calculate homogenous joint transformation matrices for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-01 09:46
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PRRRR10_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-01 09:45:56
% EndTime: 2022-02-01 09:45:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (11->11), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t180 = cos(qJ(2));
t179 = cos(qJ(3));
t178 = cos(qJ(4));
t177 = cos(qJ(5));
t176 = sin(qJ(2));
t175 = sin(qJ(3));
t174 = sin(qJ(4));
t173 = sin(qJ(5));
t172 = cos(pkin(5));
t171 = cos(pkin(6));
t170 = cos(pkin(11));
t169 = sin(pkin(5));
t168 = sin(pkin(6));
t167 = sin(pkin(11));
t1 = [t170, -t167, 0, 0; t167, t170, 0, 0; 0, 0, 1, qJ(1); t180, -t176, 0, pkin(1); t172 * t176, t172 * t180, -t169, -t169 * pkin(7); t169 * t176, t169 * t180, t172, t172 * pkin(7); t179, -t175, 0, pkin(2); t171 * t175, t171 * t179, -t168, -t168 * pkin(8); t168 * t175, t168 * t179, t171, t171 * pkin(8); t178, -t174, 0, pkin(3); 0, 0, -1, -pkin(9); t174, t178, 0, 0; t177, -t173, 0, pkin(4); 0, 0, -1, -pkin(10); t173, t177, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
