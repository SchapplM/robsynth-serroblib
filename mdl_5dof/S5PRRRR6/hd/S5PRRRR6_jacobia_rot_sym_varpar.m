% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:37
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:20
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (257->15), mult. (243->34), div. (59->8), fcn. (349->9), ass. (0->26)
	t55 = qJ(2) + qJ(3);
	t54 = cos(t55);
	t53 = sin(t55);
	t56 = sin(pkin(9));
	t61 = t56 * t53;
	t49 = atan2(-t61, -t54);
	t47 = sin(t49);
	t64 = t47 * t54;
	t51 = t53 ^ 2;
	t63 = t51 / t54 ^ 2;
	t57 = cos(pkin(9));
	t62 = t54 * t57;
	t58 = sin(qJ(4));
	t59 = cos(qJ(4));
	t46 = t56 * t58 + t59 * t62;
	t44 = 0.1e1 / t46 ^ 2;
	t45 = -t56 * t59 + t58 * t62;
	t60 = t45 ^ 2 * t44 + 0.1e1;
	t48 = cos(t49);
	t43 = 0.1e1 / t60;
	t42 = (0.1e1 + t63) * t56 / (t56 ^ 2 * t63 + 0.1e1);
	t41 = -t47 * t61 - t48 * t54;
	t40 = 0.1e1 / t41 ^ 2;
	t38 = (-t58 / t46 + t59 * t45 * t44) * t57 * t53 * t43;
	t37 = (t54 / t41 - (-t56 * t64 + t48 * t53 + (-t48 * t61 + t64) * t42) * t53 * t40) * t57 / (t57 ^ 2 * t51 * t40 + 0.1e1);
	t1 = [0, t42, t42, 0, 0; 0, t37, t37, 0, 0; 0, t38, t38, t60 * t43, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (334->16), mult. (270->34), div. (64->8), fcn. (384->9), ass. (0->28)
	t62 = qJ(2) + qJ(3);
	t60 = cos(t62);
	t58 = sin(t62);
	t63 = sin(pkin(9));
	t66 = t63 * t58;
	t53 = atan2(-t66, -t60);
	t51 = sin(t53);
	t69 = t51 * t60;
	t55 = t58 ^ 2;
	t68 = t55 / t60 ^ 2;
	t64 = cos(pkin(9));
	t67 = t60 * t64;
	t61 = qJ(4) + qJ(5);
	t57 = sin(t61);
	t59 = cos(t61);
	t50 = t63 * t57 + t59 * t67;
	t48 = 0.1e1 / t50 ^ 2;
	t49 = t57 * t67 - t63 * t59;
	t65 = t49 ^ 2 * t48 + 0.1e1;
	t52 = cos(t53);
	t47 = (0.1e1 + t68) * t63 / (t63 ^ 2 * t68 + 0.1e1);
	t46 = 0.1e1 / t65;
	t45 = -t51 * t66 - t52 * t60;
	t44 = 0.1e1 / t45 ^ 2;
	t42 = t65 * t46;
	t41 = (-t57 / t50 + t59 * t49 * t48) * t64 * t58 * t46;
	t40 = (t60 / t45 - (-t63 * t69 + t52 * t58 + (-t52 * t66 + t69) * t47) * t58 * t44) * t64 / (t64 ^ 2 * t55 * t44 + 0.1e1);
	t1 = [0, t47, t47, 0, 0; 0, t40, t40, 0, 0; 0, t41, t41, t42, t42;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end