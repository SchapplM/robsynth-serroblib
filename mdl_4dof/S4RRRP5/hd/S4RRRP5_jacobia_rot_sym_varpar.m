% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RRRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:59:41
	% EndTime: 2019-12-29 13:59:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:59:41
	% EndTime: 2019-12-29 13:59:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:59:41
	% EndTime: 2019-12-29 13:59:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:59:41
	% EndTime: 2019-12-29 13:59:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:59:41
	% EndTime: 2019-12-29 13:59:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (299->18), mult. (227->41), div. (75->9), fcn. (351->7), ass. (0->29)
	t54 = cos(qJ(1));
	t55 = t54 ^ 2;
	t52 = qJ(2) + qJ(3);
	t48 = cos(t52);
	t47 = sin(t52);
	t53 = sin(qJ(1));
	t57 = t53 * t47;
	t41 = atan2(-t57, -t48);
	t39 = sin(t41);
	t40 = cos(t41);
	t37 = -t39 * t57 - t40 * t48;
	t36 = 0.1e1 / t37 ^ 2;
	t62 = t36 * t47;
	t61 = t39 * t48;
	t44 = t47 ^ 2;
	t46 = 0.1e1 / t48 ^ 2;
	t60 = t44 * t46;
	t49 = t53 ^ 2;
	t59 = t49 / t55;
	t42 = 0.1e1 / (t49 * t60 + 0.1e1);
	t58 = t53 * t42;
	t43 = 0.1e1 / (t46 * t59 + 0.1e1);
	t56 = 0.1e1 / t54 * t46 * t43 * t57;
	t45 = 0.1e1 / t48;
	t38 = (0.1e1 + t60) * t58;
	t35 = 0.1e1 / t37;
	t34 = 0.1e1 / (t55 * t44 * t36 + 0.1e1);
	t33 = (t48 * t35 - (-t53 * t61 + t40 * t47 + (-t40 * t57 + t61) * t38) * t62) * t54 * t34;
	t1 = [t54 * t47 * t45 * t42, t38, t38, 0; (-t35 * t57 - (-t40 * t44 * t45 * t58 + (t42 - 0.1e1) * t47 * t39) * t55 * t62) * t34, t33, t33, 0; (-0.1e1 - t59) * t45 * t43, -t56, -t56, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end