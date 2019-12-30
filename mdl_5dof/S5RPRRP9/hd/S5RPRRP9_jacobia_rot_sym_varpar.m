% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP9
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
%   Wie in S5RPRRP9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRP9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:17
	% EndTime: 2019-12-29 17:27:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:17
	% EndTime: 2019-12-29 17:27:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:27:12
	% EndTime: 2019-12-29 17:27:12
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (471->18), mult. (227->41), div. (75->9), fcn. (351->7), ass. (0->29)
	t55 = cos(qJ(1));
	t56 = t55 ^ 2;
	t50 = pkin(8) + qJ(3) + qJ(4);
	t49 = cos(t50);
	t48 = sin(t50);
	t54 = sin(qJ(1));
	t58 = t54 * t48;
	t42 = atan2(-t58, -t49);
	t40 = sin(t42);
	t41 = cos(t42);
	t38 = -t40 * t58 - t41 * t49;
	t37 = 0.1e1 / t38 ^ 2;
	t63 = t37 * t48;
	t62 = t40 * t49;
	t45 = t48 ^ 2;
	t47 = 0.1e1 / t49 ^ 2;
	t61 = t45 * t47;
	t51 = t54 ^ 2;
	t60 = t51 / t56;
	t43 = 0.1e1 / (t51 * t61 + 0.1e1);
	t59 = t54 * t43;
	t44 = 0.1e1 / (t47 * t60 + 0.1e1);
	t57 = 0.1e1 / t55 * t47 * t44 * t58;
	t46 = 0.1e1 / t49;
	t39 = (0.1e1 + t61) * t59;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t56 * t45 * t37 + 0.1e1);
	t34 = (t49 * t36 - (-t54 * t62 + t41 * t48 + (-t41 * t58 + t62) * t39) * t63) * t55 * t35;
	t1 = [t55 * t48 * t46 * t43, 0, t39, t39, 0; (-t36 * t58 - (-t41 * t45 * t46 * t59 + (t43 - 0.1e1) * t48 * t40) * t56 * t63) * t35, 0, t34, t34, 0; (-0.1e1 - t60) * t46 * t44, 0, -t57, -t57, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end