% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP8
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
%   Wie in S5RPRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->0), mult. (42->0), div. (11->0), fcn. (60->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; -1, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:47:46
	% EndTime: 2019-12-31 18:47:46
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (193->18), mult. (423->46), div. (73->9), fcn. (657->9), ass. (0->34)
	t48 = cos(qJ(4));
	t62 = sin(qJ(3));
	t63 = sin(qJ(1));
	t64 = cos(qJ(3));
	t65 = cos(qJ(1));
	t38 = t65 * t62 - t63 * t64;
	t47 = sin(qJ(4));
	t56 = t38 * t47;
	t32 = atan2(t56, t48);
	t30 = sin(t32);
	t31 = cos(t32);
	t28 = t30 * t56 + t31 * t48;
	t27 = 0.1e1 / t28 ^ 2;
	t44 = t47 ^ 2;
	t37 = -t63 * t62 - t65 * t64;
	t49 = t37 ^ 2;
	t24 = 0.1e1 / (t49 * t44 * t27 + 0.1e1);
	t26 = 0.1e1 / t28;
	t34 = t38 ^ 2;
	t46 = 0.1e1 / t48 ^ 2;
	t55 = t44 * t46;
	t33 = 0.1e1 / (t34 * t55 + 0.1e1);
	t45 = 0.1e1 / t48;
	t57 = t38 * t33;
	t61 = t27 * t47;
	t68 = (t49 * t61 * (-t31 * t44 * t45 * t57 + (t33 - 0.1e1) * t47 * t30) - t26 * t56) * t24;
	t58 = t34 / t49;
	t59 = t30 * t48;
	t29 = 0.1e1 / (t46 * t58 + 0.1e1);
	t54 = t45 * t29;
	t52 = t37 * t47 * t45 * t33;
	t35 = 0.1e1 / t37;
	t25 = (0.1e1 + t55) * t57;
	t1 = [t52, 0, -t52, t25, 0; -t68, 0, t68, (-t48 * t26 + (t38 * t59 - t31 * t47 + (t31 * t56 - t59) * t25) * t61) * t37 * t24, 0; (0.1e1 + t58) * t54, 0, (-t35 * t37 - t58) * t54, t35 * t46 * t29 * t56, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end