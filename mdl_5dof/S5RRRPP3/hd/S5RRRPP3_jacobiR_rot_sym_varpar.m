% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0; t15, t15, 0, 0, 0; 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0; -t14, -t14, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:06
	% EndTime: 2019-12-29 19:40:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t26 = qJ(1) + qJ(2);
	t24 = sin(t26);
	t28 = cos(qJ(3));
	t30 = t24 * t28;
	t25 = cos(t26);
	t27 = sin(qJ(3));
	t29 = t25 * t27;
	t23 = t25 * t28;
	t22 = t24 * t27;
	t1 = [-t30, -t30, -t29, 0, 0; t23, t23, -t22, 0, 0; 0, 0, t28, 0, 0; t22, t22, -t23, 0, 0; -t29, -t29, -t30, 0, 0; 0, 0, -t27, 0, 0; t25, t25, 0, 0, 0; t24, t24, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:06
	% EndTime: 2019-12-29 19:40:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (21->6), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t63 = qJ(1) + qJ(2);
	t61 = sin(t63);
	t64 = sin(qJ(3));
	t67 = t61 * t64;
	t62 = cos(t63);
	t65 = cos(qJ(3));
	t66 = t62 * t65;
	t60 = t62 * t64;
	t59 = t61 * t65;
	t1 = [t62, t62, 0, 0, 0; t61, t61, 0, 0, 0; 0, 0, 0, 0, 0; t59, t59, t60, 0, 0; -t66, -t66, t67, 0, 0; 0, 0, -t65, 0, 0; -t67, -t67, t66, 0, 0; t60, t60, t59, 0, 0; 0, 0, t64, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->7), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t64 = qJ(1) + qJ(2);
	t62 = sin(t64);
	t65 = sin(qJ(3));
	t68 = t62 * t65;
	t66 = cos(qJ(3));
	t67 = t62 * t66;
	t63 = cos(t64);
	t61 = t63 * t66;
	t60 = t63 * t65;
	t1 = [t63, t63, 0, 0, 0; t62, t62, 0, 0, 0; 0, 0, 0, 0, 0; -t68, -t68, t61, 0, 0; t60, t60, t67, 0, 0; 0, 0, t65, 0, 0; -t67, -t67, -t60, 0, 0; t61, t61, -t68, 0, 0; 0, 0, t66, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end