% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S4PRRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(7));
	t6 = sin(pkin(7));
	t1 = [0, -t7 * t8, 0, 0; 0, -t6 * t8, 0, 0; 0, t9, 0, 0; 0, -t7 * t9, 0, 0; 0, -t6 * t9, 0, 0; 0, -t8, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:56
	% EndTime: 2019-12-29 12:29:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->11), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t15 = qJ(2) + qJ(3);
	t13 = sin(t15);
	t16 = sin(pkin(7));
	t21 = t16 * t13;
	t14 = cos(t15);
	t20 = t16 * t14;
	t17 = cos(pkin(7));
	t19 = t17 * t13;
	t18 = t17 * t14;
	t1 = [0, -t19, -t19, 0; 0, -t21, -t21, 0; 0, t14, t14, 0; 0, -t18, -t18, 0; 0, -t20, -t20, 0; 0, -t13, -t13, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:56
	% EndTime: 2019-12-29 12:29:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (36->13), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t76 = qJ(2) + qJ(3);
	t75 = cos(t76);
	t79 = sin(qJ(4));
	t87 = t75 * t79;
	t77 = sin(pkin(7));
	t86 = t77 * t79;
	t80 = cos(qJ(4));
	t85 = t77 * t80;
	t78 = cos(pkin(7));
	t84 = t78 * t79;
	t83 = t78 * t80;
	t74 = sin(t76);
	t82 = t74 * t85;
	t81 = t74 * t83;
	t73 = t75 * t80;
	t72 = t78 * t75;
	t71 = t77 * t75;
	t70 = t74 * t84;
	t69 = t74 * t86;
	t1 = [0, -t81, -t81, -t75 * t84 + t85; 0, -t82, -t82, -t75 * t86 - t83; 0, t73, t73, -t74 * t79; 0, t70, t70, -t75 * t83 - t86; 0, t69, t69, -t75 * t85 + t84; 0, -t87, -t87, -t74 * t80; 0, t72, t72, 0; 0, t71, t71, 0; 0, t74, t74, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,4);
end