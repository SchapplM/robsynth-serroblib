% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPP8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0; t13, -t16, 0, 0, 0; 0, t11, 0, 0, 0; t16, -t13, 0, 0, 0; -t15, -t14, 0, 0, 0; 0, -t9, 0, 0, 0; t12, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t60 = sin(qJ(2));
	t61 = sin(qJ(1));
	t71 = t61 * t60;
	t62 = cos(qJ(3));
	t70 = t61 * t62;
	t59 = sin(qJ(3));
	t63 = cos(qJ(2));
	t69 = t63 * t59;
	t68 = t63 * t62;
	t64 = cos(qJ(1));
	t67 = t64 * t60;
	t66 = t64 * t62;
	t65 = t64 * t63;
	t58 = t61 * t59 + t62 * t65;
	t57 = -t59 * t65 + t70;
	t56 = t64 * t59 - t61 * t68;
	t55 = t61 * t69 + t66;
	t1 = [t56, -t60 * t66, t57, 0, 0; t58, -t60 * t70, -t55, 0, 0; 0, t68, -t60 * t59, 0, 0; t55, t59 * t67, -t58, 0, 0; t57, t59 * t71, t56, 0, 0; 0, -t69, -t60 * t62, 0, 0; -t71, t65, 0, 0, 0; t67, t61 * t63, 0, 0, 0; 0, t60, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->10), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t73 = sin(qJ(2));
	t74 = sin(qJ(1));
	t84 = t74 * t73;
	t75 = cos(qJ(3));
	t83 = t74 * t75;
	t72 = sin(qJ(3));
	t76 = cos(qJ(2));
	t82 = t76 * t72;
	t81 = t76 * t75;
	t77 = cos(qJ(1));
	t80 = t77 * t73;
	t79 = t77 * t75;
	t78 = t77 * t76;
	t71 = t74 * t72 + t75 * t78;
	t70 = t72 * t78 - t83;
	t69 = -t77 * t72 + t74 * t81;
	t68 = t74 * t82 + t79;
	t1 = [-t84, t78, 0, 0, 0; t80, t74 * t76, 0, 0, 0; 0, t73, 0, 0, 0; t69, t73 * t79, t70, 0, 0; -t71, t73 * t83, t68, 0, 0; 0, -t81, t73 * t72, 0, 0; -t68, -t72 * t80, t71, 0, 0; t70, -t72 * t84, t69, 0, 0; 0, t82, t73 * t75, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:10:08
	% EndTime: 2019-12-31 21:10:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t72 = sin(qJ(2));
	t73 = sin(qJ(1));
	t83 = t73 * t72;
	t74 = cos(qJ(3));
	t82 = t73 * t74;
	t71 = sin(qJ(3));
	t75 = cos(qJ(2));
	t81 = t75 * t71;
	t80 = t75 * t74;
	t76 = cos(qJ(1));
	t79 = t76 * t72;
	t78 = t76 * t74;
	t77 = t76 * t75;
	t70 = t73 * t71 + t74 * t77;
	t69 = t71 * t77 - t82;
	t68 = -t76 * t71 + t73 * t80;
	t67 = -t73 * t81 - t78;
	t1 = [-t83, t77, 0, 0, 0; t79, t73 * t75, 0, 0, 0; 0, t72, 0, 0, 0; t67, -t71 * t79, t70, 0, 0; t69, -t71 * t83, t68, 0, 0; 0, t81, t72 * t74, 0, 0; -t68, -t72 * t78, -t69, 0, 0; t70, -t72 * t82, t67, 0, 0; 0, t80, -t72 * t71, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end