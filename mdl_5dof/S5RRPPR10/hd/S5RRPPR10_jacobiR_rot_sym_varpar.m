% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPPR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
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
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t47 = sin(qJ(2));
	t48 = sin(qJ(1));
	t55 = t48 * t47;
	t45 = sin(pkin(8));
	t49 = cos(qJ(2));
	t54 = t49 * t45;
	t46 = cos(pkin(8));
	t53 = t49 * t46;
	t50 = cos(qJ(1));
	t52 = t50 * t47;
	t51 = t50 * t49;
	t1 = [t50 * t45 - t48 * t53, -t46 * t52, 0, 0, 0; t48 * t45 + t46 * t51, -t46 * t55, 0, 0, 0; 0, t53, 0, 0, 0; t50 * t46 + t48 * t54, t45 * t52, 0, 0, 0; -t45 * t51 + t48 * t46, t45 * t55, 0, 0, 0; 0, -t54, 0, 0, 0; -t55, t51, 0, 0, 0; t52, t48 * t49, 0, 0, 0; 0, t47, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t62 = sin(qJ(2));
	t63 = sin(qJ(1));
	t70 = t63 * t62;
	t60 = sin(pkin(8));
	t64 = cos(qJ(2));
	t69 = t64 * t60;
	t61 = cos(pkin(8));
	t68 = t64 * t61;
	t65 = cos(qJ(1));
	t67 = t65 * t62;
	t66 = t65 * t64;
	t1 = [t65 * t60 - t63 * t68, -t61 * t67, 0, 0, 0; t63 * t60 + t61 * t66, -t61 * t70, 0, 0, 0; 0, t68, 0, 0, 0; -t70, t66, 0, 0, 0; t67, t63 * t64, 0, 0, 0; 0, t62, 0, 0, 0; -t65 * t61 - t63 * t69, -t60 * t67, 0, 0, 0; t60 * t66 - t63 * t61, -t60 * t70, 0, 0, 0; 0, t69, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->21), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
	t79 = sin(qJ(2));
	t76 = sin(pkin(8));
	t77 = cos(pkin(8));
	t78 = sin(qJ(5));
	t81 = cos(qJ(5));
	t87 = t76 * t81 - t77 * t78;
	t85 = t87 * t79;
	t80 = sin(qJ(1));
	t82 = cos(qJ(2));
	t91 = t80 * t82;
	t83 = cos(qJ(1));
	t90 = t83 * t82;
	t71 = t76 * t91 + t83 * t77;
	t72 = -t83 * t76 + t77 * t91;
	t89 = t71 * t81 - t72 * t78;
	t88 = -t71 * t78 - t72 * t81;
	t86 = t76 * t78 + t77 * t81;
	t84 = t86 * t79;
	t74 = t80 * t76 + t77 * t90;
	t73 = t76 * t90 - t80 * t77;
	t70 = t73 * t78 + t74 * t81;
	t69 = t73 * t81 - t74 * t78;
	t1 = [t88, -t83 * t84, 0, 0, t69; t70, -t80 * t84, 0, 0, t89; 0, t86 * t82, 0, 0, t85; -t89, -t83 * t85, 0, 0, -t70; t69, -t80 * t85, 0, 0, t88; 0, t87 * t82, 0, 0, -t84; t80 * t79, -t90, 0, 0, 0; -t83 * t79, -t91, 0, 0, 0; 0, -t79, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end