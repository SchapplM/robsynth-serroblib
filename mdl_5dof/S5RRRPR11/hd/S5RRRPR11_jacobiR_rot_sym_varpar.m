% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:25
	% EndTime: 2019-12-29 20:16:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:14
	% EndTime: 2019-12-29 20:16:14
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
	% StartTime: 2019-12-29 20:16:19
	% EndTime: 2019-12-29 20:16:19
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-12-29 20:16:19
	% EndTime: 2019-12-29 20:16:19
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-12-29 20:16:26
	% EndTime: 2019-12-29 20:16:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t77 = sin(qJ(2));
	t78 = sin(qJ(1));
	t88 = t78 * t77;
	t79 = cos(qJ(3));
	t87 = t78 * t79;
	t76 = sin(qJ(3));
	t80 = cos(qJ(2));
	t86 = t80 * t76;
	t85 = t80 * t79;
	t81 = cos(qJ(1));
	t84 = t81 * t77;
	t83 = t81 * t79;
	t82 = t81 * t80;
	t75 = t78 * t76 + t79 * t82;
	t74 = t76 * t82 - t87;
	t73 = -t81 * t76 + t78 * t85;
	t72 = -t78 * t86 - t83;
	t1 = [-t73, -t77 * t83, -t74, 0, 0; t75, -t77 * t87, t72, 0, 0; 0, t85, -t77 * t76, 0, 0; -t88, t82, 0, 0, 0; t84, t78 * t80, 0, 0, 0; 0, t77, 0, 0, 0; t72, -t76 * t84, t75, 0, 0; t74, -t76 * t88, t73, 0, 0; 0, t86, t77 * t79, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:26
	% EndTime: 2019-12-29 20:16:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (50->24), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->23)
	t93 = sin(qJ(5));
	t94 = sin(qJ(3));
	t97 = cos(qJ(5));
	t98 = cos(qJ(3));
	t102 = t93 * t94 + t97 * t98;
	t95 = sin(qJ(2));
	t109 = t102 * t95;
	t100 = cos(qJ(1));
	t96 = sin(qJ(1));
	t99 = cos(qJ(2));
	t107 = t96 * t99;
	t87 = t100 * t98 + t94 * t107;
	t88 = -t100 * t94 + t98 * t107;
	t105 = -t87 * t97 + t88 * t93;
	t106 = t100 * t99;
	t104 = t87 * t93 + t88 * t97;
	t89 = t94 * t106 - t96 * t98;
	t90 = t98 * t106 + t96 * t94;
	t82 = t89 * t97 - t90 * t93;
	t83 = t89 * t93 + t90 * t97;
	t103 = t93 * t98 - t94 * t97;
	t85 = t103 * t95;
	t1 = [-t104, -t100 * t109, -t82, 0, t82; t83, -t96 * t109, t105, 0, -t105; 0, t102 * t99, t85, 0, -t85; t105, t100 * t85, t83, 0, -t83; t82, t96 * t85, t104, 0, -t104; 0, -t103 * t99, t109, 0, -t109; t96 * t95, -t106, 0, 0, 0; -t100 * t95, -t107, 0, 0, 0; 0, -t95, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end