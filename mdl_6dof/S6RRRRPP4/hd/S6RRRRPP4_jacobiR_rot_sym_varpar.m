% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.05s
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
	t1 = [t56, -t60 * t66, t57, 0, 0, 0; t58, -t60 * t70, -t55, 0, 0, 0; 0, t68, -t60 * t59, 0, 0, 0; t55, t59 * t67, -t58, 0, 0, 0; t57, t59 * t71, t56, 0, 0, 0; 0, -t69, -t60 * t62, 0, 0, 0; -t71, t65, 0, 0, 0, 0; t67, t61 * t63, 0, 0, 0, 0; 0, t60, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (54->17), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t80 = qJ(3) + qJ(4);
	t78 = sin(t80);
	t81 = sin(qJ(2));
	t91 = t81 * t78;
	t79 = cos(t80);
	t90 = t81 * t79;
	t82 = sin(qJ(1));
	t89 = t82 * t81;
	t83 = cos(qJ(2));
	t88 = t83 * t78;
	t87 = t83 * t79;
	t84 = cos(qJ(1));
	t86 = t84 * t81;
	t85 = t84 * t83;
	t77 = t82 * t78 + t79 * t85;
	t76 = -t78 * t85 + t82 * t79;
	t75 = t84 * t78 - t82 * t87;
	t74 = t84 * t79 + t82 * t88;
	t1 = [t75, -t79 * t86, t76, t76, 0, 0; t77, -t79 * t89, -t74, -t74, 0, 0; 0, t87, -t91, -t91, 0, 0; t74, t78 * t86, -t77, -t77, 0, 0; t76, t78 * t89, t75, t75, 0, 0; 0, -t88, -t90, -t90, 0, 0; -t89, t85, 0, 0, 0, 0; t86, t82 * t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (88->17), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t87 = qJ(3) + qJ(4) + pkin(10);
	t85 = sin(t87);
	t88 = sin(qJ(2));
	t98 = t88 * t85;
	t86 = cos(t87);
	t97 = t88 * t86;
	t89 = sin(qJ(1));
	t96 = t89 * t88;
	t90 = cos(qJ(2));
	t95 = t90 * t85;
	t94 = t90 * t86;
	t91 = cos(qJ(1));
	t93 = t91 * t88;
	t92 = t91 * t90;
	t84 = t89 * t85 + t86 * t92;
	t83 = -t85 * t92 + t89 * t86;
	t82 = t91 * t85 - t89 * t94;
	t81 = t91 * t86 + t89 * t95;
	t1 = [t82, -t86 * t93, t83, t83, 0, 0; t84, -t86 * t96, -t81, -t81, 0, 0; 0, t94, -t98, -t98, 0, 0; t81, t85 * t93, -t84, -t84, 0, 0; t83, t85 * t96, t82, t82, 0, 0; 0, -t95, -t97, -t97, 0, 0; -t96, t92, 0, 0, 0, 0; t93, t89 * t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (87->15), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t111 = qJ(3) + qJ(4) + pkin(10);
	t109 = sin(t111);
	t112 = sin(qJ(2));
	t121 = t112 * t109;
	t113 = sin(qJ(1));
	t120 = t113 * t112;
	t114 = cos(qJ(2));
	t119 = t114 * t109;
	t110 = cos(t111);
	t118 = t114 * t110;
	t115 = cos(qJ(1));
	t117 = t115 * t112;
	t116 = t115 * t114;
	t107 = t112 * t110;
	t106 = t113 * t109 + t110 * t116;
	t105 = t109 * t116 - t113 * t110;
	t104 = -t115 * t109 + t113 * t118;
	t103 = -t115 * t110 - t113 * t119;
	t1 = [-t104, -t110 * t117, -t105, -t105, 0, 0; t106, -t110 * t120, t103, t103, 0, 0; 0, t118, -t121, -t121, 0, 0; -t120, t116, 0, 0, 0, 0; t117, t113 * t114, 0, 0, 0, 0; 0, t112, 0, 0, 0, 0; t103, -t109 * t117, t106, t106, 0, 0; t105, -t109 * t120, t104, t104, 0, 0; 0, t119, t107, t107, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end