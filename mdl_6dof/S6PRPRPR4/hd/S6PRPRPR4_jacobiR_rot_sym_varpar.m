% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (27->18), div. (0->0), fcn. (42->8), ass. (0->14)
	t55 = sin(pkin(6));
	t60 = cos(qJ(2));
	t63 = t55 * t60;
	t58 = cos(pkin(6));
	t59 = sin(qJ(2));
	t62 = t58 * t59;
	t61 = t58 * t60;
	t57 = cos(pkin(10));
	t56 = cos(pkin(11));
	t54 = sin(pkin(10));
	t53 = sin(pkin(11));
	t52 = -t54 * t61 - t57 * t59;
	t51 = -t54 * t59 + t57 * t61;
	t1 = [0, t52 * t56, 0, 0, 0, 0; 0, t51 * t56, 0, 0, 0, 0; 0, t56 * t63, 0, 0, 0, 0; 0, -t52 * t53, 0, 0, 0, 0; 0, -t51 * t53, 0, 0, 0, 0; 0, -t53 * t63, 0, 0, 0, 0; 0, -t54 * t62 + t57 * t60, 0, 0, 0, 0; 0, t54 * t60 + t57 * t62, 0, 0, 0, 0; 0, t55 * t59, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t72 = sin(pkin(10));
	t73 = sin(pkin(6));
	t83 = t72 * t73;
	t74 = cos(pkin(10));
	t82 = t73 * t74;
	t76 = sin(qJ(2));
	t81 = t73 * t76;
	t77 = cos(qJ(2));
	t80 = t73 * t77;
	t75 = cos(pkin(6));
	t79 = t75 * t76;
	t78 = t75 * t77;
	t71 = pkin(11) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = -t72 * t79 + t74 * t77;
	t67 = -t72 * t78 - t74 * t76;
	t66 = t72 * t77 + t74 * t79;
	t65 = -t72 * t76 + t74 * t78;
	t1 = [0, t67 * t70, 0, -t68 * t69 + t70 * t83, 0, 0; 0, t65 * t70, 0, -t66 * t69 - t70 * t82, 0, 0; 0, t70 * t80, 0, -t69 * t81 + t75 * t70, 0, 0; 0, -t67 * t69, 0, -t68 * t70 - t69 * t83, 0, 0; 0, -t65 * t69, 0, -t66 * t70 + t69 * t82, 0, 0; 0, -t69 * t80, 0, -t75 * t69 - t70 * t81, 0, 0; 0, t68, 0, 0, 0, 0; 0, t66, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (61->20), mult. (107->52), div. (0->0), fcn. (156->10), ass. (0->27)
	t105 = pkin(11) + qJ(4);
	t104 = cos(t105);
	t106 = sin(pkin(12));
	t121 = t104 * t106;
	t109 = cos(pkin(12));
	t120 = t104 * t109;
	t113 = cos(qJ(2));
	t119 = t104 * t113;
	t107 = sin(pkin(10));
	t108 = sin(pkin(6));
	t118 = t107 * t108;
	t110 = cos(pkin(10));
	t117 = t108 * t110;
	t112 = sin(qJ(2));
	t116 = t108 * t112;
	t111 = cos(pkin(6));
	t115 = t111 * t112;
	t114 = t111 * t113;
	t103 = sin(t105);
	t102 = -t107 * t115 + t110 * t113;
	t101 = -t107 * t114 - t110 * t112;
	t100 = t107 * t113 + t110 * t115;
	t99 = -t107 * t112 + t110 * t114;
	t98 = -t103 * t116 + t111 * t104;
	t97 = -t102 * t103 + t104 * t118;
	t96 = -t100 * t103 - t104 * t117;
	t1 = [0, t101 * t120 + t102 * t106, 0, t97 * t109, 0, 0; 0, t100 * t106 + t99 * t120, 0, t96 * t109, 0, 0; 0, (t106 * t112 + t109 * t119) * t108, 0, t98 * t109, 0, 0; 0, -t101 * t121 + t102 * t109, 0, -t97 * t106, 0, 0; 0, t100 * t109 - t99 * t121, 0, -t96 * t106, 0, 0; 0, (-t106 * t119 + t109 * t112) * t108, 0, -t98 * t106, 0, 0; 0, t101 * t103, 0, t102 * t104 + t103 * t118, 0, 0; 0, t99 * t103, 0, t100 * t104 - t103 * t117, 0, 0; 0, t108 * t113 * t103, 0, t111 * t103 + t104 * t116, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (123->29), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->32)
	t123 = pkin(12) + qJ(6);
	t119 = sin(t123);
	t124 = pkin(11) + qJ(4);
	t122 = cos(t124);
	t139 = t119 * t122;
	t121 = cos(t123);
	t138 = t121 * t122;
	t130 = cos(qJ(2));
	t137 = t122 * t130;
	t125 = sin(pkin(10));
	t126 = sin(pkin(6));
	t136 = t125 * t126;
	t127 = cos(pkin(10));
	t135 = t126 * t127;
	t129 = sin(qJ(2));
	t134 = t126 * t129;
	t133 = t126 * t130;
	t128 = cos(pkin(6));
	t132 = t128 * t129;
	t131 = t128 * t130;
	t120 = sin(t124);
	t117 = -t125 * t132 + t127 * t130;
	t116 = t125 * t131 + t127 * t129;
	t115 = t125 * t130 + t127 * t132;
	t114 = t125 * t129 - t127 * t131;
	t113 = t128 * t120 + t122 * t134;
	t112 = -t120 * t134 + t128 * t122;
	t111 = t117 * t122 + t120 * t136;
	t110 = -t117 * t120 + t122 * t136;
	t109 = t115 * t122 - t120 * t135;
	t108 = -t115 * t120 - t122 * t135;
	t1 = [0, -t116 * t138 + t117 * t119, 0, t110 * t121, 0, -t111 * t119 + t116 * t121; 0, -t114 * t138 + t115 * t119, 0, t108 * t121, 0, -t109 * t119 + t114 * t121; 0, (t119 * t129 + t121 * t137) * t126, 0, t112 * t121, 0, -t113 * t119 - t121 * t133; 0, t116 * t139 + t117 * t121, 0, -t110 * t119, 0, -t111 * t121 - t116 * t119; 0, t114 * t139 + t115 * t121, 0, -t108 * t119, 0, -t109 * t121 - t114 * t119; 0, (-t119 * t137 + t121 * t129) * t126, 0, -t112 * t119, 0, -t113 * t121 + t119 * t133; 0, -t116 * t120, 0, t111, 0, 0; 0, -t114 * t120, 0, t109, 0, 0; 0, t120 * t133, 0, t113, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end