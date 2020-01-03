% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPPP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:38
	% EndTime: 2019-12-31 19:25:38
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
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
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
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (16->12), mult. (65->33), div. (0->0), fcn. (96->8), ass. (0->21)
	t68 = sin(pkin(5));
	t72 = sin(qJ(1));
	t86 = t68 * t72;
	t70 = cos(pkin(5));
	t71 = sin(qJ(2));
	t85 = t70 * t71;
	t84 = t71 * t68;
	t83 = t72 * t70;
	t67 = sin(pkin(8));
	t73 = cos(qJ(2));
	t82 = t73 * t67;
	t69 = cos(pkin(8));
	t81 = t73 * t69;
	t74 = cos(qJ(1));
	t80 = t73 * t74;
	t79 = t74 * t70;
	t78 = t67 * t71 - t70 * t81;
	t77 = -t69 * t71 - t70 * t82;
	t76 = -t71 * t79 + t86;
	t75 = t68 * t74 + t71 * t83;
	t1 = [t75 * t67 - t72 * t81, t77 * t74, 0, 0, 0; t76 * t67 + t69 * t80, t77 * t72, 0, 0, 0; 0, -t67 * t85 + t81, 0, 0, 0; t75 * t69 + t72 * t82, t78 * t74, 0, 0, 0; -t67 * t80 + t76 * t69, t78 * t72, 0, 0, 0; 0, -t69 * t85 - t82, 0, 0, 0; -t72 * t84 + t79, t68 * t80, 0, 0, 0; t74 * t84 + t83, t73 * t86, 0, 0, 0; 0, t84, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (16->12), mult. (65->33), div. (0->0), fcn. (96->8), ass. (0->21)
	t85 = sin(pkin(5));
	t89 = sin(qJ(1));
	t103 = t85 * t89;
	t87 = cos(pkin(5));
	t88 = sin(qJ(2));
	t102 = t87 * t88;
	t101 = t88 * t85;
	t100 = t89 * t87;
	t84 = sin(pkin(8));
	t90 = cos(qJ(2));
	t99 = t90 * t84;
	t86 = cos(pkin(8));
	t98 = t90 * t86;
	t91 = cos(qJ(1));
	t97 = t90 * t91;
	t96 = t91 * t87;
	t95 = -t84 * t88 + t87 * t98;
	t94 = t86 * t88 + t87 * t99;
	t93 = t88 * t96 - t103;
	t92 = -t88 * t100 - t85 * t91;
	t1 = [-t89 * t101 + t96, t85 * t97, 0, 0, 0; t91 * t101 + t100, t90 * t103, 0, 0, 0; 0, t101, 0, 0, 0; t92 * t84 + t89 * t98, t94 * t91, 0, 0, 0; t93 * t84 - t86 * t97, t94 * t89, 0, 0, 0; 0, t84 * t102 - t98, 0, 0, 0; t92 * t86 - t89 * t99, t95 * t91, 0, 0, 0; t84 * t97 + t93 * t86, t95 * t89, 0, 0, 0; 0, t86 * t102 + t99, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->11), mult. (65->33), div. (0->0), fcn. (96->8), ass. (0->21)
	t95 = sin(pkin(5));
	t99 = sin(qJ(1));
	t113 = t95 * t99;
	t97 = cos(pkin(5));
	t98 = sin(qJ(2));
	t112 = t97 * t98;
	t111 = t98 * t95;
	t110 = t99 * t97;
	t100 = cos(qJ(2));
	t109 = t100 * t97;
	t108 = t100 * t99;
	t101 = cos(qJ(1));
	t107 = t101 * t95;
	t94 = sin(pkin(8));
	t96 = cos(pkin(8));
	t106 = t100 * t94 + t96 * t112;
	t105 = t100 * t96 - t94 * t112;
	t104 = t96 * t109 - t94 * t98;
	t103 = -t94 * t109 - t96 * t98;
	t102 = t98 * t110 + t107;
	t1 = [t101 * t97 - t99 * t111, t100 * t107, 0, 0, 0; t98 * t107 + t110, t95 * t108, 0, 0, 0; 0, t111, 0, 0, 0; -t102 * t96 - t94 * t108, t104 * t101, 0, 0, 0; t106 * t101 - t96 * t113, t104 * t99, 0, 0, 0; 0, t106, 0, 0, 0; t102 * t94 - t96 * t108, t103 * t101, 0, 0, 0; t105 * t101 + t94 * t113, t103 * t99, 0, 0, 0; 0, t105, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end