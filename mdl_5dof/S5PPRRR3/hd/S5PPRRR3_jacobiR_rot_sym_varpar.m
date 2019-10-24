% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:21
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:18
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:18
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t22 = sin(pkin(8));
	t25 = sin(qJ(3));
	t30 = t22 * t25;
	t26 = cos(qJ(3));
	t29 = t22 * t26;
	t24 = cos(pkin(8));
	t28 = t24 * t25;
	t27 = t24 * t26;
	t23 = cos(pkin(9));
	t21 = sin(pkin(9));
	t1 = [0, 0, -t23 * t28 + t29, 0, 0; 0, 0, -t23 * t30 - t27, 0, 0; 0, 0, -t21 * t25, 0, 0; 0, 0, -t23 * t27 - t30, 0, 0; 0, 0, -t23 * t29 + t28, 0, 0; 0, 0, -t21 * t26, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:18
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (19->13), mult. (57->29), div. (0->0), fcn. (88->8), ass. (0->20)
	t68 = sin(pkin(9));
	t72 = sin(qJ(4));
	t82 = t68 * t72;
	t74 = cos(qJ(4));
	t81 = t68 * t74;
	t75 = cos(qJ(3));
	t80 = t68 * t75;
	t69 = sin(pkin(8));
	t73 = sin(qJ(3));
	t79 = t69 * t73;
	t78 = t69 * t75;
	t71 = cos(pkin(8));
	t77 = t71 * t73;
	t76 = t71 * t75;
	t70 = cos(pkin(9));
	t67 = t70 * t76 + t79;
	t66 = -t70 * t77 + t78;
	t65 = t70 * t78 - t77;
	t64 = -t70 * t79 - t76;
	t1 = [0, 0, t66 * t74, -t67 * t72 + t71 * t81, 0; 0, 0, t64 * t74, -t65 * t72 + t69 * t81, 0; 0, 0, -t73 * t81, -t70 * t74 - t72 * t80, 0; 0, 0, -t66 * t72, -t67 * t74 - t71 * t82, 0; 0, 0, -t64 * t72, -t65 * t74 - t69 * t82, 0; 0, 0, t73 * t82, t70 * t72 - t74 * t80, 0; 0, 0, t67, 0, 0; 0, 0, t65, 0, 0; 0, 0, t80, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:18
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->14), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->28)
	t97 = sin(pkin(9));
	t98 = sin(pkin(8));
	t110 = t97 * t98;
	t100 = cos(pkin(8));
	t109 = t100 * t97;
	t101 = sin(qJ(3));
	t108 = t101 * t97;
	t102 = cos(qJ(3));
	t107 = t97 * t102;
	t106 = t98 * t101;
	t105 = t98 * t102;
	t104 = t100 * t101;
	t103 = t100 * t102;
	t99 = cos(pkin(9));
	t96 = qJ(4) + qJ(5);
	t95 = cos(t96);
	t94 = sin(t96);
	t93 = t99 * t103 + t106;
	t92 = -t99 * t104 + t105;
	t91 = t99 * t105 - t104;
	t90 = -t99 * t106 - t103;
	t89 = -t95 * t107 + t99 * t94;
	t88 = -t94 * t107 - t99 * t95;
	t87 = -t94 * t109 - t93 * t95;
	t86 = t95 * t109 - t93 * t94;
	t85 = -t94 * t110 - t91 * t95;
	t84 = t95 * t110 - t91 * t94;
	t1 = [0, 0, t92 * t95, t86, t86; 0, 0, t90 * t95, t84, t84; 0, 0, -t95 * t108, t88, t88; 0, 0, -t92 * t94, t87, t87; 0, 0, -t90 * t94, t85, t85; 0, 0, t94 * t108, t89, t89; 0, 0, t93, 0, 0; 0, 0, t91, 0, 0; 0, 0, t107, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end