% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
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
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(9));
	t7 = sin(pkin(9));
	t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(9) + qJ(3);
	t12 = sin(t14);
	t15 = sin(qJ(1));
	t20 = t15 * t12;
	t13 = cos(t14);
	t19 = t15 * t13;
	t16 = cos(qJ(1));
	t18 = t16 * t12;
	t17 = t16 * t13;
	t1 = [-t19, 0, -t18, 0, 0, 0; t17, 0, -t20, 0, 0, 0; 0, 0, t13, 0, 0, 0; t20, 0, -t17, 0, 0, 0; -t18, 0, -t19, 0, 0, 0; 0, 0, -t12, 0, 0, 0; t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t71 = sin(qJ(4));
	t72 = sin(qJ(1));
	t78 = t72 * t71;
	t73 = cos(qJ(4));
	t77 = t72 * t73;
	t74 = cos(qJ(1));
	t76 = t74 * t71;
	t75 = t74 * t73;
	t70 = pkin(9) + qJ(3);
	t69 = cos(t70);
	t68 = sin(t70);
	t67 = t69 * t75 + t78;
	t66 = -t69 * t76 + t77;
	t65 = -t69 * t77 + t76;
	t64 = t69 * t78 + t75;
	t1 = [t65, 0, -t68 * t75, t66, 0, 0; t67, 0, -t68 * t77, -t64, 0, 0; 0, 0, t69 * t73, -t68 * t71, 0, 0; t64, 0, t68 * t76, -t67, 0, 0; t66, 0, t68 * t78, t65, 0, 0; 0, 0, -t69 * t71, -t68 * t73, 0, 0; -t72 * t68, 0, t74 * t69, 0, 0, 0; t74 * t68, 0, t72 * t69, 0, 0, 0; 0, 0, t68, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:56
	% EndTime: 2019-10-10 01:14:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t81 = qJ(4) + pkin(10);
	t77 = sin(t81);
	t82 = sin(qJ(1));
	t89 = t82 * t77;
	t80 = pkin(9) + qJ(3);
	t78 = cos(t80);
	t88 = t82 * t78;
	t79 = cos(t81);
	t87 = t82 * t79;
	t83 = cos(qJ(1));
	t86 = t83 * t77;
	t85 = t83 * t78;
	t84 = t83 * t79;
	t76 = sin(t80);
	t75 = t78 * t84 + t89;
	t74 = -t77 * t85 + t87;
	t73 = -t78 * t87 + t86;
	t72 = t77 * t88 + t84;
	t1 = [t73, 0, -t76 * t84, t74, 0, 0; t75, 0, -t76 * t87, -t72, 0, 0; 0, 0, t78 * t79, -t76 * t77, 0, 0; t72, 0, t76 * t86, -t75, 0, 0; t74, 0, t76 * t89, t73, 0, 0; 0, 0, -t78 * t77, -t76 * t79, 0, 0; -t82 * t76, 0, t85, 0, 0, 0; t83 * t76, 0, t88, 0, 0, 0; 0, 0, t76, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:56
	% EndTime: 2019-10-10 01:14:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t96 = qJ(4) + pkin(10);
	t92 = sin(t96);
	t97 = sin(qJ(1));
	t104 = t97 * t92;
	t95 = pkin(9) + qJ(3);
	t93 = cos(t95);
	t103 = t97 * t93;
	t94 = cos(t96);
	t102 = t97 * t94;
	t98 = cos(qJ(1));
	t101 = t98 * t92;
	t100 = t98 * t93;
	t99 = t98 * t94;
	t91 = sin(t95);
	t90 = t93 * t99 + t104;
	t89 = t92 * t100 - t102;
	t88 = t93 * t102 - t101;
	t87 = -t92 * t103 - t99;
	t1 = [-t88, 0, -t91 * t99, -t89, 0, 0; t90, 0, -t91 * t102, t87, 0, 0; 0, 0, t93 * t94, -t91 * t92, 0, 0; -t97 * t91, 0, t100, 0, 0, 0; t98 * t91, 0, t103, 0, 0, 0; 0, 0, t91, 0, 0, 0; t87, 0, -t91 * t101, t90, 0, 0; t89, 0, -t91 * t104, t88, 0, 0; 0, 0, t93 * t92, t91 * t94, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end