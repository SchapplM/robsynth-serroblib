% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
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
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(9);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t12 = cos(pkin(10));
	t11 = sin(pkin(10));
	t10 = qJ(1) + pkin(9);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [-t8 * t12, 0, 0, 0, 0, 0; t9 * t12, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t8 * t11, 0, 0, 0, 0, 0; -t9 * t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t18 = pkin(10) + qJ(4);
	t14 = sin(t18);
	t19 = qJ(1) + pkin(9);
	t15 = sin(t19);
	t23 = t15 * t14;
	t16 = cos(t18);
	t22 = t15 * t16;
	t17 = cos(t19);
	t21 = t17 * t14;
	t20 = t17 * t16;
	t1 = [-t22, 0, 0, -t21, 0, 0; t20, 0, 0, -t23, 0, 0; 0, 0, 0, t16, 0, 0; t23, 0, 0, -t20, 0, 0; -t21, 0, 0, -t22, 0, 0; 0, 0, 0, -t14, 0, 0; t17, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t76 = qJ(1) + pkin(9);
	t72 = sin(t76);
	t77 = sin(qJ(5));
	t84 = t72 * t77;
	t78 = cos(qJ(5));
	t83 = t72 * t78;
	t75 = pkin(10) + qJ(4);
	t73 = cos(t75);
	t82 = t73 * t77;
	t81 = t73 * t78;
	t74 = cos(t76);
	t80 = t74 * t77;
	t79 = t74 * t78;
	t71 = sin(t75);
	t70 = t73 * t79 + t84;
	t69 = -t73 * t80 + t83;
	t68 = -t72 * t81 + t80;
	t67 = t72 * t82 + t79;
	t1 = [t68, 0, 0, -t71 * t79, t69, 0; t70, 0, 0, -t71 * t83, -t67, 0; 0, 0, 0, t81, -t71 * t77, 0; t67, 0, 0, t71 * t80, -t70, 0; t69, 0, 0, t71 * t84, t68, 0; 0, 0, 0, -t82, -t71 * t78, 0; -t72 * t71, 0, 0, t74 * t73, 0, 0; t74 * t71, 0, 0, t72 * t73, 0, 0; 0, 0, 0, t71, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t90 = qJ(1) + pkin(9);
	t86 = sin(t90);
	t91 = sin(qJ(5));
	t98 = t86 * t91;
	t92 = cos(qJ(5));
	t97 = t86 * t92;
	t89 = pkin(10) + qJ(4);
	t87 = cos(t89);
	t96 = t87 * t91;
	t95 = t87 * t92;
	t88 = cos(t90);
	t94 = t88 * t91;
	t93 = t88 * t92;
	t85 = sin(t89);
	t84 = t87 * t93 + t98;
	t83 = t87 * t94 - t97;
	t82 = t86 * t95 - t94;
	t81 = -t86 * t96 - t93;
	t1 = [-t82, 0, 0, -t85 * t93, -t83, 0; t84, 0, 0, -t85 * t97, t81, 0; 0, 0, 0, t95, -t85 * t91, 0; -t86 * t85, 0, 0, t88 * t87, 0, 0; t88 * t85, 0, 0, t86 * t87, 0, 0; 0, 0, 0, t85, 0, 0; t81, 0, 0, -t85 * t94, t84, 0; t83, 0, 0, -t85 * t98, t82, 0; 0, 0, 0, t96, t85 * t92, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end