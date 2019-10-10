% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
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
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(10);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(10);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0, 0; t16, 0, -t19, 0, 0, 0; 0, 0, t15, 0, 0, 0; t19, 0, -t16, 0, 0, 0; -t17, 0, -t18, 0, 0, 0; 0, 0, -t14, 0, 0, 0; t12, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t46 = qJ(1) + pkin(10);
	t44 = sin(t46);
	t47 = sin(qJ(3));
	t50 = t44 * t47;
	t45 = cos(t46);
	t48 = cos(qJ(3));
	t49 = t45 * t48;
	t43 = t45 * t47;
	t42 = t44 * t48;
	t1 = [t45, 0, 0, 0, 0, 0; t44, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t42, 0, t43, 0, 0, 0; -t49, 0, t50, 0, 0, 0; 0, 0, -t48, 0, 0, 0; -t50, 0, t49, 0, 0, 0; t43, 0, t42, 0, 0, 0; 0, 0, t47, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t67 = sin(qJ(5));
	t68 = sin(qJ(3));
	t74 = t68 * t67;
	t69 = cos(qJ(5));
	t73 = t68 * t69;
	t70 = cos(qJ(3));
	t72 = t70 * t67;
	t71 = t70 * t69;
	t66 = qJ(1) + pkin(10);
	t65 = cos(t66);
	t64 = sin(t66);
	t63 = -t64 * t74 + t65 * t69;
	t62 = t64 * t73 + t65 * t67;
	t61 = t64 * t69 + t65 * t74;
	t60 = -t64 * t67 + t65 * t73;
	t1 = [t63, 0, t65 * t72, 0, t60, 0; t61, 0, t64 * t72, 0, t62, 0; 0, 0, t74, 0, -t71, 0; -t62, 0, t65 * t71, 0, -t61, 0; t60, 0, t64 * t71, 0, t63, 0; 0, 0, t73, 0, t72, 0; -t64 * t70, 0, -t65 * t68, 0, 0, 0; t65 * t70, 0, -t64 * t68, 0, 0, 0; 0, 0, t70, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (83->14), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
	t91 = qJ(5) + qJ(6);
	t88 = sin(t91);
	t92 = sin(qJ(3));
	t96 = t92 * t88;
	t89 = cos(t91);
	t95 = t92 * t89;
	t93 = cos(qJ(3));
	t85 = t93 * t88;
	t94 = t93 * t89;
	t90 = qJ(1) + pkin(10);
	t87 = cos(t90);
	t86 = sin(t90);
	t84 = -t86 * t96 + t87 * t89;
	t83 = t86 * t95 + t87 * t88;
	t82 = t86 * t89 + t87 * t96;
	t81 = -t86 * t88 + t87 * t95;
	t1 = [t84, 0, t87 * t85, 0, t81, t81; t82, 0, t86 * t85, 0, t83, t83; 0, 0, t96, 0, -t94, -t94; -t83, 0, t87 * t94, 0, -t82, -t82; t81, 0, t86 * t94, 0, t84, t84; 0, 0, t95, 0, t85, t85; -t86 * t93, 0, -t87 * t92, 0, 0, 0; t87 * t93, 0, -t86 * t92, 0, 0, 0; 0, 0, t93, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end