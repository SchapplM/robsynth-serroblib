% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
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
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
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
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
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
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t47 = pkin(9) + qJ(3);
	t45 = sin(t47);
	t48 = sin(qJ(1));
	t51 = t48 * t45;
	t46 = cos(t47);
	t49 = cos(qJ(1));
	t50 = t49 * t46;
	t44 = t49 * t45;
	t43 = t48 * t46;
	t1 = [t49, 0, 0, 0, 0, 0; t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t43, 0, t44, 0, 0, 0; -t50, 0, t51, 0, 0, 0; 0, 0, -t46, 0, 0, 0; -t51, 0, t50, 0, 0, 0; t44, 0, t43, 0, 0, 0; 0, 0, t45, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t72 = sin(qJ(5));
	t73 = sin(qJ(1));
	t79 = t73 * t72;
	t74 = cos(qJ(5));
	t78 = t73 * t74;
	t75 = cos(qJ(1));
	t77 = t75 * t72;
	t76 = t75 * t74;
	t71 = pkin(9) + qJ(3);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = -t69 * t79 + t76;
	t67 = t69 * t78 + t77;
	t66 = t69 * t77 + t78;
	t65 = t69 * t76 - t79;
	t1 = [t68, 0, t70 * t77, 0, t65, 0; t66, 0, t70 * t79, 0, t67, 0; 0, 0, t69 * t72, 0, -t70 * t74, 0; -t67, 0, t70 * t76, 0, -t66, 0; t65, 0, t70 * t78, 0, t68, 0; 0, 0, t69 * t74, 0, t70 * t72, 0; -t73 * t70, 0, -t75 * t69, 0, 0, 0; t75 * t70, 0, -t73 * t69, 0, 0, 0; 0, 0, t70, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t76 = sin(qJ(5));
	t77 = sin(qJ(1));
	t83 = t77 * t76;
	t78 = cos(qJ(5));
	t82 = t77 * t78;
	t79 = cos(qJ(1));
	t81 = t79 * t76;
	t80 = t79 * t78;
	t75 = pkin(9) + qJ(3);
	t74 = cos(t75);
	t73 = sin(t75);
	t72 = -t73 * t83 + t80;
	t71 = t73 * t82 + t81;
	t70 = t73 * t81 + t82;
	t69 = t73 * t80 - t83;
	t1 = [t72, 0, t74 * t81, 0, t69, 0; t70, 0, t74 * t83, 0, t71, 0; 0, 0, t73 * t76, 0, -t74 * t78, 0; -t71, 0, t74 * t80, 0, -t70, 0; t69, 0, t74 * t82, 0, t72, 0; 0, 0, t73 * t78, 0, t74 * t76, 0; -t77 * t74, 0, -t79 * t73, 0, 0, 0; t79 * t74, 0, -t77 * t73, 0, 0, 0; 0, 0, t74, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end