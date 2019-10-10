% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
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
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
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
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(9);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t23 = t18 * t15;
	t16 = cos(t17);
	t22 = t18 * t16;
	t19 = cos(qJ(1));
	t21 = t19 * t15;
	t20 = t19 * t16;
	t1 = [-t22, -t21, 0, 0, 0, 0; t20, -t23, 0, 0, 0, 0; 0, t16, 0, 0, 0, 0; t23, -t20, 0, 0, 0, 0; -t21, -t22, 0, 0, 0, 0; 0, -t15, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t50 = qJ(2) + pkin(9);
	t48 = sin(t50);
	t51 = sin(qJ(1));
	t54 = t51 * t48;
	t49 = cos(t50);
	t52 = cos(qJ(1));
	t53 = t52 * t49;
	t47 = t52 * t48;
	t46 = t51 * t49;
	t1 = [t52, 0, 0, 0, 0, 0; t51, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t46, t47, 0, 0, 0, 0; -t53, t54, 0, 0, 0, 0; 0, -t49, 0, 0, 0, 0; -t54, t53, 0, 0, 0, 0; t47, t46, 0, 0, 0, 0; 0, t48, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t75 = sin(qJ(5));
	t76 = sin(qJ(1));
	t82 = t76 * t75;
	t77 = cos(qJ(5));
	t81 = t76 * t77;
	t78 = cos(qJ(1));
	t80 = t78 * t75;
	t79 = t78 * t77;
	t74 = qJ(2) + pkin(9);
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = -t72 * t82 + t79;
	t70 = t72 * t81 + t80;
	t69 = t72 * t80 + t81;
	t68 = t72 * t79 - t82;
	t1 = [t71, t73 * t80, 0, 0, t68, 0; t69, t73 * t82, 0, 0, t70, 0; 0, t72 * t75, 0, 0, -t73 * t77, 0; -t70, t73 * t79, 0, 0, -t69, 0; t68, t73 * t81, 0, 0, t71, 0; 0, t72 * t77, 0, 0, t73 * t75, 0; -t76 * t73, -t78 * t72, 0, 0, 0, 0; t78 * t73, -t76 * t72, 0, 0, 0, 0; 0, t73, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t79 = sin(qJ(5));
	t80 = sin(qJ(1));
	t86 = t80 * t79;
	t81 = cos(qJ(5));
	t85 = t80 * t81;
	t82 = cos(qJ(1));
	t84 = t82 * t79;
	t83 = t82 * t81;
	t78 = qJ(2) + pkin(9);
	t77 = cos(t78);
	t76 = sin(t78);
	t75 = -t76 * t86 + t83;
	t74 = t76 * t85 + t84;
	t73 = t76 * t84 + t85;
	t72 = t76 * t83 - t86;
	t1 = [t75, t77 * t84, 0, 0, t72, 0; t73, t77 * t86, 0, 0, t74, 0; 0, t76 * t79, 0, 0, -t77 * t81, 0; -t74, t77 * t83, 0, 0, -t73, 0; t72, t77 * t85, 0, 0, t75, 0; 0, t76 * t81, 0, 0, t77 * t79, 0; -t80 * t77, -t82 * t76, 0, 0, 0, 0; t82 * t77, -t80 * t76, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end